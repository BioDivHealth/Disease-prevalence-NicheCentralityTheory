## Comprehensive Niche Model Suite
## ------------------------------------------------------------
## Systematic testing of all single-effect and targeted pairwise 
## combination models with gold-standard base structure:
## - Beta-Binomial likelihood
## - Spatial GP (gp(lon_z, lat_z, k=20))
## - Species RE ((1 | host_species))
##
## Run from the project root:
##   Rscript "code/2.4_comprehensive_models.R"
##

# 1. Packages -------------------------------------------------------------
pacman::p_load(tidyverse, brms, loo, here)

# Spatial GP fits are most stable with cores = 1 (avoids SIGPIPE/socket issues).
options(mc.cores = 1)

# Control flags
fit_all_screening <- FALSE  # Set to TRUE to fit the entire ~40 model suite (screening)
run_loo <- TRUE            # Compute LOO for all loaded models
production_models <- c("m1centrality_within_spatial",           
                       "m2suit_marg_raw_additive_spatial",      
                       "m2suit_marg_within_additive_spatial",   
                       "m2suit_marg_raw_interaction_spatial",   
                       "m2suit_marg_within_interaction_spatial")   # Models to fit with high-precision settings (listed at end)
# Helpers to satisfy linters when using pacman::p_load() in scripts
bf <- brms::bf
brm <- brms::brm
beta_binomial <- brms::beta_binomial
save_pars <- brms::save_pars
here <- here::here

# 2. Load data and prepare -----------------------------------------------
dat <- readRDS(here("Data", "dat_clean_agg.rds"))

scaled_covars <- dat |>
  mutate(
    assay_group = factor(assay_group),
    host_family = factor(host_family),
    pathogen_family = factor(pathogen_family),
    host_species = factor(host_species),
    pathogen_species_cleaned = factor(pathogen_species_cleaned),
    obs_id = row_number(),
    lon_z = as.numeric(scale(longitude)),
    lat_z = as.numeric(scale(latitude))
  )

stopifnot(all(is.finite(scaled_covars$lon_z)), all(is.finite(scaled_covars$lat_z)))

# 2.1 Extend within/between decomposition to all niche variables ---------
# Extend from 2.3 to include all 6 niche variables
vars_for_within_between <- c(
  "prob_occur_orderNorm",
  "Suitability_orderNorm",
  "Marginality_orderNorm",
  "Specificity_orderNorm",
  "Centroid_d_orderNorm",
  "Boundary_d_orderNorm"
)

missing_vars <- setdiff(vars_for_within_between, names(scaled_covars))
if (length(missing_vars) > 0) {
  stop("Missing required predictors in dat_clean_agg.rds: ", paste(missing_vars, collapse = ", "))
}

# Calculate species-level means (between-species component)
species_means_raw <- scaled_covars |>
  group_by(host_species) |>
  summarise(
    across(all_of(vars_for_within_between), ~ mean(.x, na.rm = TRUE), .names = "{.col}_between_raw"),
    .groups = "drop"
  )

# Scale the species means
species_means_scaled <- species_means_raw |>
  mutate(
    across(
      ends_with("_between_raw"),
      ~ as.numeric(scale(.x)),
      .names = "{.col}_z"
    )
  )

# Join back and calculate within-species deviations
scaled_covars <- scaled_covars |>
  left_join(species_means_scaled, by = "host_species") |>
  mutate(
    # Within-species deviations (raw) then z-score across rows for stable priors
    prob_occur_within_raw = prob_occur_orderNorm - prob_occur_orderNorm_between_raw,
    Suitability_within_raw = Suitability_orderNorm - Suitability_orderNorm_between_raw,
    Marginality_within_raw = Marginality_orderNorm - Marginality_orderNorm_between_raw,
    Specificity_within_raw = Specificity_orderNorm - Specificity_orderNorm_between_raw,
    Centroid_d_within_raw = Centroid_d_orderNorm - Centroid_d_orderNorm_between_raw,
    Boundary_d_within_raw = Boundary_d_orderNorm - Boundary_d_orderNorm_between_raw,
    # Re-standardize for stable priors (SD = 1)
    prob_occur_within_z = as.numeric(scale(prob_occur_within_raw)),
    Suitability_within_z = as.numeric(scale(Suitability_within_raw)),
    Marginality_within_z = as.numeric(scale(Marginality_within_raw)),
    Specificity_within_z = as.numeric(scale(Specificity_within_raw)),
    Centroid_d_within_z = as.numeric(scale(Centroid_d_within_raw)),
    Boundary_d_within_z = as.numeric(scale(Boundary_d_within_raw))
  )

# Centrality score (composite: Suitability - Marginality)
scaled_covars <- scaled_covars |>
  mutate(
    centrality_score_raw = Suitability_orderNorm - Marginality_orderNorm,
    centrality_between_raw = Suitability_orderNorm_between_raw - Marginality_orderNorm_between_raw,
    centrality_within_raw = centrality_score_raw - centrality_between_raw,
    centrality_within_z = as.numeric(scale(centrality_within_raw)),
    centrality_between_z = as.numeric(scale(centrality_between_raw))
  )

# Verify all variables are finite
stopifnot(
  all(is.finite(scaled_covars$prob_occur_within_z)),
  all(is.finite(scaled_covars$Suitability_within_z)),
  all(is.finite(scaled_covars$Marginality_within_z)),
  all(is.finite(scaled_covars$Specificity_within_z)),
  all(is.finite(scaled_covars$Centroid_d_within_z)),
  all(is.finite(scaled_covars$Boundary_d_within_z)),
  all(is.finite(scaled_covars$centrality_within_z)),
  all(is.finite(scaled_covars$centrality_between_z))
)

# 3. Priors ---------------------------------------------------------------
priors_beta_binom <- prior(normal(0, 1.5), class = "Intercept") +
  prior(normal(0, 1), class = "b") +
  prior(exponential(1), class = "sd") +
  prior(exponential(1), class = "phi")

priors_beta_binom_spatial <- priors_beta_binom +
  prior(exponential(1), class = "sdgp") +
  prior(exponential(1), class = "lscale")

# 4. Helper functions ----------------------------------------------------
load_if_exists <- function(path_no_ext) {
  rds_path <- paste0(path_no_ext, ".rds")
  if (file.exists(rds_path)) return(readRDS(rds_path))
  NULL
}

# Generate single-effect formula
# Note: brms::bf() accepts formula objects; as.formula() parses the string
make_single_formula <- function(predictor) {
  formula_str <- paste0(
    "number_positive | trials(number_tested) ~ 1 + ", predictor, 
    " + assay_group + (1 | host_species) + gp(lon_z, lat_z, k = 20)"
  )
  bf(as.formula(formula_str))
}

# Generate pairwise formula (additive or interaction)
make_pair_formula <- function(var1, var2, interaction = FALSE) {
  if (interaction) {
    formula_str <- paste0(
      "number_positive | trials(number_tested) ~ 1 + ", var1, " * ", var2,
      " + assay_group + (1 | host_species) + gp(lon_z, lat_z, k = 20)"
    )
  } else {
    formula_str <- paste0(
      "number_positive | trials(number_tested) ~ 1 + ", var1, " + ", var2,
      " + assay_group + (1 | host_species) + gp(lon_z, lat_z, k = 20)"
    )
  }
  bf(as.formula(formula_str))
}

# Model fitting wrapper with screening settings
fit_model_screen <- function(formula, model_name, data = scaled_covars) {
  model_path <- here("Results", "models_v4", model_name)
  
  # Check if model already exists
  existing <- load_if_exists(model_path)
  if (!is.null(existing)) {
    message("Loading existing model: ", model_name)
    return(existing)
  }
  
  message("Fitting model: ", model_name)
  brm(
    formula = formula,
    family = beta_binomial(link = "logit"),
    data = data,
    prior = priors_beta_binom_spatial,
    iter = 1000,
    warmup = 500,
    chains = 2,
    cores = 1,  # Prevents SIGPIPE crashes with GP
    init = 0,
    init_r = 0.1,
    control = list(adapt_delta = 0.97, max_treedepth = 12),
    seed = 123,
    file = model_path
  )
}

# Model fitting wrapper with production settings (high precision)
fit_model_production <- function(formula, model_name, data = scaled_covars) {
  # Production models are saved with a _final suffix
  final_name <- paste0(model_name, "_final")
  model_path <- here("Results", "models_v4", final_name)
  
  # Check if model already exists
  existing <- load_if_exists(model_path)
  if (!is.null(existing)) {
    message("Loading existing production model: ", final_name)
    return(existing)
  }
  
  message("Fitting production model: ", final_name)
  brm(
    formula = formula,
    family = beta_binomial(link = "logit"),
    data = data,
    prior = priors_beta_binom_spatial,
    iter = 2000,
    warmup = 1000,
    chains = 4,
    cores = 1,  # Keep cores=1 for stability with GP
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    seed = 123,
    file = model_path
  )
}

# Wrapper that chooses between screen and production
fit_model_auto <- function(formula, model_name, data = scaled_covars) {
  if (model_name %in% production_models) {
    return(fit_model_production(formula, model_name, data))
  } else {
    return(fit_model_screen(formula, model_name, data))
  }
}

# Helper to load a model, checking for _final version first
load_model_auto <- function(model_name) {
  # Try loading the final version first
  final_path <- here("Results", "models_v4", paste0(model_name, "_final"))
  mod <- load_if_exists(final_path)
  
  if (is.null(mod)) {
    # Fall back to the screen version
    screen_path <- here("Results", "models_v4", model_name)
    mod <- load_if_exists(screen_path)
  }
  
  return(mod)
}

# 5. Define model specifications ------------------------------------------
# Registry to store all potential model formulas
model_formulas <- list()

# 5.0 Spatial baseline
m0_name <- "m0_spatial"
f_m0_spatial <- bf(
  number_positive | trials(number_tested) ~ 1 + assay_group + (1 | host_species) +
    gp(lon_z, lat_z, k = 20)
)
model_formulas[[m0_name]] <- f_m0_spatial

# 5.1 Single-effect models (raw _orderNorm)
single_raw_vars <- c(
  "prob_occur_orderNorm",
  "Suitability_orderNorm",
  "Marginality_orderNorm",
  "Specificity_orderNorm",
  "Centroid_d_orderNorm",
  "Boundary_d_orderNorm"
)
single_raw_names <- c(
  "m1prob_raw_spatial", "m1suit_raw_spatial", "m1marg_raw_spatial",
  "m1spec_raw_spatial", "m1centroid_raw_spatial", "m1boundary_raw_spatial"
)
for (i in seq_along(single_raw_vars)) {
  model_formulas[[single_raw_names[i]]] <- make_single_formula(single_raw_vars[i])
}

# 5.2 Single-effect models (within-species)
single_within_vars <- c(
  "prob_occur_within_z", "Suitability_within_z", "Marginality_within_z",
  "Specificity_within_z", "Centroid_d_within_z", "Boundary_d_within_z"
)
single_within_names <- c(
  "m1prob_within_spatial", "m1suit_within_spatial", "m1marg_within_spatial",
  "m1spec_within_spatial", "m1centroid_within_spatial", "m1boundary_within_spatial"
)
for (i in seq_along(single_within_vars)) {
  model_formulas[[single_within_names[i]]] <- make_single_formula(single_within_vars[i])
}

# 5.3 Centrality models
model_formulas[["m1centrality_within_spatial"]] <- make_single_formula("centrality_within_z")
model_formulas[["m1centrality_between_spatial"]] <- make_single_formula("centrality_between_z")

# 5.4 Targeted pairwise combinations
targeted_pairs <- list(
  # Centrality-focused
  list("Suitability_orderNorm", "Marginality_orderNorm", "suit_marg"),
  list("Suitability_orderNorm", "Specificity_orderNorm", "suit_spec"),
  list("Marginality_orderNorm", "Specificity_orderNorm", "marg_spec"),
  # Geographic + Niche
  list("Centroid_d_orderNorm", "Suitability_orderNorm", "centroid_suit"),
  list("Boundary_d_orderNorm", "Marginality_orderNorm", "boundary_marg"),
  list("Centroid_d_orderNorm", "Boundary_d_orderNorm", "centroid_boundary"),
  # SDM + Niche (avoiding prob_occur × Suitability)
  list("prob_occur_orderNorm", "Marginality_orderNorm", "prob_marg"),
  list("prob_occur_orderNorm", "Specificity_orderNorm", "prob_spec")
)

for (pair in targeted_pairs) {
  name_base <- pair[[3]]
  # Raw
  model_formulas[[paste0("m2", name_base, "_raw_additive_spatial")]] <- 
    make_pair_formula(pair[[1]], pair[[2]], interaction = FALSE)
  model_formulas[[paste0("m2", name_base, "_raw_interaction_spatial")]] <- 
    make_pair_formula(pair[[1]], pair[[2]], interaction = TRUE)
}

# Within-species pairs
targeted_pairs_within <- list(
  list("Suitability_within_z", "Marginality_within_z", "suit_marg"),
  list("Suitability_within_z", "Specificity_within_z", "suit_spec"),
  list("Marginality_within_z", "Specificity_within_z", "marg_spec"),
  list("Centroid_d_within_z", "Suitability_within_z", "centroid_suit"),
  list("Boundary_d_within_z", "Marginality_within_z", "boundary_marg"),
  list("Centroid_d_within_z", "Boundary_d_within_z", "centroid_boundary"),
  list("prob_occur_within_z", "Marginality_within_z", "prob_marg"),
  list("prob_occur_within_z", "Specificity_within_z", "prob_spec")
)

for (pair in targeted_pairs_within) {
  name_base <- pair[[3]]
  model_formulas[[paste0("m2", name_base, "_within_additive_spatial")]] <- 
    make_pair_formula(pair[[1]], pair[[2]], interaction = FALSE)
  model_formulas[[paste0("m2", name_base, "_within_interaction_spatial")]] <- 
    make_pair_formula(pair[[1]], pair[[2]], interaction = TRUE)
}

# Create models directory
models_dir <- here("Results", "models_v4")
if (!dir.exists(models_dir)) dir.create(models_dir, recursive = TRUE)

# 6. Load/Fit models -----------------------------------------------------------
all_models <- list()

# 6.0 Spatial baseline
message("\n========== Loading spatial baseline ==========")
if (fit_all_screening) {
    all_models[[m0_name]] <- fit_model_auto(model_formulas[[m0_name]], m0_name)
} else {
  all_models[[m0_name]] <- load_model_auto(m0_name)
}

# 6.1-6.7 Bulk Loading/Screening Loop
message("\n========== Loading/Screening all spatial models ==========")
for (m_name in names(model_formulas)) {
  if (m_name == m0_name) next # already handled
  
  if (fit_all_screening) {
    all_models[[m_name]] <- fit_model_auto(model_formulas[[m_name]], m_name)
  } else {
    mod <- load_model_auto(m_name)
    if (!is.null(mod)) all_models[[m_name]] <- mod
  }
}

# Remove NULL entries (models that don't exist yet)
all_models <- all_models[!sapply(all_models, is.null)]

message("\n========== Model fitting complete ==========")
message("Total models loaded/fitted: ", length(all_models))

# 7. Load non-spatial models from 2.2 for comprehensive comparison --------
message("\n========== Loading non-spatial models from 2.2_models.R ==========")

non_spatial_models <- list()

# Load baselines from models_v1
non_spatial_models$m0_base_bb <- load_if_exists(here("Results", "models_v1", "m0_base_bb"))
non_spatial_models$m0_beta_binom <- load_if_exists(here("Results", "models_v1", "m0_beta_binom"))
non_spatial_models$m0_controls <- load_if_exists(here("Results", "models_v1", "m0_controls"))
non_spatial_models$m0_base <- load_if_exists(here("Results", "models_v1", "m0_base"))

# Load M1 models from models_v1 (with controls)
non_spatial_models$m1prob <- load_if_exists(here("Results", "models_v1", "m1prob"))
non_spatial_models$m1suit <- load_if_exists(here("Results", "models_v1", "m1suit"))
non_spatial_models$m1marg <- load_if_exists(here("Results", "models_v1", "m1marg"))
non_spatial_models$m1spec <- load_if_exists(here("Results", "models_v1", "m1spec"))
non_spatial_models$m1centr <- load_if_exists(here("Results", "models_v1", "m1centr"))
non_spatial_models$m1bound <- load_if_exists(here("Results", "models_v1", "m1bound"))

# Load M1 models from models_v2 (minimal baseline)
non_spatial_models$m1prob_v2 <- load_if_exists(here("Results", "models_v2", "m1prob_v2"))
non_spatial_models$m1suit_v2 <- load_if_exists(here("Results", "models_v2", "m1suit_v2"))
non_spatial_models$m1marg_v2 <- load_if_exists(here("Results", "models_v2", "m1marg_v2"))
non_spatial_models$m1spec_v2 <- load_if_exists(here("Results", "models_v2", "m1spec_v2"))
non_spatial_models$m1centr_v2 <- load_if_exists(here("Results", "models_v2", "m1centr_v2"))
non_spatial_models$m1bound_v2 <- load_if_exists(here("Results", "models_v2", "m1bound_v2"))

# Load M2 models from models_v2
non_spatial_models$m2a_v2 <- load_if_exists(here("Results", "models_v2", "m2a_v2"))
non_spatial_models$m2b_v2 <- load_if_exists(here("Results", "models_v2", "m2b_v2"))
non_spatial_models$m2c_v2 <- load_if_exists(here("Results", "models_v2", "m2c_v2"))
non_spatial_models$m2d_v2 <- load_if_exists(here("Results", "models_v2", "m2d_v2"))

# Remove NULL entries
non_spatial_models <- non_spatial_models[!sapply(non_spatial_models, is.null)]

message("Loaded ", length(non_spatial_models), " non-spatial models from 2.2_models.R")

# 8. Model comparison with LOO -------------------------------------------
available_loos <- list()
if (isTRUE(run_loo)) {
  message("\n========== Computing LOO for all models ==========")
  
  # Compute LOO for spatial models (from 2.4)
  for (model_name in names(all_models)) {
    if (!is.null(all_models[[model_name]])) {
      message("Computing LOO for spatial model: ", model_name)
      tryCatch({
        available_loos[[model_name]] <- loo(all_models[[model_name]])
      }, error = function(e) {
        message("  Error computing LOO for ", model_name, ": ", conditionMessage(e))
      })
    }
  }
  
  # Compute LOO for non-spatial models (from 2.2)
  for (model_name in names(non_spatial_models)) {
    if (!is.null(non_spatial_models[[model_name]])) {
      message("Computing LOO for non-spatial model: ", model_name)
      tryCatch({
        available_loos[[model_name]] <- loo(non_spatial_models[[model_name]])
      }, error = function(e) {
        message("  Error computing LOO for ", model_name, ": ", conditionMessage(e))
      })
    }
  }

  message("\n========== COMPREHENSIVE MODEL COMPARISON ==========")
  message("Total models with LOO: ", length(available_loos))

  if (length(available_loos) >= 2) {
    comp <- loo_compare(available_loos)
    print(comp)

    # Print top 10 models
    message("\n========== TOP 10 MODELS ==========")
    n_top <- min(10, nrow(comp))
    print(comp[1:n_top, ])

    # Save comparison to file
    write.csv(comp, here("Results", "models_v4", "loo_comparison.csv"))
    message("\nComparison saved to: Results/models_v4/loo_comparison.csv")
  } else {
    message("Need at least 2 models with LOO for comparison")
  }
} else {
  message("\nSkipping LOO computation (run_loo = FALSE).")
}

if (isTRUE(run_loo) && exists("comp")) {
  write.csv(round(comp, 3), here("Results","models_v4","loo_comparison_rounded.csv"))
}

# 9. Summary of models fitted ---------------------------------------------
# Combine spatial and non-spatial models for summary
all_models_combined <- c(all_models, non_spatial_models)

model_summary <- tibble(
  model_name = names(all_models_combined),
  model_type = if_else(names(all_models_combined) %in% names(all_models), "spatial", "non_spatial"),
  is_production = sapply(names(all_models_combined), function(nm) {
    file.exists(here("Results", "models_v4", paste0(nm, "_final.rds")))
  }),
  fitted = !sapply(all_models_combined, is.null),
  has_loo = names(all_models_combined) %in% names(available_loos)
)

print(model_summary)
write.csv(model_summary, here("Results", "models_v4", "model_summary.csv"))
message("\nModel summary saved to: Results/models_v4/model_summary.csv")

# Summary by model type
message("\n========== Model Summary by Type ==========")
summary_by_type <- model_summary |>
  group_by(model_type) |>
  summarise(
    total = n(),
    fitted = sum(fitted),
    production = sum(is_production),
    has_loo = sum(has_loo),
    .groups = "drop"
  )
print(summary_by_type)

message("\n========== Script complete ==========")
message("The comparison includes both spatial (2.4) and non-spatial (2.2) models")
message("This allows you to see the improvement that spatial GP provides (~30 ELPD)")

# 10. Targeted Fitting (Refit Selected Models) ---------------------------
# Use this section to run ONLY the models you choose.
# This is useful for:
#   a) Fitting a few models for the first time without running the whole suite.
#   b) Upgrading a 'screen' model to 'production' by adding it to production_models.

# To upgrade a model to production:
# 1. Add its name to production_models (e.g., production_models <- c("m1centrality_within_spatial"))
# 2. Add its name to models_to_run below.

models_to_run <- production_models

if (length(models_to_run) > 0) {
  message("\n========== Running Targeted Models ==========")
  for (m_name in models_to_run) {
    if (is.null(model_formulas[[m_name]])) {
      warning("Model name not found in registry: ", m_name)
      next
    }
    
    # Fit or refit (fit_model_auto handles checking if already exists)
    # If m_name is in production_models, it fits/loads [m_name]_final.rds
    # If not, it fits/loads [m_name].rds
    all_models[[m_name]] <- fit_model_auto(model_formulas[[m_name]], m_name)
    
    # Update/Compute LOO for this specific model
    message("Computing LOO for: ", m_name)
    available_loos[[m_name]] <- loo(all_models[[m_name]])
  }
  
  # Final re-comparison if targeted models were run
  if (length(available_loos) >= 2) {
    message("\n========== UPDATED MODEL COMPARISON ==========")
    comp_updated <- loo_compare(available_loos)
    print(comp_updated)
    write.csv(comp_updated, here("Results", "models_v4", "loo_comparison_updated.csv"))
  }
}

# 11. Detailed Exploration of Top Models ---------------------------------
message("\n========== Section 11: Exploring Top Models ==========")

# Ensure we have the production models loaded
# production_models are defined at the top
top_models_list <- list()
for (m_name in production_models) {
  mod <- load_model_auto(m_name)
  if (!is.null(mod)) {
    top_models_list[[m_name]] <- mod
  }
}

if (length(top_models_list) > 0) {
  # 11.1 Extract Coefficients and Parameters -------------------------------
  model_metrics <- list()
  
  for (m_name in names(top_models_list)) {
    message("Analyzing model: ", m_name)
    mod <- top_models_list[[m_name]]
    
    # Bayes R2
    message("  Computing Bayes R2...")
    r2 <- tryCatch(bayes_R2(mod), error = function(e) NA)
    
    # Fixed effects
    fe <- fixef(mod) |> 
      as.data.frame() |> 
      rownames_to_column("parameter") |>
      mutate(model = m_name, type = "fixed_effect")
    
    # GP and RE parameters
    post_sum <- posterior_summary(mod) |> 
      as.data.frame() |> 
      rownames_to_column("parameter") |>
      filter(grepl("^sd_|^phi|^sdgp|^lscale", parameter)) |>
      mutate(model = m_name, type = "hyperparameter")
    
    # Combine
    model_metrics[[m_name]] <- bind_rows(fe, post_sum) |>
      mutate(r2_mean = if(is.matrix(r2)) r2[1,1] else r2)
  }
  
  summary_table <- bind_rows(model_metrics) |>
    select(model, parameter, type, Estimate, Est.Error, Q2.5, Q97.5, r2_mean)
  
  write.csv(summary_table, here("Results", "models_v4", "top_models_exploration.csv"), row.names = FALSE)
  summary_table = read.csv(here("Results", "models_v4", "top_models_exploration.csv"))
  print(summary_table)
  message("Summary table saved to: Results/models_v4/top_models_exploration.csv")
  
  # 11.2 Visual Comparison: Key Predictors ----------------------------------
  # Identify niche-related parameters
  message("Generating coefficient comparison plot...")
  
  plot_data <- summary_table |>
    filter(type == "fixed_effect") |>
    filter(!grepl("Intercept|assay_group", parameter))
    
  if (nrow(plot_data) > 0) {
    coef_plot <- ggplot(plot_data, aes(x = Estimate, y = parameter, color = model)) +
      geom_pointrange(aes(xmin = Q2.5, xmax = Q97.5), 
                      position = position_dodge(width = 0.5)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      labs(
        title = "Niche Effect Comparison across Top Models",
        subtitle = "Estimates with 95% Credibility Intervals (logit scale)",
        x = "Estimate (log-odds)",
        y = "Predictor",
        color = "Model"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "vertical"
      )
    
    ggsave(here("Results", "models_v4", "top_models_coef_comparison.png"), 
           coef_plot, width = 10, height = 7)
    message("Comparison plot saved to: Results/models_v4/top_models_coef_comparison.png")
  }
  
  # 11.3 Performance Summary ------------------------------------------------
  performance_summary <- summary_table |>
    group_by(model) |>
    summarise(Bayes_R2 = first(r2_mean), .groups = "drop") |>
    arrange(desc(Bayes_R2))
  
  message("\nPerformance Summary (Bayes R2):")
  print(performance_summary)
  
  # 11.4 Effect Plots (Prevalence Scale) ------------------------------------
  message("\nGenerating prevalence-scale effect plots...")
  dir.create(here("Results", "models_v4", "effects"), showWarnings = FALSE)
  
  for (m_name in names(top_models_list)) {
    mod <- top_models_list[[m_name]]
    
    # Get the main predictors (excluding Intercept, assay_group)
    main_preds <- summary_table |> 
      filter(model == m_name, type == "fixed_effect", !grepl("Intercept|assay_group", parameter)) |>
      pull(parameter)
    
    if (length(main_preds) > 0) {
      # For interaction models, we want to plot the interaction
      if (any(grepl(":", main_preds))) {
        inter_pred <- main_preds[grepl(":", main_preds)]
        # brms conditional_effects needs the base variables, not the colon string
        base_vars <- unlist(strsplit(inter_pred, ":"))
        eff_plot <- conditional_effects(mod, effects = paste0(base_vars[1], ":", base_vars[2]), 
                                       re_formula = NA, method = "posterior_epred")
      } else {
        # Single predictors
        eff_plot <- conditional_effects(mod, effects = main_preds, 
                                       re_formula = NA, method = "posterior_epred")
      }
      
      # Save individual model effects
      p_combined <- plot(eff_plot, plot = FALSE)
      if (length(p_combined) > 0) {
        # Style the plots
        styled_plots <- lapply(p_combined, function(p) {
          p + theme_minimal() + labs(y = "Predicted Prevalence")
        })
        
        final_eff_plot <- do.call(grid.arrange, c(styled_plots, ncol = 1))
        ggsave(here("Results", "models_v4", "effects", paste0(m_name, "_effects.png")), 
               final_eff_plot, width = 8, height = 5 * length(styled_plots))
      }
    }
  }
  
  # 11.5 Random Effects Caterpillar Plot ------------------------------------
  message("\nGenerating random effects caterpillar plot...")
  # Use the top model (first in performance summary)
  best_model_name <- performance_summary$model[1]
  mod_best <- top_models_list[[best_model_name]]
  
  re <- ranef(mod_best, summary = TRUE)$host_species |>
    as.data.frame() |>
    rownames_to_column("host_species")
  
  if (nrow(re) > 0) {
    re_plot <- ggplot(re, aes(x = reorder(host_species, Estimate.Intercept), y = Estimate.Intercept)) +
      geom_pointrange(aes(ymin = Q2.5.Intercept, ymax = Q97.5.Intercept), 
                      color = "steelblue", alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(
        title = paste("Species Random Effects (Intercepts):", best_model_name),
        subtitle = "Species-level deviations in log-odds prevalence",
        x = "Host Species",
        y = "Estimate (log-odds deviation)"
      ) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 7))
    
    ggsave(here("Results", "models_v4", "species_random_effects.png"), 
           re_plot, width = 10, height = max(6, nrow(re)/4))
    message("Random effects plot saved to: Results/models_v4/species_random_effects.png")
  }
  
  # 11.6 LOO Diagnostic Plots -----------------------------------------------
  message("\nGenerating LOO diagnostic plots...")
  dir.create(here("Results", "models_v4", "loo_diagnostics"), showWarnings = FALSE)
  
  for (m_name in names(top_models_list)) {
    mod <- top_models_list[[m_name]]
    # Compute LOO if not already in available_loos
    loo_res <- available_loos[[m_name]]
    if (is.null(loo_res)) {
       try(loo_res <- loo(mod))
    }
    
    if (!is.null(loo_res)) {
      png(here("Results", "models_v4", "loo_diagnostics", paste0(m_name, "_pareto_k.png")), 
          width = 800, height = 600)
      plot(loo_res, main = paste("Pareto-k Diagnostic:", m_name))
      dev.off()
    }
  }
  
  # 11.7 Model Prediction Comparison ----------------------------------------
  message("\nGenerating model prediction comparison plot...")
  # We want to see how predictions change across a range of a common niche variable
  # Let's use Suitability (within or raw) as it's common in many models
  
  # Define a range for Suitability (z-scaled)
  suit_range <- seq(-3, 3, length.out = 100)
  
  # Create a prediction grid
  # We'll fix other vars at 0 (mean) and assay_group at its first level
  # Use data.frame instead of expand.grid to avoid 100^3 combinatorial explosion
  pred_grid_base <- data.frame(
    Suitability_within_z = suit_range,
    Suitability_orderNorm = suit_range,
    centrality_within_z = suit_range,
    number_tested = 1, # Required for trials() models to return probabilities
    assay_group = factor(levels(scaled_covars$assay_group)[1], levels = levels(scaled_covars$assay_group)),
    host_species = scaled_covars$host_species[1], # Will be ignored if re_formula=NA
    lon_z = 0,
    lat_z = 0,
    Marginality_within_z = 0,
    Marginality_orderNorm = 0,
    # Add other possible niche vars at 0 to satisfy models that might include them
    Specificity_within_z = 0,
    Specificity_orderNorm = 0,
    Centroid_d_within_z = 0,
    Centroid_d_orderNorm = 0,
    Boundary_d_within_z = 0,
    Boundary_d_orderNorm = 0
  )
  
  all_preds <- list()
  for (m_name in names(top_models_list)) {
    mod <- top_models_list[[m_name]]
    
    # Check if model has suitability or centrality
    has_suit <- any(grepl("Suitability|centrality", variables(mod)))
    if (!has_suit) next
    
    # Determine which predictor to use for the X-axis in this model
    if (any(grepl("centrality_within_z", variables(mod)))) {
      x_var <- "centrality_within_z"
    } else if (any(grepl("Suitability_within_z", variables(mod)))) {
      x_var <- "Suitability_within_z"
    } else {
      x_var <- "Suitability_orderNorm"
    }
    
    # Predict
    preds <- tryCatch({
      as.data.frame(fitted(mod, newdata = pred_grid_base, re_formula = NA)) |>
        mutate(x = suit_range, model = m_name, predictor = x_var)
    }, error = function(e) NULL)
    
    if (!is.null(preds)) all_preds[[m_name]] <- preds
  }
  
  if (length(all_preds) > 0) {
    comp_pred_data <- bind_rows(all_preds)
    
    comp_pred_plot <- ggplot(comp_pred_data, aes(x = x, y = Estimate, color = model)) +
      geom_line(size = 1.2) +
      geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = model), alpha = 0.1, color = NA) +
      labs(
        title = "Model Prediction Comparison",
        subtitle = "Predicted prevalence across range of niche variable (Suitability/Centrality)",
        x = "Predictor Value (SD)",
        y = "Predicted Prevalence",
        fill = "Model",
        color = "Model"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.direction = "vertical")
    
    ggsave(here("Results", "models_v4", "model_prediction_comparison.png"), 
           comp_pred_plot, width = 10, height = 7)
    message("Model comparison plot saved to: Results/models_v4/model_prediction_comparison.png")
  }

} else {
  message("No production models found for exploration. Ensure they are fitted and saved in Results/models_v4/")
}



