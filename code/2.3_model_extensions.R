## Model extensions to reduce ties / sharpen model distinctions
## ------------------------------------------------------------
## Motivation (from results.md + code/2.2_models.R):
## - When ΔELPD is small relative to SE, models are effectively tied under PSIS-LOO.
## - This often means either (a) the true predictive gain is genuinely tiny, or
##   (b) the evaluation unit is mismatched to the scientific question (e.g.,
##       predicting new *rows* vs new *species*), or
##   (c) important structure is missing (spatial clustering, pathogen heterogeneity),
##       inflating pointwise variability and SE.
##
## This script focuses on changes that can make comparisons more informative:
## - Spatial effects using longitude/latitude
## - Group-wise cross-validation (e.g., leave-one-host-species-out)
## - Additional grouping structure (pathogen species / host–pathogen pair)
## - Optional non-linear effects (splines)
##
## Run from the project root:
##   Rscript "code/2.3_model_extensions.R"
##

# 1. Packages -------------------------------------------------------------
pacman::p_load(tidyverse, brms, loo, here)

options(mc.cores = max(1, parallel::detectCores()))

# 2. Load data (same standardized dataset used in 2.2) -------------------
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

# 2.1 Within- vs between-host-species decomposition (central to the hypothesis) ----
# Your core question is essentially within-host-species:
# "At sites that are more central/preferred for a given host species, is prevalence different?"
#
# A simple way to target that estimand is to decompose predictors into:
# - between-species component: each species' mean predictor value across sampled sites
# - within-species component: deviation from that species mean at each site
#
# This helps avoid conflating "species that tend to occur in X environments" with
# "within a species, more-central sites have different prevalence".

vars_for_within_between <- c(
  "prob_occur_orderNorm",
  "Suitability_orderNorm",
  "Marginality_orderNorm"
)

missing_vars <- setdiff(vars_for_within_between, names(scaled_covars))
if (length(missing_vars) > 0) {
  stop("Missing required predictors in dat_clean_agg.rds: ", paste(missing_vars, collapse = ", "))
}

species_means_raw <- scaled_covars |>
  group_by(host_species) |>
  summarise(
    across(all_of(vars_for_within_between), ~ mean(.x, na.rm = TRUE), .names = "{.col}_between_raw"),
    .groups = "drop"
  )

species_means_scaled <- species_means_raw |>
  mutate(
    across(
      ends_with("_between_raw"),
      ~ as.numeric(scale(.x)),
      .names = "{.col}_z"
    )
  )

scaled_covars <- scaled_covars |>
  left_join(species_means_scaled, by = "host_species") |>
  mutate(
    # Within-species deviations (raw) then z-score across rows for stable priors
    prob_occur_within_raw = prob_occur_orderNorm - prob_occur_orderNorm_between_raw,
    Suitability_within_raw = Suitability_orderNorm - Suitability_orderNorm_between_raw,
    Marginality_within_raw = Marginality_orderNorm - Marginality_orderNorm_between_raw,
    prob_occur_within_z = as.numeric(scale(prob_occur_within_raw)),
    Suitability_within_z = as.numeric(scale(Suitability_within_raw)),
    Marginality_within_z = as.numeric(scale(Marginality_within_raw))
  )

# One-dimensional "centrality score" (optional) to match the hypothesis more directly.
# Higher score = more central (high suitability AND low marginality).
scaled_covars <- scaled_covars |>
  mutate(
    centrality_score_raw = Suitability_orderNorm - Marginality_orderNorm,
    centrality_between_raw = Suitability_orderNorm_between_raw - Marginality_orderNorm_between_raw,
    centrality_within_raw = centrality_score_raw - centrality_between_raw,
    centrality_within_z = as.numeric(scale(centrality_within_raw)),
    centrality_between_z = as.numeric(scale(centrality_between_raw))
  )

stopifnot(
  all(is.finite(scaled_covars$centrality_within_z)),
  all(is.finite(scaled_covars$centrality_between_z))
)

# 3. Priors ---------------------------------------------------------------
# Keep same baseline priors as 2.2, plus GP priors when we use spatial terms.
priors_beta_binom <- prior(normal(0, 1.5), class = "Intercept") +
  prior(normal(0, 1), class = "b") +
  prior(exponential(1), class = "sd") +
  prior(exponential(1), class = "phi")

priors_beta_binom_spatial <- priors_beta_binom +
  # GP hyperparameters (used only when gp() terms exist)
  prior(exponential(1), class = "sdgp") +
  prior(exponential(1), class = "lscale")

# 4. Helpers --------------------------------------------------------------
load_if_exists <- function(path_no_ext) {
  rds_path <- paste0(path_no_ext, ".rds")
  if (file.exists(rds_path)) return(readRDS(rds_path))
  NULL
}

# Create fold IDs for group-wise CV (e.g., each host_species held out together).
# Returns an integer vector of length nrow(data) mapping each row to a fold.
make_group_folds <- function(data, group_var) {
  group_vec <- data[[group_var]]
  stopifnot(length(group_vec) == nrow(data))
  group_levels <- unique(group_vec)
  group_to_fold <- setNames(seq_along(group_levels), as.character(group_levels))
  as.integer(group_to_fold[as.character(group_vec)])
}

run_group_kfold <- function(model, data, group_var) {
  folds <- make_group_folds(data, group_var)
  brms::kfold(model, folds = folds, save_fits = TRUE)
}

# 5. Candidate model formulas --------------------------------------------
# NOTE: we keep assay_group as a fixed effect (few levels) and host_species as RE.
f_m0 <- bf(
  number_positive | trials(number_tested) ~ 1 + assay_group + (1 | host_species)
)

# Add pathogen species heterogeneity (often large; can reduce unexplained variation).
f_m0_pathogen <- bf(
  number_positive | trials(number_tested) ~ 1 + assay_group +
    (1 | host_species) + (1 | pathogen_species_cleaned)
)

# Add host–pathogen pair intercepts (captures pair-specific baseline prevalence).
# This can be a big structural improvement if many host–pathogen combinations repeat.
f_m0_pair <- bf(
  number_positive | trials(number_tested) ~ 1 + assay_group +
    (1 | host_species) + (1 | pathogen_species_cleaned) +
    (1 | host_species:pathogen_species_cleaned)
)

# Spatial baseline via low-rank GP on standardized lon/lat.
# k controls approximation rank; smaller is faster but less flexible.
f_m0_spatial <- bf(
  number_positive | trials(number_tested) ~ 1 + assay_group + (1 | host_species) +
    gp(lon_z, lat_z, k = 50)
)

# Example niche model (replace prob_occur_orderNorm with any niche metric).
f_m1prob <- bf(
  number_positive | trials(number_tested) ~ 1 + prob_occur_orderNorm +
    assay_group + (1 | host_species)
)

# Spatial niche model: niche predictor + spatial term.
f_m1prob_spatial <- bf(
  number_positive | trials(number_tested) ~ 1 + prob_occur_orderNorm +
    assay_group + (1 | host_species) + gp(lon_z, lat_z, k = 50)
)

# Within-species versions (targets the centrality hypothesis more directly)
f_m1centrality_within <- bf(
  number_positive | trials(number_tested) ~ 1 + centrality_within_z +
    assay_group + (1 | host_species)
)

# Add between-species component as an explanatory covariate for species baselines.
# This can explain part of the host_species intercept variation.
f_m1centrality_within_between <- bf(
  number_positive | trials(number_tested) ~ 1 + centrality_within_z + centrality_between_z +
    assay_group + (1 | host_species)
)

f_m1centrality_within_spatial <- bf(
  number_positive | trials(number_tested) ~ 1 + centrality_within_z +
    assay_group + (1 | host_species) + gp(lon_z, lat_z, k = 50)
)

# Interaction between within-species suitability and marginality
f_interaction_within <- bf(
  number_positive | trials(number_tested) ~ 1 + 
    Suitability_within_z * Marginality_within_z +
    assay_group + (1 | host_species) + gp(lon_z, lat_z, k=20)
)

# 6. Fit / load models ----------------------------------------------------
# We default to loading existing fits; only fit new ones if you flip fit_new_models.

models_dir <- here("Results", "models_v3")
if (!dir.exists(models_dir)) dir.create(models_dir, recursive = TRUE)

m0_base_bb <- load_if_exists(here("Results", "models_v1", "m0_base_bb"))
m1prob_v2 <- load_if_exists(here("Results", "models_v2", "m1prob_v2"))

if (isTRUE(fit_new_models)) {
  if (is.null(m0_base_bb)) {
    m0_base_bb <- brm(
      formula = f_m0,
      family = beta_binomial(link = "logit"),
      data = scaled_covars,
      prior = priors_beta_binom,
      iter = 2000,
      save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      chains = 4,
      cores = 4,
      seed = 123,
      file = here("Results", "models_v1", "m0_base_bb")
    )
  }
  if (is.null(m1prob_v2)) {
    m1prob_v2 <- brm(
      formula = f_m1prob,
      family = beta_binomial(link = "logit"),
      data = scaled_covars,
      prior = priors_beta_binom,
      iter = 2000,
      save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      chains = 4,
      cores = 4,
      seed = 123,
      file = here("Results", "models_v2", "m1prob_v2")
    )
  }
}

# New extension models (optional to fit; safe to leave fit_new_models = FALSE)
m0_spatial <- load_if_exists(here("Results", "models_v3", "m0_base_bb_spatial"))
m1prob_spatial <- load_if_exists(here("Results", "models_v3", "m1prob_v2_spatial"))
m0_pathogen <- load_if_exists(here("Results", "models_v3", "m0_base_bb_pathogen"))
m0_pair <- load_if_exists(here("Results", "models_v3", "m0_base_bb_pair"))

centrality_within <- load_if_exists(here("Results", "models_v3", "m1centrality_within"))
centrality_within_between <- load_if_exists(here("Results", "models_v3", "m1centrality_within_between"))
centrality_within_spatial <- load_if_exists(here("Results", "models_v3", "m1centrality_within_spatial"))
interaction_within <- load_if_exists(here("Results", "models_v3", "interaction_within"))

  if (is.null(m0_spatial)) {
    # m0_spatial <- brm(
    #   formula = f_m0_spatial,
    #   family = beta_binomial(link = "logit"),
    #   data = scaled_covars,
    #   prior = priors_beta_binom_spatial,
    #   iter = 2000,
    #   save_pars = save_pars(all = TRUE),
    #   control = list(adapt_delta = 0.99, max_treedepth = 15),
    #   chains = 4,
    #   cores = 4,
    #   seed = 123,
    #   file = here("Results", "models_v3", "m0_base_bb_spatial")
    # )
    # cheaper GP rank
    f_m0_spatial <- bf(
      number_positive | trials(number_tested) ~ 1 + assay_group + (1 | host_species) +
        gp(lon_z, lat_z, k = 20)
    )
    
    m0_spatial <- brm(
      formula = f_m0_spatial,
      family = beta_binomial(link = "logit"),
      data = scaled_covars,
      prior = priors_beta_binom_spatial,
      iter = 1000,
      warmup = 500,
      save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.97, max_treedepth = 12),
      chains = 2,
      cores = 1,          # prevents SIGPIPE/socket-worker issues
      init = 0,
      init_r = 0.1,
      seed = 123,
      file = here("Results", "models_v3", "m0_base_bb_spatial_k20_screen")
    )  
  
  }
  if (is.null(m1prob_spatial)) {
    m1prob_spatial <- brm(
      formula = f_m1prob_spatial,
      family = beta_binomial(link = "logit"),
      data = scaled_covars,
      prior = priors_beta_binom_spatial,
      iter = 2000,
      save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      chains = 4,
      cores = 4,
      seed = 123,
      file = here("Results", "models_v3", "m1prob_v2_spatial")
    )
  }
  if (is.null(m0_pathogen)) {
    m0_pathogen <- brm(
      formula = f_m0_pathogen,
      family = beta_binomial(link = "logit"),
      data = scaled_covars,
      prior = priors_beta_binom,
      iter = 2000,
      save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      chains = 4,
      cores = 4,
      seed = 123,
      file = here("Results", "models_v3", "m0_base_bb_pathogen")
    )
  }
  if (is.null(m0_pair)) {
    m0_pair <- brm(
      formula = f_m0_pair,
      family = beta_binomial(link = "logit"),
      data = scaled_covars,
      prior = priors_beta_binom,
      iter = 2000,
      save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      chains = 4,
      cores = 4,
      seed = 123,
      file = here("Results", "models_v3", "m0_base_bb_pair")
    )
  }

  if (is.null(centrality_within)) {
    centrality_within <- brm(
      formula = f_m1centrality_within,
      family = beta_binomial(link = "logit"),
      data = scaled_covars,
      prior = priors_beta_binom,
      iter = 2000,
      save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      chains = 4,
      cores = 4,
      seed = 123,
      file = here("Results", "models_v3", "m1centrality_within")
    )
  }

  if (is.null(centrality_within_between)) {
    centrality_within_between <- brm(
      formula = f_m1centrality_within_between,
      family = beta_binomial(link = "logit"),
      data = scaled_covars,
      prior = priors_beta_binom,
      iter = 2000,
      save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      chains = 4,
      cores = 4,
      seed = 123,
      file = here("Results", "models_v3", "m1centrality_within_between")
    )
  }

  # Centrality after spatial control 
  if (is.null(centrality_within_spatial)) {
    centrality_within_spatial <- brm(
      formula = f_m1centrality_within_spatial,
      family = beta_binomial(link = "logit"),
      data = scaled_covars,
      prior = priors_beta_binom_spatial,
      iter = 2000,
      save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      chains = 4,
      cores = 4,
      seed = 123,
      file = here("Results", "models_v3", "m1centrality_within_spatial")
    )
  }
  # Cheaper version of centrality after spatial control
  # Make sure the formula uses the cheaper GP rank (k=20)
f_m1centrality_within_spatial <- bf(
  number_positive | trials(number_tested) ~ 1 + centrality_within_z +
    assay_group + (1 | host_species) + gp(lon_z, lat_z, k = 20)
)

centrality_within_spatial <- brm(
  formula = f_m1centrality_within_spatial,
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom_spatial,
  iter = 1000,
  warmup = 500,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.97, max_treedepth = 12),
  chains = 2,
  cores = 1,
  init = 0,
  init_r = 0.1,
  seed = 123,
  file = here("Results", "models_v3", "m1centrality_within_spatial_k20_screen")
)

centrality_within_spatial_refit <- update(
  centrality_within_spatial,
  iter = 2000,
  warmup = 1000,
  chains = 2,
  cores = 1
)


if (is.null(interaction_within)) {
  interaction_within <- brm(
    formula = f_interaction_within,
    family = beta_binomial(link = "logit"),
    data = scaled_covars,
    prior = priors_beta_binom_spatial,
    iter = 2000,
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    chains = 4,
    cores = 4,
    seed = 123,
    file = here("Results", "models_v3", "interaction_within")
  )
}
loo_interaction_within <- loo(interaction_within)



loo_cent_spat2 <- loo(centrality_within_spatial_refit)
loo_spat2 <- loo(m0_spatial)  # (or reuse if unchanged)
loo_compare(loo_spat2, loo_cent_spat2)

loo_cent_spat <- loo(centrality_within_spatial)
loo_compare(loo_spat, loo_cent_spat)
fixef(centrality_within_spatial)["centrality_within_z", ]

# 7. Comparisons ----------------------------------------------------------
# 7.1 Load models from 2.2 (models_v1 and models_v2) for comprehensive comparison
message("\nLoading models from code/2.2_models.R for comparison...\n")

# Load baseline models from v1
m0_controls_v1 <- load_if_exists(here("Results", "models_v1", "m0_controls"))
m0_beta_binom_v1 <- load_if_exists(here("Results", "models_v1", "m0_beta_binom"))

# Load M1 models from v1 (with controls: assay + host_family + pathogen_family)
m1prob_v1 <- load_if_exists(here("Results", "models_v1", "m1prob"))
m1suit_v1 <- load_if_exists(here("Results", "models_v1", "m1suit"))
m1marg_v1 <- load_if_exists(here("Results", "models_v1", "m1marg"))
m1spec_v1 <- load_if_exists(here("Results", "models_v1", "m1spec"))
m1centr_v1 <- load_if_exists(here("Results", "models_v1", "m1centr"))
m1bound_v1 <- load_if_exists(here("Results", "models_v1", "m1bound"))

# Load M1 models from v2 (minimal baseline: assay + host_species only)
m1prob_v2_old <- load_if_exists(here("Results", "models_v2", "m1prob_v2"))
m1suit_v2 <- load_if_exists(here("Results", "models_v2", "m1suit_v2"))
m1marg_v2 <- load_if_exists(here("Results", "models_v2", "m1marg_v2"))
m1spec_v2 <- load_if_exists(here("Results", "models_v2", "m1spec_v2"))
m1centr_v2 <- load_if_exists(here("Results", "models_v2", "m1centr_v2"))
m1bound_v2 <- load_if_exists(here("Results", "models_v2", "m1bound_v2"))

# Load M2 models from v2 (multivariate models)
m2a_v2 <- load_if_exists(here("Results", "models_v2", "m2a_v2"))
m2b_v2 <- load_if_exists(here("Results", "models_v2", "m2b_v2"))
m2c_v2 <- load_if_exists(here("Results", "models_v2", "m2c_v2"))
m2d_v2 <- load_if_exists(here("Results", "models_v2", "m2d_v2"))

# 7.2 PSIS-LOO (row-wise). This answers: "predict a new row like the observed rows".
available_loos <- list()

# Models from 2.3 (current script)
if (!is.null(m0_base_bb)) available_loos$m0_base_bb <- loo(m0_base_bb)
if (!is.null(m1prob_v2)) available_loos$m1prob_v2 <- loo(m1prob_v2)
if (!is.null(m0_spatial)) available_loos$m0_spatial <- loo(m0_spatial)
if (!is.null(m1prob_spatial)) available_loos$m1prob_spatial <- loo(m1prob_spatial)
if (!is.null(m0_pathogen)) available_loos$m0_pathogen <- loo(m0_pathogen)
if (!is.null(m0_pair)) available_loos$m0_pair <- loo(m0_pair)
if (!is.null(centrality_within)) available_loos$m1centrality_within <- loo(centrality_within)
if (!is.null(centrality_within_between)) {
  available_loos$m1centrality_within_between <- loo(centrality_within_between)
}
if (!is.null(centrality_within_spatial)) {
  available_loos$m1centrality_within_spatial <- loo(centrality_within_spatial)
  available_loos$m1centrality_within_spatial_refit <- loo_cent_spat2
}
if (!is.null(interaction_within)) {
  available_loos$interaction_within <- loo(interaction_within)
}

# Add models from 2.2 (v1 and v2)
message("Computing LOO for models from code/2.2_models.R...")
if (!is.null(m0_controls_v1)) available_loos$m0_controls_v1 <- loo(m0_controls_v1)
if (!is.null(m0_beta_binom_v1)) available_loos$m0_beta_binom_v1 <- loo(m0_beta_binom_v1)

# M1 v1 (with controls)
if (!is.null(m1prob_v1)) available_loos$m1prob_v1 <- loo(m1prob_v1)
if (!is.null(m1suit_v1)) available_loos$m1suit_v1 <- loo(m1suit_v1)
if (!is.null(m1marg_v1)) available_loos$m1marg_v1 <- loo(m1marg_v1)
if (!is.null(m1spec_v1)) available_loos$m1spec_v1 <- loo(m1spec_v1)
if (!is.null(m1centr_v1)) available_loos$m1centr_v1 <- loo(m1centr_v1)
if (!is.null(m1bound_v1)) available_loos$m1bound_v1 <- loo(m1bound_v1)

# M1 v2 (minimal baseline)
if (!is.null(m1prob_v2_old)) available_loos$m1prob_v2_old <- loo(m1prob_v2_old)
if (!is.null(m1suit_v2)) available_loos$m1suit_v2 <- loo(m1suit_v2)
if (!is.null(m1marg_v2)) available_loos$m1marg_v2 <- loo(m1marg_v2)
if (!is.null(m1spec_v2)) available_loos$m1spec_v2 <- loo(m1spec_v2)
if (!is.null(m1centr_v2)) available_loos$m1centr_v2 <- loo(m1centr_v2)
if (!is.null(m1bound_v2)) available_loos$m1bound_v2 <- loo(m1bound_v2)

# M2 v2 (multivariate)
if (!is.null(m2a_v2)) available_loos$m2a_v2 <- loo(m2a_v2)
if (!is.null(m2b_v2)) available_loos$m2b_v2 <- loo(m2b_v2)
if (!is.null(m2c_v2)) available_loos$m2c_v2 <- loo(m2c_v2)
if (!is.null(m2d_v2)) available_loos$m2d_v2 <- loo(m2d_v2)

message("\n========== COMPREHENSIVE MODEL COMPARISON ==========")
message("Comparing all models from 2.2 and 2.3...")
message("Total models in comparison: ", length(available_loos), "\n")

if (length(available_loos) >= 2) {
  comp <- loo_compare(available_loos)
  print(comp)
  
  # Also print a summary of the top 5 models
  message("\n========== TOP 5 MODELS ==========")
  if (nrow(comp) >= 5) {
    print(comp[1:5, ])
  }
}

# 7.3 Group-wise CV. This answers: "predict new host species" (often more meaningful).
# This can change the ranking substantially vs row-wise LOO.
run_groupwise_cv <- FALSE

if (isTRUE(run_groupwise_cv) && !is.null(m0_base_bb) && !is.null(m1prob_v2)) {
  kf_m0_species <- run_group_kfold(m0_base_bb, scaled_covars, "host_species")
  kf_m1_species <- run_group_kfold(m1prob_v2, scaled_covars, "host_species")
  print(kfold_compare(kf_m0_species, kf_m1_species))
}

# 7.4 Model averaging (when ties persist): stacking weights.
run_stacking <- FALSE

if (isTRUE(run_stacking) && length(available_loos) >= 2) {
  w <- loo_model_weights(available_loos, method = "stacking")
  print(w)
}

# 8. Practical guidance (what usually helps) -----------------------------
message(
  "\nNext steps to reduce ties / SE:\n",
  "- Try group-wise CV (host_species) to match your generalization target.\n",
  "- Add spatial structure: compare m0_base_bb vs m0_spatial; then m1prob_v2 vs m1prob_spatial.\n",
  "- Target the centrality hypothesis directly: use within-species centrality (m1centrality_within).\n",
  "- Optional: add centrality_between_z to explain species baselines (m1centrality_within_between).\n",
  "- Add pathogen heterogeneity (pathogen_species_cleaned) and/or host–pathogen pair intercepts.\n",
  "- If ties remain, report stacking weights instead of a single 'winner'.\n"
)


