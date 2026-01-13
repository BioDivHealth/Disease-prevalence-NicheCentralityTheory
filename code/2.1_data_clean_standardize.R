# Data clean & standardize

#--------------------------------------#
# Predict zoonotic pathogen prevalence #
#--------------------------------------#

# CODEX CONVERSATIONS:
# Variable scaling: codex resume 019a5e7a-b5ff-7ec1-849a-8980cdc55294
# Priors discussions: codex resume 019a5f74-5ed5-7c52-94a5-f635e3e4c5c1


# 0. Script purpose ----
# - Load SDM estimates, host & pathogen data
# - Clean & format data
# - Formulate Bayesian GLMs to predict prevalence given covariates
# - Visualise outputs


# 1. Load packages ----
pacman::p_load(sf, tidyverse, brms, bayesplot, bestNormalize)


# 2. Load data ------ ----------------------------------------
arha      <- readRDS("Data/Project_ArHa_database_2025-09-18.rds")  # ArHA
arha_path <- arha$pathogen
arha_host <- arha$host
dat       <- read.csv("Data/Full_data.csv") # SDM estimates + ArHA data


# 3. Wrangle data --------------------------------------------

# Join data: get additional columns from ArHA data
dat_joined <- dat %>% 
  # Rename columns
  rename(decimalLatitude = bp_Y, decimalLongitude = bp_X, prob_occur = mean) %>% 
  # Format date
  mutate(start_date = ymd(start_date), end_date = ymd(end_date)) %>% 
  # Columns for host & study ID
  tidyr::separate_wider_delim(cols = id.s, delim = " ", names = c("host_record_id", "study_id")) %>% 
  # Pathogen & host columns
  left_join(arha_path) %>%
  left_join(arha_host)

# Clean data
dat_clean <- dat_joined %>% 
  filter(
    # Remove NAs in relevant columns
    if_all(c(number_positive, number_tested, host_family, pathogen_family, 
             assay, prob_occur, Marginality, Specificity, Suitability, 
             Centroid_d, Boundary_d), ~ !is.na(.)),
    tolower(assay) != "missing",
    number_tested  != 0,
    # Keep desired coordinate resolution
    coordinate_resolution_processed %in% c("site", "village", "town", "city", "adm3")
  ) %>% 
  mutate(
    # Transform no. tested & no. positive to integers
    number_positive = as.integer(number_positive),
    number_tested   = as.integer(number_tested),
    # Group by assay type
    assay_group = case_when(
      assay == "Serology" | assay == "Western Blot" ~ "Serology", 
      assay != "Serology" ~ "Culture/molecular"), 
    # Group by temporal resolution - ** can use this later for including temporal effects **
    date_interval = interval(start_date, end_date),
    month_diff = date_interval %/% months(1), 
    days_diff = date_interval %/% days(1),
    month_diff_group = case_when(
      month_diff < 1                    ~ "<1",
      month_diff >= 1 & month_diff < 3  ~ "1-3",
      month_diff >= 3 & month_diff < 6  ~ "3-6",
      month_diff >= 6 & month_diff < 12 ~ "6-12",
      month_diff >= 12                  ~ ">12"),
    month_diff_group = factor(month_diff_group, levels = c("<1", "1-3", "3-6", "6-12", ">12"))
  ) %>% 
  reframe(number_positive, number_tested, number_negative,number_inconclusive,host_family, pathogen_family,
          assay_group, prob_occur, Marginality, Specificity, Suitability, month_diff_group, month_diff, date_interval, days_diff, start_date, end_date,
          Centroid_d, Boundary_d, X.1, host_record_id, host_species, study_id, pathogen_record_id, pathogen_species_cleaned, longitude, latitude) %>% 
  # Exclude rows where its group has fewer than 20 observations
  group_by(host_family, pathogen_family, assay_group) %>% 
  filter(n() >= 20) %>% 
  ungroup()


# 3.1 Investigation & Aggregation of Duplicate Records -------------------------

# Aggregate rows that have the same lon/lat, host, pathogen and covariates
# Sum number tested, positive, negative, inconclusive
dat_clean_agg = dat_clean %>% 
  dplyr::select(host_family, pathogen_family, assay_group, host_species, pathogen_species_cleaned,
                                              prob_occur, Marginality, Specificity, Suitability, Centroid_d, Boundary_d,
                                              longitude, latitude, date_interval, days_diff,
                                              number_tested, number_positive, number_negative, number_inconclusive,
                                              host_record_id, pathogen_record_id, study_id) %>% 
  dplyr::group_by(host_family, pathogen_family, assay_group, host_species, pathogen_species_cleaned,
           prob_occur, Marginality, Specificity, Suitability, Centroid_d, Boundary_d,
           longitude, latitude, date_interval, days_diff) %>% 
  dplyr::summarise(number_tested = sum(number_tested, na.rm = TRUE),
            number_positive = sum(number_positive, na.rm = TRUE),
            number_negative = sum(number_negative, na.rm = TRUE),
            number_inconclusive = sum(number_inconclusive, na.rm = TRUE),
            rows = n(),
            .groups = "drop") %>% 
  group_by(host_species) %>% # Keep species for which more than X rows ARTUR
  filter(n() >= 15) %>%
  ungroup()


# 3.1.2 Temporal Resolution Diagnostics ----------------------------------------
# These plots evaluate the distribution of sampling intervals (days_diff).
# Since aggregation now preserves unique date intervals, we check for:
# - Common sampling durations (e.g., 1 week, 1 month)
# - Extreme outliers (multi-year studies)
# - Consistency across the dataset to inform temporal covariate scaling

dim(dat_clean_agg)
# 2538 by default

# Summary of sampling duration in days
summary(dat_clean_agg$days_diff)


df <- dat_clean_agg %>%
  filter(!is.na(days_diff))


# Plot 1: Distribution of shorter studies (< 1 year) using density
p_zoom_density <- ggplot(df, aes(x = days_diff)) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 7, boundary = 0,
    colour = "white", linewidth = 0.3
  ) +
  coord_cartesian(xlim = c(0, 370)) +
  scale_x_continuous(
    breaks = seq(0, 365, by = 30),
    minor_breaks = seq(0, 365, by = 7),
    expand = c(0, 0)
  ) +
  labs(
    title = "Time between start and end dates (zoomed to 1 year)",
    subtitle = "Histogram (7-day bins); y-axis is density",
    x = "Days (end − start)",
    y = "Density"
  ) +
  theme_classic(base_size = 13)

p_zoom_density

# Plot 2: Distribution of shorter studies (< 1 year) using counts
p_zoom_count <- ggplot(df, aes(x = days_diff)) +
  geom_histogram(
    binwidth = 7, boundary = 0,
    colour = "white", linewidth = 0.3, fill = "coral3"
  ) +
  coord_cartesian(xlim = c(0, 365)) +
  scale_x_continuous(
    breaks = seq(0, 365, by = 30),
    minor_breaks = seq(0, 365, by = 7),
    expand = c(0, 0)
  ) +
  labs(
    title = "Time between start and end dates (0–365 days)",
    subtitle = "Histogram (7-day bins)",
    x = "Days (end − start)",
    y = "Count"
  ) +
  theme_classic(base_size = 13)

p_zoom_count

# Plot 3: Full range distribution with yearly markers to identify multi-year studies
max_days <- max(df$days_diff, na.rm = TRUE)
year_breaks <- seq(0, ceiling(max_days / 365) * 365, by = 365)

p_full_linear <- ggplot(df, aes(x = days_diff)) +
  geom_histogram(
    binwidth = 60, boundary = 0,
    colour = "white", linewidth = 0.5, fill = "coral4"
  ) +
  geom_vline(
    xintercept = year_breaks,
    linewidth = 0.3,
    linetype = "dashed",
    alpha = 0.3
  ) +
  scale_x_continuous(
    breaks = year_breaks,
    labels = paste0(year_breaks / 365, "y"),
    expand = c(0, 0)
  ) +
  labs(
    title = "Time between start and end dates (full range)",
    subtitle = "Dashed lines mark each 365-day interval",
    x = "Days (end − start)",
    y = "Count"
  ) +
  theme_classic(base_size = 13)

p_full_linear

table(dat_clean_agg$host_species)
dim(table(dat_clean_agg$host_species))
# Should sum to 11399 - as that was initial total after cleaning
sum(dat_clean_agg$rows)
table(dat_clean_agg$host_species)

# 3.2 Adjust Suitability values that are exactly 0 -------------------------
# Issue: Many rows have Suitability = 0 exactly due to precision loss in 
#        calculating 1 - pchisq() instead of pchisq(..., lower.tail = FALSE).
#        These exact 0s are treated as ties by orderNorm transformation.
# Solution: Add tiny jitter around 1e-16 to values that are exactly 0.

n_zero_suit <- sum(dat_clean_agg$Suitability == 0)
cat(sprintf("Found %d rows with Suitability == 0 (%.2f%% of data)\n", 
            n_zero_suit, 100 * n_zero_suit / nrow(dat_clean_agg)))

if (n_zero_suit > 0) {
  set.seed(123)  # Reproducibility
  dat_clean_agg <- dat_clean_agg %>%
    mutate(
      Suitability = if_else(
        Suitability == 0,
        runif(n(), min = 1e-17, max = 1e-15),  # Tiny jitter around 1e-16
        Suitability
      )
    )
  cat("Adjusted exact 0s to tiny positive values (range: 1e-17 to 1e-15)\n")
}

# Summarise numeric constraints for key covariates (helps pick transformations)
covariate_vars <- c("prob_occur", "Marginality", "Specificity",
                    "Suitability", "Centroid_d", "Boundary_d")

summarise_numeric_constraints <- function(data, vars) {
  purrr::map_dfr(vars, function(var_name) {
    values <- data[[var_name]]
    rng <- range(values, na.rm = TRUE)
    tibble::tibble(
      variable = var_name,
      min = rng[1],
      max = rng[2],
      crosses_zero = rng[1] < 0 & rng[2] > 0,
      bounded_between_0_1 = rng[1] >= 0 & rng[2] <= 1,
      non_negative = rng[1] >= 0,
      strictly_positive = rng[1] > 0
    )
  })
}

covariate_constraints <- summarise_numeric_constraints(dat_clean_agg, covariate_vars)
print(covariate_constraints)

# sum(dat_clean$number_tested)
# [1] 83877

# 4. Scale covariates ----------------------------------------------------------

# Function to apply unit scaling to covariates (mean of 0, SD of 1) #ARTUR
unitScale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Give covariates unit scaling
scaled_covars <- dat_clean_agg %>% 
  mutate(mean_prob_occur = mean(prob_occur),
         mean_marg = mean(Marginality),
         sd_marg = sd(Marginality),
         mean_spec = mean(Specificity),
         sd_spec = sd(Specificity),
         mean_suit = mean(Suitability),
         sd_suit = sd(Suitability),
         mean_centroid_dist = mean(Centroid_d),
         sd_centroid_dist   = sd(Centroid_d),
         mean_boundary_dist = mean(Boundary_d),
         sd_boundary_dist   = sd(Boundary_d),
         across(c(prob_occur, Marginality, Specificity, Suitability, 
                  Centroid_d, Boundary_d), unitScale))

# Explore scales and distribution of covariates
covar_hist_data <- bind_rows(
  dat_clean_agg %>% 
    select(prob_occur, Marginality, Specificity, Suitability, Centroid_d, Boundary_d) %>% 
    mutate(scale = "Unscaled"),
  scaled_covars %>% 
    select(prob_occur, Marginality, Specificity, Suitability, Centroid_d, Boundary_d) %>% 
    mutate(scale = "Scaled")
) %>% 
  pivot_longer(cols = -scale, names_to = "variable", values_to = "value") %>% 
  mutate(
    scale = factor(scale, levels = c("Unscaled", "Scaled")),
    variable = factor(variable, levels = c("prob_occur", "Marginality", "Specificity", "Suitability", "Centroid_d", "Boundary_d"))
  )

covar_hist_plot <- covar_hist_data %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "#1f78b4", colour = "white") +
  facet_grid(rows = vars(scale), cols = vars(variable), scales = "free", switch = "y") +
  labs(
    x = NULL,
    y = "Count",
    title = "Distributions of unscaled and scaled covariates"
  ) +
  theme_bw() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank()
  )

print(covar_hist_plot)

# 4.1 New scaling functions - ARTUR ----------------------------------------------------
transform_labels <- c(
  arcsinh_x = "Arcsinh",
  boxcox = "Box-Cox",
  center_scale = "Center+Scale",
  double_reverse_log = "Double Reversed Log_b",
  exp_x = "Exp",
  log_x = "Log_b",
  sqrt_x = "Sqrt",
  yeojohnson = "Yeo-Johnson",
  orderNorm = "OrderNorm"
)

# Plot histogram centered around zero with consistent styling
hist_centered = function(x, main, fill = "#bdbdbd"){
  hist(x,
       xlim = c(-max(abs(x)), max(abs(x))),
       main = main,
       col = fill,
       border = "white")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Safe z-score even when SD is zero/NA
safe_scale_vec <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- stats::sd(x, na.rm = TRUE)
  if (is.na(sd_x) || sd_x == 0) {
    return(rep(0, length(x)))
  }
  (x - mean_x) / sd_x
}

# Logit transform with clipping away from 0/1
logit_clipped <- function(x, eps = 1e-6) {
  x_clipped <- pmin(pmax(x, eps), 1 - eps)
  stats::qlogis(x_clipped)
}

# Build deterministic transforms tailored to each variable
compute_custom_transforms <- function(var_name, values) {
  transforms <- list()
  add_transform <- function(name, vector, label) {
    transforms[[name]] <<- list(
      x.t = vector,
      label = label,
      fill = "#90c987"
    )
  }
  
  var_lower <- stringr::str_to_lower(var_name)
  
  if (var_lower %in% c("prob_occur", "suitability")) {
    add_transform(
      "logit_clipped",
      safe_scale_vec(logit_clipped(values)),
      "Logit (clipped) + z-score"
    )
    add_transform(
      "zscore_manual",
      safe_scale_vec(values),
      "Z-score"
    )
  }
  
  if (var_lower %in% c("centroid_d", "boundary_d")) {
    add_transform(
      "log1p_zscore",
      safe_scale_vec(log1p(pmax(values, 0))),
      "log1p + z-score"
    )
  }
  
  if (var_lower %in% c("marginality", "specificity")) {
    add_transform(
      "zscore_manual",
      safe_scale_vec(values),
      "Z-score"
    )
  }
  
  transforms
}

# Close plotting device if requested (useful for scripted runs)
close_plot_device <- function(auto_close = FALSE) {
  if (isTRUE(auto_close) && dev.cur() > 1) dev.off()
}


# Canonical string for matching transform names (case/punct insensitive)
canonicalize_transform_name <- function(name) {
  stringr::str_replace_all(stringr::str_to_lower(name %||% ""), "[^a-z0-9]", "")
}

# Determine best method name to label plots
get_best_method_name <- function(bn_obj) {
  if (!is.null(bn_obj$norm_stats)) {
    idx <- which.min(bn_obj$norm_stats)
    if (length(idx) > 0) {
      candidate <- names(bn_obj$norm_stats)[idx][1]
      if (!is.na(candidate) && candidate != "") return(candidate)
    }
  }
  
  chosen_class <- class(bn_obj$chosen_transform)
  if (length(chosen_class) > 0) {
    return(chosen_class[1])
  }
  
  method_string <- bn_obj$method
  if (!is.null(method_string) && nzchar(method_string)) {
    return(method_string)
  }
  
  "Unknown"
}

# Retrieve custom transform metadata if present
get_custom_transform <- function(bn_obj, transform_name) {
  custom <- bn_obj$custom_transforms %||% list()
  custom[[transform_name]]
}

get_transform_values <- function(bn_obj, transform_name) {
  if (is.null(transform_name) || length(transform_name) == 0 || is.na(transform_name)) {
    return(NULL)
  }
  
  canonical_name <- canonicalize_transform_name(transform_name)
  best_name <- canonicalize_transform_name(get_best_method_name(bn_obj))
  chosen_classes <- canonicalize_transform_name(class(bn_obj$chosen_transform))
  method_name <- canonicalize_transform_name(bn_obj$method)
  
  chosen_match <- canonical_name %in% c("best", "chosen") ||
    canonical_name == best_name ||
    canonical_name %in% chosen_classes ||
    (nzchar(method_name) && canonical_name == method_name)
  
  if (chosen_match) {
    return(bn_obj$chosen_transform$x.t)
  }
  
  other_names <- names(bn_obj$other_transforms)
  if (length(other_names) == 0) {
    return(NULL)
  }
  
  other_canon <- canonicalize_transform_name(other_names)
  idx <- which(other_canon == canonical_name)
  if (length(idx) > 0) {
    target <- other_names[idx[1]]
    transform_obj <- bn_obj$other_transforms[[target]]
    if (!is.null(transform_obj)) {
      return(transform_obj$x.t)
    }
  }
  
  custom_transform <- get_custom_transform(bn_obj, transform_name)
  if (!is.null(custom_transform)) {
    return(custom_transform$x.t)
  }
  
  NULL
}

friendly_transform_name <- function(name, bn_obj = NULL) {
  if (length(name) == 0 || is.na(name)) {
    return("Unknown")
  }
  
  if (!is.null(bn_obj)) {
    custom_transform <- get_custom_transform(bn_obj, name)
    if (!is.null(custom_transform)) {
      return(custom_transform$label %||% name)
    }
  }
  
  if (name %in% names(transform_labels)) {
    return(transform_labels[[name]])
  }
  
  label <- stringr::str_replace_all(name, "_", " ")
  stringr::str_to_title(label)
}

# Is this transform the best-ranked (for labelling/fill)
is_best_transform <- function(bn_obj, transform_name) {
  canonicalize_transform_name(transform_name) ==
    canonicalize_transform_name(get_best_method_name(bn_obj))
}

# Decide fill colour based on role (best/custom/default)
get_transform_fill <- function(bn_obj, transform_name, default_fill = "#bdbdbd") {
  if (is_best_transform(bn_obj, transform_name)) {
    return("#fb9a99")
  }
  custom_transform <- get_custom_transform(bn_obj, transform_name)
  if (!is.null(custom_transform)) {
    return(custom_transform$fill %||% "#90c987")
  }
  default_fill
}

# Draw grid of histograms for all available transforms (built-in + custom)
plot_transform_grid <- function(bn_obj, var_name,
                                grid_transforms = NULL,
                                layout = NULL,
                                include_original = TRUE,
                                auto_close_plots = FALSE) {
  if (is.null(grid_transforms)) {
    grid_transforms <- c(
      names(bn_obj$norm_stats),
      names(bn_obj$custom_transforms %||% list())
    )
  }
  
  grid_transforms <- unique(stats::na.omit(grid_transforms))
  n_panels <- length(grid_transforms) + as.integer(include_original)
  
  if (is.null(layout)) {
    ncol <- 3
    nrow <- ceiling(n_panels / ncol)
    layout <- c(max(1, nrow), ncol)
  }
  
  par(mfrow = layout, oma = c(0, 0, 3, 0))
  
  if (isTRUE(include_original)) {
    hist_centered(
      bn_obj$x,
      main = "Original data",
      fill = "#a6cee3"
    )
  }
  
  purrr::walk(grid_transforms, ~{
    transform_values <- get_transform_values(bn_obj, .x)
    if (is.null(transform_values)) return(invisible(NULL))
    
    best_suffix <- if (is_best_transform(bn_obj, .x)) " (best)" else ""
    panel_title <- glue::glue("{friendly_transform_name(.x, bn_obj)}{best_suffix}")
    fill_color <- get_transform_fill(bn_obj, .x)
    hist_centered(
      transform_values,
      main = panel_title,
      fill = fill_color
    )
  })
  mtext(glue::glue("{var_name}: transformation comparison"), outer = TRUE, cex = 1.1, font = 2)
  close_plot_device(auto_close_plots)
}

# Paired comparison (original vs selected transform)
plot_transform_pair <- function(bn_obj, var_name, transform_name, auto_close_plots = FALSE) {
  transform_values <- get_transform_values(bn_obj, transform_name)
  if (is.null(transform_values)) return(invisible(NULL))
  
  resolved_name <- transform_name
  if (canonicalize_transform_name(transform_name) %in% c("best", "chosen")) {
    resolved_name <- get_best_method_name(bn_obj)
    transform_name <- resolved_name
  }
  best_suffix <- if (is_best_transform(bn_obj, transform_name)) " (best)" else ""
  
  par(mfrow = c(2, 1))
  hist_centered(bn_obj$x, main = glue::glue("{var_name}: original data"), fill = "#a6cee3")
  hist_centered(
    transform_values,
    main = glue::glue("{var_name}: {friendly_transform_name(resolved_name, bn_obj)} transformed data{best_suffix}"),
    fill = get_transform_fill(bn_obj, resolved_name)
  )
  close_plot_device(auto_close_plots)
}

# Wrapper: fit bestNormalize, add custom transforms, and plot diagnostics
best_normalize_workflow <- function(values,
                                    var_name,
                                    jitter_width = NULL,
                                    jitter_seed = 123,
                                    use_jitter = FALSE,
                                    grid_transforms = NULL,
                                    pairwise_transforms = NULL,
                                    auto_close_plots = FALSE) {
  vector_to_transform <- values
  
  if (isTRUE(use_jitter) && !is.null(jitter_width) && jitter_width > 0) {
    set.seed(jitter_seed)
    vector_to_transform <- vector_to_transform + runif(length(vector_to_transform),
                                                       -jitter_width,
                                                       jitter_width)
  }
  
  bn_obj <- bestNormalize::bestNormalize(vector_to_transform)
  bn_obj$custom_transforms <- compute_custom_transforms(var_name, vector_to_transform)
  
  plot_transform_grid(
    bn_obj = bn_obj,
    var_name = var_name,
    grid_transforms = grid_transforms,
    auto_close_plots = auto_close_plots
  )
  
  if (is.null(pairwise_transforms)) {
    pairwise_transforms <- "best"
  }
  
  purrr::walk(
    pairwise_transforms,
    ~plot_transform_pair(
      bn_obj,
      var_name,
      .x,
      auto_close_plots = auto_close_plots
    )
  )
  
  bn_obj
}

best_normalize_results <- list()


# 4.2 Apply normalization workflows to covariates --------------------------
## prob_occur --------------------------------------------------------------
best_normalize_results$prob_occur <- best_normalize_workflow(
  values = dat_clean_agg$prob_occur,
  var_name = "prob_occur",
  use_jitter = FALSE
)

# Example back-transform workflow: original values → transformed → recovered
prob_occur_bt_demo <- tibble(prob_occur = dat_clean_agg$prob_occur) %>%
  mutate(
    prob_occur_transformed = predict(
      best_normalize_results$prob_occur,
      newdata = prob_occur,
      inverse = FALSE
    ),
    prob_occur_back = predict(
      best_normalize_results$prob_occur,
      newdata = prob_occur_transformed,
      inverse = TRUE
    ),
    diff = prob_occur_back - prob_occur
  )

# Quick check: maximum absolute difference should be ~0 (floating point noise)
prob_occur_bt_demo %>%
  summarise(max_abs_diff = max(abs(diff), na.rm = TRUE)) %>%
  print()

# Show first few rows as a concrete example
prob_occur_bt_demo %>%
  slice_head(n = 6) %>%
  print()

## Marginality -------------------------------------------------------------------
best_normalize_results$Marginality <- best_normalize_workflow(
  values = dat_clean_agg$Marginality,
  var_name = "Marginality",
  use_jitter = FALSE
)

## Specificity ------------------------------------------------------------------
best_normalize_results$Specificity <- best_normalize_workflow(
  values = dat_clean_agg$Specificity,
  var_name = "Specificity",
  use_jitter = FALSE
)

## Suitability ------------------------------------------------------------------
best_normalize_results$Suitability <- best_normalize_workflow(
  values = dat_clean_agg$Suitability,
  var_name = "Suitability",
  use_jitter = FALSE
)

table(dat_clean_agg$Suitability)[order(table(dat_clean_agg$Suitability), decreasing = TRUE)]

## Centroid_d ------------------------------------------------------------------
best_normalize_results$Centroid_d <- best_normalize_workflow(
  values = dat_clean_agg$Centroid_d,
  var_name = "Centroid_d",
  use_jitter = FALSE
)

## Boundary_d ------------------------------------------------------------------
best_normalize_results$Boundary_d <- best_normalize_workflow(
  values = dat_clean_agg$Boundary_d,
  var_name = "Boundary_d",
  use_jitter = FALSE
)

## <> Save best transforms as new columns in dat_clean_agg ------------------------
# Check name of chosen transforms
transform_types = lapply(covariate_vars, function(var) {
  class(best_normalize_results[[var]]$chosen_transform)
})

if(all(transform_types=="orderNorm")){
dat_clean_agg <- dat_clean_agg %>%
  mutate(
    across(
      .cols = c(prob_occur, Marginality, Specificity, Suitability, Centroid_d, Boundary_d),
      .fns = ~ predict(best_normalize_results[[cur_column()]], newdata = .x, inverse = FALSE),
      .names = "{.col}_orderNorm"
    )
  )
}

table(dat_clean_agg$Suitability)[order(table(dat_clean_agg$Suitability), decreasing = TRUE)][1:50]
table(dat_clean_agg$Suitability_orderNorm)[order(table(dat_clean_agg$Suitability_orderNorm), decreasing = TRUE)][1:30]


a = table(dat_clean_agg$host_species)[order(table(dat_clean_agg$host_species), decreasing = TRUE)]

table(dat_clean_agg$prob_occur)[order(table(dat_clean_agg$prob_occur), decreasing = TRUE)][1:30]

# 4.3 Save transformed data and objects -----------------------
saveRDS(dat_clean_agg, "Data/dat_clean_agg2.rds")
saveRDS(best_normalize_results, "Data/best_normalize_results2.rds")
