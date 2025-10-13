#--------------------------------------#
# Predict zoonotic pathogen prevalence #
#--------------------------------------#

# 0. Script purpose ----
# - Load SDM estimates, host & pathogen data
# - Clean & format data
# - Formulate Bayesian GLMs to predict prevalence given covariates
# - Visualise outputs


# 1. Load packages ----
pacman::p_load(sf, tidyverse, brms, bayesplot)


# 2. Load data ----
arha      <- readRDS("Data/Project_ArHa_database_2025-09-18.rds")  # ArHA
arha_path <- arha$pathogen
arha_host <- arha$host
dat       <- read.csv("Data/Full_data.csv") # SDM estimates + ArHA data


# 3. Wrangle data ----

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
    month_diff_group = case_when(
      month_diff < 1                    ~ "<1",
      month_diff >= 1 & month_diff < 3  ~ "1-3",
      month_diff >= 3 & month_diff < 6  ~ "3-6",
      month_diff >= 6 & month_diff < 12 ~ "6-12",
      month_diff >= 12                  ~ ">12"),
    month_diff_group = factor(month_diff_group, levels = c("<1", "1-3", "3-6", "6-12", ">12"))
  ) %>% 
  reframe(number_positive, number_tested, host_family, pathogen_family, 
          assay_group, prob_occur, Marginality, Specificity, Suitability, 
          Centroid_d, Boundary_d) %>% 
  # Exclude rows where its group has fewer than 20 observations
  group_by(host_family, pathogen_family, assay_group) %>% 
  filter(n() >= 20) %>% 
  ungroup()

# sum(dat_clean$number_tested)
# [1] 83877


# 4. Scale covariates ----

# Function to apply unit scaling to covariates (mean of 0, SD of 1)
unitScale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Give covariates unit scaling
scaled_covars <- dat_clean %>% 
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


# 5. Formulate & fit model ----

## (a) Find suitable priors for all predictors ----

# # Check distributions of all covariates
# scaled_covars %>% 
#   pivot_longer(cols = all_of(c("prob_occur", "Marginality", "Specificity", "Suitability", 
#                                "Centroid_d", "Boundary_d"))) %>%
#   ggplot(aes(value)) + 
#   geom_histogram() + 
#   facet_wrap(~name, scales = "free")

# Set priors - ** ideally would want to change these as poor priors for some predictors **
my_prior <- prior(normal(0, 2), class = b) +
  prior(normal(0, 2), class = "Intercept")

# Validate (& check what defaults are used)
validate_prior(prior = my_prior, 
               formula = number_positive | trials(number_tested) ~ (1|pathogen_family) + (1|assay_group) + prob_occur + Marginality + Specificity + Suitability,
               data = scaled_covars, family = binomial())

## Prior predictive checks (check that priors roughly approximate distribution of data)
# ** to do here **


## (b) Fit Bayesian binomial GLM ----

# NB could take a long time to fit if lots of data / complex model structure
mod <- brm(number_positive | trials(number_tested) ~  # Formulated to properly capture error structure
             (1|host_family) + (1|pathogen_family) + (1|assay_group) + # random intercepts
             prob_occur + Marginality + Specificity + Suitability, # Other predictors
           family = "binomial"(link = "logit"), # Specify binomial model
           data = scaled_covars, # Data
           prior = my_prior,     # Prior distributions
           iter = 2000, # No. MCMC iterations (2000 usually fine)
           chains = 4,  # No. MCMC chains (4 is good)
           cores = 4, , # No. cores to run parallel chains (set to number of chains)
           seed = 1)    # Ensure repeatable sampling by MCMC


## (c) Assess model performance ----

# Approximate leave-one-out cross-validation on fitted model
loo_mod1 <- loo(mod)
# loo_mod2 <- loo(mo2)  # could define alternative models above and compare here
# 
# # Compare different models
# loo_compare(loo_mod1, loo_mod2)


# 6. Plot posteriors ----

# Parameter correlations
pairs(mod)

# Niche metrics
niche_metrics <- c("prob_occur", "Marginality", "Specificity", "Suitability", "Centroid_d", "Boundary_d")
posterior <- as.matrix(mod)
plot_title <- ggtitle("Posterior distributions", "with medians and 80% credible intervals")
mcmc_areas(posterior, pars = paste0("b_", niche_metrics), prob = 0.8) + plot_title

# Host variables
niche_metrics <- c("host_family", "pathogen_family", "assay_group")
posterior <- as.matrix(mod)
plot_title <- ggtitle("Posterior distributions", "with medians and 80% credible intervals")
mcmc_areas(posterior, pars = paste0("b_", niche_metrics), prob = 0.8) + plot_title

