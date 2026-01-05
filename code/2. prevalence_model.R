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


# 1. Load packages -------------------------------------------------
pacman::p_load(sf, tidyverse, brms, bayesplot, bestNormalize, tidybayes, here)

# 2. Load standardised data -------------------------------------------------
data_dir <- here("Data")
dat = readRDS(here("Data", "dat_clean_agg.rds")) # loads 'prevalence_data' dataframe
names(dat)
# r$> names(dat)
#  [1] "host_family"              "pathogen_family"          "assay_group"              "host_species"             "pathogen_species_cleaned" "prob_occur"               "Marginality"             
#  [8] "Specificity"              "Suitability"              "Centroid_d"               "Boundary_d"               "longitude"                "latitude"                 "number_tested"           
# [15] "number_positive"          "number_negative"          "number_inconclusive"      "rows"                     "prob_occur_orderNorm"     "Marginality_orderNorm"    "Specificity_orderNorm"   
# [22] "Suitability_orderNorm"    "Centroid_d_orderNorm"     "Boundary_d_orderNorm"  

scaled_covars <- dat
hist(scaled_covars$Suitability)

# 3. Formulate & fit model --------------------------------------------------

## (a1) Explore covariate distributions ----

# Check distributions of all covariates
# Look at raw values 
scaled_covars %>%
  pivot_longer(cols = all_of(c("prob_occur", "Marginality", "Specificity", "Suitability",
                               "Centroid_d", "Boundary_d"))) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free") +
  labs(title = "Raw covariate distributions")

# Look at orderNorm-transformed values (should be approximately N(0,1))
scaled_covars %>%
  pivot_longer(cols = all_of(c("prob_occur_orderNorm", "Marginality_orderNorm", 
                               "Specificity_orderNorm", "Suitability_orderNorm",
                               "Centroid_d_orderNorm", "Boundary_d_orderNorm"))) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free") +
  labs(title = "orderNorm-transformed covariates (should be ~N(0,1))")


## (a2) Set priors ----
# Since predictors are orderNorm-transformed (Mean = 0, SD = 1), we use:
# - normal(0, 1.5) for intercept (allows wide range of baseline prevalence)
# - normal(0, 1) for fixed effects (weakly regularizing: OR ~ 2.7 per SD)
# - exponential(1) for random effect SDs (standard weakly informative choice)

priors_m0 <- prior(normal(0, 1.5), class = "Intercept") +
  prior(exponential(1), class = "sd")

priors_m1 <- prior(normal(0, 1.5), class = "Intercept") +
  prior(normal(0, 1), class = "b") +
  prior(exponential(1), class = "sd")


## (a2.1) Visualize Priors (New Section) ----

# 1. Visualize coefficients on the Log-Odds scale (what the model sees)
prior_samples <- data.frame(
  log_odds = rnorm(10000, mean = 0, sd = 1) # simulating our normal(0, 1) prior
)

p1 <- ggplot(prior_samples, aes(x = log_odds)) +
  geom_density(fill = "skyblue", alpha = 0.7) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "red") +
  labs(title = "Prior on Log-Odds Scale: normal(0, 1)",
       subtitle = "Red lines: ±2 SD (covering 95% of probability mass)",
       x = "Log-Odds Coefficient (beta)") +
  theme_minimal()

# 2. Visualize coefficients on the Odds Ratio scale (Biological interpretation)
# Exponentiating the log-odds gives the multiplicative effect on odds
prior_samples$odds_ratio <- exp(prior_samples$log_odds)

p2 <- ggplot(prior_samples, aes(x = odds_ratio)) +
  geom_density(fill = "orange", alpha = 0.7) +
  coord_cartesian(xlim = c(0, 10)) + # Zoom in to relevant range
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title = "Prior on Odds Ratio Scale: exp(normal(0, 1))",
       subtitle = "Implies most effects multiply odds by 0.3x to 3x",
       x = "Odds Ratio (OR)") +
  theme_minimal()

# Combine and display
# gridExtra::grid.arrange(p1, p2, ncol = 2) 
print(p1)
print(p2)


## (a3) Prior predictive checks for M0 (Null Model) ----
# Goal: Ensure priors don't predict biologically impossible scenarios
# before the model sees data

qa_dir <- file.path("Results", "QA")
if (!dir.exists(qa_dir)) dir.create(qa_dir, recursive = TRUE)

# Observed prevalence (proportion) for reference
obs_prev <- with(scaled_covars, number_positive / number_tested)
stopifnot(all(is.finite(obs_prev)), all(obs_prev >= 0), all(obs_prev <= 1))

# Fit prior-only M0 model (intercept + random effects only)
m0_prior <- brm(
  number_positive | trials(number_tested) ~ 1 +
    (1 | host_species) + 
    (1 | host_family) + 
    (1 | pathogen_family) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m0,
  sample_prior = "only",
  iter = 2000,
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results/m0_prior")
)

# Draw prior-predictive replicated counts
set.seed(123)
yrep_prior <- brms::posterior_predict(m0_prior, ndraws = 500)

# Convert to prevalence proportions per replicate
yrep_prev <- sweep(yrep_prior, 2, scaled_covars$number_tested, "/")

# Check the prior predictive for a simple statistic: mean prevalence
p_prior_stat <- bayesplot::ppc_stat(
  y = scaled_covars$number_positive,
  yrep = yrep_prior,
  stat = function(y) mean(y / scaled_covars$number_tested)
) +
  labs(title = "Prior predictive check: mean prevalence (M0)")

print(p_prior_stat)

ggsave(filename = file.path(qa_dir, "m0_prior_ppc_mean_prevalence.png"),
       plot = p_prior_stat, width = 7, height = 5, dpi = 300)

# Compare the full distribution of prevalence (overlay density)
p_prior_dens <- bayesplot::ppc_dens_overlay(
  y = obs_prev,
  yrep = yrep_prev[1:50, , drop = FALSE]
) +
  labs(title = "Prior predictive check: prevalence distribution (M0)")

print(p_prior_dens)

ggsave(filename = file.path(qa_dir, "m0_prior_ppc_prevalence_density.png"),
       plot = p_prior_dens, width = 7, height = 5, dpi = 300)

# Inspect expected prevalence under the prior (E[y]/n)
mu_counts <- posterior_epred(m0_prior, ndraws = 500)
mu_prev   <- sweep(mu_counts, 2, scaled_covars$number_tested, "/")

mu_prev_df <- tibble(value = as.vector(mu_prev))

p_mu_prev <- mu_prev_df %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 40, fill = "#fb9a99", colour = "white") +
  labs(
    title = "Expected prevalence under the prior (M0)",
    x = "E[prevalence]",
    y = "Count"
  ) +
  theme_minimal()

ggsave(filename = file.path(qa_dir, "m0_prior_expected_prevalence_hist.png"),
       plot = p_mu_prev, width = 7, height = 5, dpi = 300)

# Print to interactive device if running interactively
print(p_prior_stat)
print(p_prior_dens)
print(p_mu_prev)


## (b) Fit M0: Null Model (intercept + random effects only) ----
# Establishes baseline variance explained by grouping structure

m0 <- brm(
  number_positive | trials(number_tested) ~ 1 +
    (1 | host_species) + 
    (1 | host_family) + 
    (1 | pathogen_family) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m0,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results/m0")
)

# Summarize M0
summary(m0)
# > summary(m0)
#  Family: binomial 
#   Links: mu = logit 
# Formula: number_positive | trials(number_tested) ~ 1 + (1 | host_species) + (1 | host_family) + (1 | pathogen_family) + (1 | assay_group) 
#    Data: scaled_covars (Number of observations: 2538) 
#   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#          total post-warmup draws = 4000

# Multilevel Hyperparameters:
# ~assay_group (Number of levels: 2) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.37      0.60     0.00     2.17 1.01      976     1418

# ~host_family (Number of levels: 3) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.55      0.59     0.01     2.13 1.00     1085     1845

# ~host_species (Number of levels: 29) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.21      0.20     0.88     1.65 1.00      955     1705

# ~pathogen_family (Number of levels: 2) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.79      0.76     0.11     2.93 1.00     2038     2885

# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -2.58      1.00    -3.98    -0.13 1.01     1308     2390

# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


# Posterior predictive check for M0
brms::pp_check(m0, ndraws = 50) + 
  labs(title = "Posterior predictive check: M0 (Null)")

# Full range of posterior predictive checks
pp_check(m0, ndraws = 50, type = "bars")

# Zoom in to the range of interest
pp_check(m0, ndraws = 50, type = "bars", ) +
  coord_cartesian(xlim = c(0, 10))

# Relative prevalence 
y  <- scaled_covars$number_positive
n  <- scaled_covars$number_tested
p_obs <- y / n

yrep     <- posterior_predict(m0, ndraws = 200)
p_rep    <- sweep(yrep, 2, n, "/")

bayesplot::ppc_dens_overlay(p_obs, p_rep[1:50, ]) +
  labs(title = "PPC: prevalence distribution (M0)")

pp_check(m0, ndraws = 200, type = "stat",
         stat = function(y) mean(y / n))

pp_check(m0, ndraws = 200, type = "stat",
         stat = function(y) mean(y == 0))

# Loo pit QQ plots
pp_check(m0,type = "loo_pit_qq", ndraws = 200)

# Explained variance 
bayes_R2(m0)
posterior_summary(m0, pars="sd_")
## (c) Fit M1 Series: Single Predictor Models ----
# Test individual predictive power of each niche metric

# M1a: Probability of occurrence
m1a <- brm(
  number_positive | trials(number_tested) ~ 1 + prob_occur_orderNorm +
    (1 | host_species) + 
    (1 | host_family) + 
    (1 | pathogen_family) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1a"
)

# M1b: Suitability
m1b <- brm(
  number_positive | trials(number_tested) ~ 1 + Suitability_orderNorm +
    (1 | host_species) + 
    (1 | host_family) + 
    (1 | pathogen_family) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1b"
)

# M1c: Marginality
m1c <- brm(
  number_positive | trials(number_tested) ~ 1 + Marginality_orderNorm +
    (1 | host_species) + 
    (1 | host_family) + 
    (1 | pathogen_family) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1c"
)

# M1d: Specificity
m1d <- brm(
  number_positive | trials(number_tested) ~ 1 + Specificity_orderNorm +
    (1 | host_species) + 
    (1 | host_family) + 
    (1 | pathogen_family) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1d"
)

# M1e: Centroid distance
m1e <- brm(
  number_positive | trials(number_tested) ~ 1 + Centroid_d_orderNorm +
    (1 | host_species) + 
    (1 | host_family) + 
    (1 | pathogen_family) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1e"
)

# M1f: Boundary distance
m1f <- brm(
  number_positive | trials(number_tested) ~ 1 + Boundary_d_orderNorm +
    (1 | host_species) + 
    (1 | host_family) + 
    (1 | pathogen_family) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1f"
)


### (c1) Model comparison: M0 vs M1 series ----
# Use LOO-CV to compare out-of-sample predictive accuracy

# Compute LOO for all models
loo_m0  <- loo(m0, moment_match = TRUE)
loo_m1a <- loo(m1a)
loo_m1b <- loo(m1b)
loo_m1c <- loo(m1c)
loo_m1d <- loo(m1d)
loo_m1e <- loo(m1e)
loo_m1f <- loo(m1f)

# Compare all models
loo_compare(loo_m0, loo_m1a, loo_m1b, loo_m1c, loo_m1d, loo_m1e, loo_m1f)

# Store comparison results
loo_comparison <- loo_compare(loo_m0, loo_m1a, loo_m1b, loo_m1c, loo_m1d, loo_m1e, loo_m1f)
print(loo_comparison)


## (c simple) Fit M0 & M1 simple ------------------

### (M0 simple) Fit M0 simple ------------------
m0_simple <- brm(
  number_positive | trials(number_tested) ~ 1 + assay_group +
    (1 | host_species),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m0,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results/m0_simple")
)

m0_alternative <- brm(
  number_positive | trials(number_tested) ~ 1 + assay_group + (1 | host_species) + (1 | host_family) + (1 | pathogen_family),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m0,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results/m0_alternative")
)

loo_m0 = loo(m0, moment_match = TRUE)
loo_m0_simple = loo(m0_simple, moment_match = TRUE)
loo_m0_alternative = loo(m0_alternative, moment_match = TRUE)

loo_compare(loo_m0, loo_m0_simple, loo_m0_alternative)



# Posterior predictive check for M0
brms::pp_check(m0_simple, ndraws = 50) + 
  labs(title = "Posterior predictive check: M0 (Null)")

# Full range of posterior predictive checks
pp_check(m0_simple, ndraws = 50, type = "bars")

# Zoom in to the range of interest
pp_check(m0_simple, ndraws = 50, type = "bars", ) +
  coord_cartesian(xlim = c(0, 10))

# Relative prevalence 
y  <- scaled_covars$number_positive
n  <- scaled_covars$number_tested
p_obs <- y / n

yrep     <- posterior_predict(m0_simple, ndraws = 200)
p_rep    <- sweep(yrep, 2, n, "/")

bayesplot::ppc_dens_overlay(p_obs, p_rep[1:50, ]) +
  labs(title = "PPC: prevalence distribution (M0)")

pp_check(m0_simple, ndraws = 200, type = "stat",
         stat = function(y) mean(y / n))

pp_check(m0_simple, ndraws = 200, type = "stat",
         stat = function(y) mean(y == 0))

# Loo pit QQ plots
pp_check(m0_simple,type = "loo_pit_qq", ndraws = 200)

# Explained variance 
bayes_R2(m0_simple)
posterior_summary(m0_simple, pars="sd_")

### (M1 simple) Fit M1 simple Series: Single Predictor Models -------
# Test individual predictive power of each niche metric

# M1a: Probability of occurrence
m1a_simple <- brm(
  number_positive | trials(number_tested) ~ 1 + prob_occur_orderNorm +
    (1 | host_species) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1a_simple"
)

# M1b: Suitability
m1b_simple <- brm(
  number_positive | trials(number_tested) ~ 1 + Suitability_orderNorm +
    (1 | host_species) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1b_simple"
)

# M1c: Marginality
m1c_simple <- brm(
  number_positive | trials(number_tested) ~ 1 + Marginality_orderNorm +
    (1 | host_species) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1c_simple"
)

# M1d: Specificity
m1d_simple <- brm(
  number_positive | trials(number_tested) ~ 1 + Specificity_orderNorm +
    (1 | host_species) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1d_simple"
)

# M1e: Centroid distance
m1e_simple <- brm(
  number_positive | trials(number_tested) ~ 1 + Centroid_d_orderNorm +
    (1 | host_species) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1e_simple"
)

# M1f: Boundary distance
m1f_simple <- brm(
  number_positive | trials(number_tested) ~ 1 + Boundary_d_orderNorm +
    (1 | host_species) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m1f_simple"
)


### (c1 simple) Model comparison: M0 vs M1 series ----
# Use LOO-CV to compare out-of-sample predictive accuracy

# Compute LOO for all models
# loo_m0 <- loo(m0)
# loo_m0_simple  <- loo(m0_simple)
# loo_m1a_simple <- loo(m1a_simple)
# loo_m1b_simple <- loo(m1b_simple)
# loo_m1c_simple <- loo(m1c_simple)
# loo_m1d_simple <- loo(m1d_simple)
# loo_m1e_simple <- loo(m1e_simple)
# loo_m1f_simple <- loo(m1f_simple)
# Compute LOO for all models (with moment matching to fix unstable k values)
loo_m0         <- loo(m0, moment_match = TRUE)
loo_m0_simple  <- loo(m0_simple, moment_match = TRUE)
loo_m1a_simple <- loo(m1a_simple, moment_match = TRUE)
loo_m1b_simple <- loo(m1b_simple, moment_match = TRUE)
loo_m1c_simple <- loo(m1c_simple, moment_match = TRUE)
loo_m1d_simple <- loo(m1d_simple, moment_match = TRUE)
loo_m1e_simple <- loo(m1e_simple, moment_match = TRUE)
loo_m1f_simple <- loo(m1f_simple, moment_match = TRUE)

# Compare all models
loo_compare(loo_m0, loo_m0_simple, loo_m1a_simple, loo_m1b_simple, loo_m1c_simple, loo_m1d_simple, loo_m1e_simple, loo_m1f_simple)

# Store comparison results
loo_comparison <- loo_compare(loo_m0, loo_m0_simple, loo_m1a_simple, loo_m1b_simple, loo_m1c_simple, loo_m1d_simple, loo_m1e_simple, loo_m1f_simple)
print(loo_comparison)


## (c2) Check colinarity etc. --------
# Before fitting multivariate models, check for correlation between predictors
cor_matrix <- scaled_covars %>%
  select(ends_with("_orderNorm")) %>%
  cor()

print("Correlation matrix of predictors:")
print(round(cor_matrix, 2))

# Visual check
# pairs(scaled_covars %>% select(ends_with("_orderNorm")), 
#       main = "Correlation between transformed predictors")

# Check specifically Suitability vs Marginality
cor_suit_marg <- cor(scaled_covars$Suitability_orderNorm, scaled_covars$Marginality_orderNorm)
print(paste("Correlation Suitability vs Marginality:", round(cor_suit_marg, 3)))


## (d) Fit M2 models ---------------------------------------------------------
# Multivariate models combining key predictors

# M2a: Additive (Suitability + Marginality)
# Tests if both matter independently
m2a <- brm(
  number_positive | trials(number_tested) ~ 1 + Suitability_orderNorm + Marginality_orderNorm +
    (1 | host_species) + 
    (1 | host_family) + 
    (1 | pathogen_family) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1, # Using same priors for fixed effects
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m2a"
)

# M2b: Interaction (Suitability * Marginality)
# Tests if effect of Suitability depends on Marginality
m2b <- brm(
  number_positive | trials(number_tested) ~ 1 + Suitability_orderNorm * Marginality_orderNorm +
    (1 | host_species) + 
    (1 | host_family) + 
    (1 | pathogen_family) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1, 
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m2b"
)

# Compute LOO for M2 models
loo_m2a <- loo(m2a)
loo_m2b <- loo(m2b)

# Compare M2 models against M0 and best M1 (m1c from previous run)
loo_comparison_m2 <- loo_compare(loo_m0, loo_m1c, loo_m2a, loo_m2b)
print(loo_comparison_m2)

## (e) Fit alternative models ---------------------------------------------------------
# Exploring strategies to reduce LOO standard errors

# M2b_simple: Simplified Random Effects
# Hypothesis: Removing weak grouping factors (assay_group, pathogen_family)
# might reduce noise in the LOO estimation.
m2b_simple <- brm(
  number_positive | trials(number_tested) ~ 1 + Suitability_orderNorm * Marginality_orderNorm +
    (1 | host_species) + 
    (1 | host_family),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1, 
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m2b_simple"
)

# M2c: Alternative Hypothesis (Suitability * Boundary Distance)
# Testing if "distance to edge" is a cleaner signal than "Marginality"
m2c <- brm(
  number_positive | trials(number_tested) ~ 1 + Suitability_orderNorm * Boundary_d_orderNorm +
    (1 | host_species) + 
    (1 | host_family) + 
    (1 | pathogen_family) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1, 
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m2c"
)

# M2d: Interaction with minimal random effects
# Testing assay_group + host_species only (dropping host_family & pathogen_family)
m2d <- brm(
  number_positive | trials(number_tested) ~ 1 + Suitability_orderNorm * Marginality_orderNorm +
    (1 | host_species) + 
    (1 | assay_group),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_m1, 
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = "Results/m2d"
)

# Compute LOO for new alternatives
loo_m2b_simple <- loo(m2b_simple)
loo_m2c        <- loo(m2c)
loo_m2d        <- loo(m2d)

# Compare all interaction models + best single predictor + null
loo_comparison_alt <- loo_compare(loo_m0, loo_m1c, loo_m2b, loo_m2b_simple, loo_m2c, loo_m2d)
print(loo_comparison_alt)

# Check for influential observations in the best model so far (m2b)
plot(loo_m2b, label_points = TRUE)



## (f) Visualize Interaction for Best Model (m2b_simple) ----------------
# We need to understand the Suitability * Marginality interaction.
# We will plot the predicted prevalence vs Suitability at 3 levels of Marginality:
# - Low Marginality (Central niche)
# - Mean Marginality
# - High Marginality (Edge niche)

# 1. Define specific values for Marginality_orderNorm to condition on
# Since it's N(0,1), we choose -1 (Low), 0 (Mean), +1 (High)
cond_effects <- conditional_effects(
  m2b_simple, 
  effects = "Suitability_orderNorm:Marginality_orderNorm", 
  int_conditions = list(Marginality_orderNorm = c(-1, 0, 1))#,
  #spaghetti = TRUE
)

# 2. Customize the plot
p_interaction <- plot(cond_effects, plot = FALSE)[[1]] +
  labs(
    title = "Interaction Effect: Suitability * Marginality",
    subtitle = "Does the effect of Suitability depend on how 'marginal' the habitat is?",
    x = "Suitability (Standardized)",
    y = "Predicted Prevalence Probability",
    fill = "Marginality",
    color = "Marginality"
  ) +
  scale_fill_viridis_d(labels = c("-1 SD (Central)", "Mean", "+1 SD (Marginal)"), alpha = 0.2) +
  scale_color_viridis_d(labels = c("-1 SD (Central)", "Mean", "+1 SD (Marginal)")) +
  theme_minimal()

print(p_interaction)

# Save the plot
ggsave(here("Results/Figures/interaction_m2b_simple.png"), p_interaction, width = 8, height = 6)


# 3. Data Density Check (New)
# Plot the raw data points in the Suitability/Marginality space
# Colored by Observed Prevalence to see if the "Stress Peak" has data support.

p_data_density <- scaled_covars %>%
  mutate(obs_prevalence = number_positive / number_tested) %>%
  ggplot(aes(x = Suitability_orderNorm, y = Marginality_orderNorm)) +
  geom_point(aes(color = obs_prevalence, size = number_tested), alpha = 0.7) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0.5, # Adjust based on your data distribution, or use mean(obs_prev)
    name = "Observed\nPrevalence"
  ) +
  # Alternatively for 0-to-Max gradient:
  # scale_color_gradient(low = "blue", high = "red", name = "Observed\nPrevalence") + 
  scale_size_continuous(range = c(1, 5), name = "Sample Size") +
  geom_hline(yintercept = c(-1, 0, 1), linetype = "dashed", color = "gray50") +
  labs(
    title = "Data Density: Suitability vs Marginality",
    subtitle = "Dashed lines show the 'Low', 'Mean', 'High' levels from the interaction plot.\nCheck the top-left corner (Low Suit/High Marg) for data.",
    x = "Suitability (Standardized)",
    y = "Marginality (Standardized)"
  ) +
  theme_minimal()

print(p_data_density)
ggsave(here("Results/Figures/data_density_check.png"), p_data_density, width = 7, height = 6)


# 4. Visualize posteriors & diagnostics ----

## (a) Posterior predictive checks for M1 models ----

# M1a
pp_check(m1c, ndraws = 1000) + 
  labs(title = "Posterior predictive check: M1a (prob_occur)") +
  coord_cartesian(xlim = c(0, 5))

# M1b  
pp_check(m1b, ndraws = 50) + 
  labs(title = "Posterior predictive check: M1b (Suitability)")

# M1c
pp_check(m1c, ndraws = 50) + 
  labs(title = "Posterior predictive check: M1c (Marginality)")

# Compare best-performing model (based on LOO comparison above)
# Replace 'm1a' with the actual best model from loo_compare() results


## (b) Plot coefficient estimates for M1 series ----
# Extract posterior draws for all M1 models

m1_models <- list(
  m1a = m1a,
  m1b = m1b,
  m1c = m1c,
  m1d = m1d,
  m1e = m1e,
  m1f = m1f
)

# Extract fixed effects coefficients (exclude intercept)
m1_coefs <- map_dfr(names(m1_models), function(model_name) {
  mod <- m1_models[[model_name]]
  as_draws_df(mod) %>%
    select(starts_with("b_"), -b_Intercept) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(model = model_name)
})

# Plot posterior distributions
m1_coefs %>%
  ggplot(aes(x = value, y = model, fill = model)) +
  stat_halfeye(.width = c(0.66, 0.95)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Fixed effect coefficients: M1 series",
    subtitle = "Single-predictor models compared to M0 null",
    x = "Coefficient (log-odds scale)",
    y = "Model"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


## (c) Parameter diagnostics for best model ----
# Replace 'm1a' with the best model from LOO comparison

# Trace plots (check convergence)
mcmc_trace(m1c, pars = c("b_Intercept", "b_prob_occur_orderNorm"))

# Parameter correlations
pairs(m1a)

# Full summary
summary(m1a)
