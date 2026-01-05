## Testing models: baseline options for prevalence models
## -----------------------------------------------------
## Goal: define a few stable "baseline" models to use as the foundation for
## M1/M2 niche-predictor models, and to troubleshoot LOO stability.

# 1. Load packages -------------------------------------------------
pacman::p_load(sf, tidyverse, brms, bayesplot, bestNormalize, tidybayes, here)

# 2. Load standardised data ----------------------------------------
dat <- readRDS(here("Data", "dat_clean_agg.rds"))
scaled_covars <- dat |>
  mutate(
    assay_group = as.factor(assay_group),
    host_family = as.factor(host_family),
    pathogen_family = as.factor(pathogen_family),
    host_species = as.factor(host_species),
    obs_id = row_number()
  )

# 3. Priors ---------------------------------------------------------
# Baselines almost always include at least one fixed effect (e.g., assay_group),
# so we define a "baseline prior" that covers intercept, fixed effects, and RE SDs.
priors_base <- prior(normal(0, 1.5), class = "Intercept") +
  prior(normal(0, 1), class = "b") +
  prior(exponential(1), class = "sd")

# If you use beta-binomial, add a weakly informative prior for overdispersion.
# (brms parameterization uses 'phi' for beta_binomial2)
priors_beta_binom <- priors_base +
  prior(exponential(1), class = "phi")

# 4. Baseline model candidates -------------------------------------
# Important design principle from our diagnostics:
# - Avoid random effects with very few levels (e.g., assay_group with 2 levels).
#   Use fixed effects for those instead.

# B0a: Minimal robust baseline (controls + host species heterogeneity)
m0_base <- brm(
  number_positive | trials(number_tested) ~ 1 + assay_group + (1 | host_species),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_base,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v1","m0_base")
)

# B0a_bb: Beta-binomial version of minimal baseline
# (same structure as m0_base but with beta-binomial family)
m0_base_bb <- brm(
  number_positive | trials(number_tested) ~ 1 + assay_group + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v1","m0_base_bb")
)

# B0b: Add "broad biology" controls as fixed effects
# (these had too-few levels to behave well as random effects)
m0_controls <- brm(
  number_positive | trials(number_tested) ~ 1 + assay_group + host_family +
    pathogen_family + (1 | host_species),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_base,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v1","m0_controls")
)

# B0c: Observation-level random intercept (OLRE) for extra-binomial variation
# Use when PPC/LOO suggests overdispersion or influential observations.
m0_olre <- brm(
  number_positive | trials(number_tested) ~ 1 + assay_group + host_family +
    pathogen_family + (1 | host_species) + (1 | obs_id),
  family = binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_base,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v1","m0_olre")
)

# B0d: Beta-binomial baseline (explicit overdispersion)
m0_beta_binom <- brm(
  number_positive | trials(number_tested) ~ 1 + assay_group + host_family +
    pathogen_family + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v1","m0_beta_binom")
)

# 5. Model comparison helpers --------------------------------------
# LOO can be unstable when there are many high Pareto-k observations.
# - Use moment_match when possible.
# - If you still see many k > 0.7, prefer k-fold CV for the short list.

loo_m0_base <- loo(m0_base, moment_match = TRUE)
loo_m0_base_bb <- loo(m0_base_bb, moment_match = TRUE)
loo_m0_controls <- loo(m0_controls, moment_match = TRUE)
#loo_m0_olre <- loo(m0_olre, moment_match = TRUE)
loo_m0_beta_binom <- loo(m0_beta_binom, moment_match = TRUE)

# Compare baseline models
# This helps isolate:
# - m0_base vs m0_base_bb → effect of beta-binomial likelihood on minimal structure
# - m0_base_bb vs m0_beta_binom → effect of adding host_family + pathogen_family controls
# - m0_controls vs m0_beta_binom → effect of beta-binomial likelihood with controls
loo_compare(loo_m0_base, loo_m0_base_bb, loo_m0_controls, loo_m0_beta_binom)

#-------------------------------------------------------------------------------------/

# Common objects
y  <- scaled_covars$number_positive
n  <- scaled_covars$number_tested
p_obs <- y / n

# Posterior predictive draws (adjust ndraws if you want more/less)
yrep_ctrl <- posterior_predict(m0_controls, ndraws = 200)
yrep_bb   <- posterior_predict(m0_beta_binom, ndraws = 200)
yrep_base_bb   <- posterior_predict(m0_base_bb, ndraws = 200)

p_rep_ctrl <- sweep(yrep_ctrl, 2, n, "/")
p_rep_bb   <- sweep(yrep_bb,   2, n, "/")
p_rep_base_bb   <- sweep(yrep_base_bb,   2, n, "/")

# 1) Overall prevalence distribution
bayesplot::ppc_dens_overlay(p_obs, p_rep_ctrl[1:50, ]) +
  ggtitle("Binomial (m0_controls): prevalence distribution")

bayesplot::ppc_dens_overlay(p_obs, p_rep_bb[1:50, ]) +
  ggtitle("Beta-binomial (m0_beta_binom): prevalence distribution")

bayesplot::ppc_dens_overlay(p_obs, p_rep_base_bb[1:50, ]) +
  ggtitle("Beta-binomial (m0_base_bb): prevalence distribution") 

# 1a) Zoom into low prevalence region (0 to 0.2)
bayesplot::ppc_dens_overlay(p_obs, p_rep_ctrl[1:50, ]) +
  ggtitle("Binomial: low prevalence (0-0.2)") +
  coord_cartesian(xlim = c(0, 0.025))

bayesplot::ppc_dens_overlay(p_obs, p_rep_bb[1:50, ]) +
  ggtitle("Beta-binomial: low prevalence (0-0.2)") +
  coord_cartesian(xlim = c(0, 0.025))

bayesplot::ppc_dens_overlay(p_obs, p_rep_base_bb[1:50, ]) +
  ggtitle("Beta-binomial (m0_base_bb): prevalence distribution") +
  coord_cartesian(xlim = c(0, 0.025))

# 1b) Zoom into high prevalence region (0.8 to 1.0)
bayesplot::ppc_dens_overlay(p_obs, p_rep_ctrl[1:50, ]) +
  ggtitle("Binomial: high prevalence (0.8-1.0)") +
  coord_cartesian(xlim = c(0.999, 1.0))

bayesplot::ppc_dens_overlay(p_obs, p_rep_bb[1:50, ]) +
  ggtitle("Beta-binomial: high prevalence (0.8-1.0)") +
  coord_cartesian(xlim = c(0.999, 1.0))

bayesplot::ppc_dens_overlay(p_obs, p_rep_base_bb[1:50, ]) +
  ggtitle("Beta-binomial (m0_base_bb): high prevalence (0.8-1.0)") +
  coord_cartesian(xlim = c(0.999, 1.0))

# 1c) Zoom into middle range (0.3 to 0.7) - adjust as needed
bayesplot::ppc_dens_overlay(p_obs, p_rep_ctrl[1:50, ]) +
  ggtitle("Binomial: middle prevalence (0.3-0.7)") +
  coord_cartesian(xlim = c(0.4975, 0.5025))

bayesplot::ppc_dens_overlay(p_obs, p_rep_bb[1:50, ]) +
  ggtitle("Beta-binomial: middle prevalence (0.3-0.7)") +
  coord_cartesian(xlim = c(0.4975, 0.5025))

bayesplot::ppc_dens_overlay(p_obs, p_rep_base_bb[1:50, ]) +
  ggtitle("Beta-binomial (m0_base_bb): middle prevalence (0.3-0.7)") +
  coord_cartesian(xlim = c(0.4975, 0.5025))

# 2) Mean prevalence and zero-mass stats
pp_check(m0_controls, ndraws = 200, type = "stat",
         stat = function(y) mean(y / n)) +
  ggtitle("Mean prevalence: binomial")

pp_check(m0_beta_binom, ndraws = 200, type = "stat",
         stat = function(y) mean(y / n)) +
  ggtitle("Mean prevalence: beta-binomial")

pp_check(m0_base_bb, ndraws = 200, type = "stat",
         stat = function(y) mean(y / n)) +
  ggtitle("Mean prevalence: beta-binomial (m0_base_bb)")

pp_check(m0_controls, ndraws = 200, type = "stat",
         stat = function(y) mean(y == 0)) +
  ggtitle("Zero mass: binomial")

pp_check(m0_beta_binom, ndraws = 200, type = "stat",
         stat = function(y) mean(y == 0)) +
  ggtitle("Zero mass: beta-binomial")

pp_check(m0_base_bb, ndraws = 200, type = "stat",
         stat = function(y) mean(y == 0)) +
  ggtitle("Zero mass: beta-binomial (m0_base_bb)")

# 3) LOO-PIT calibration
pp_check(m0_controls, type = "loo_pit_qq")
pp_check(m0_beta_binom, type = "loo_pit_qq")
pp_check(m0_base_bb, type = "loo_pit_qq")

# 4) Grouped fit (e.g., assay)
pp_check(m0_controls, type = "stat_grouped",
         stat = "mean", group = "assay_group") +
  ggtitle("Grouped mean prevalence by assay: binomial")

pp_check(m0_beta_binom, type = "stat_grouped",
         stat = "mean", group = "assay_group") +
  ggtitle("Grouped mean prevalence by assay: beta-binomial")

pp_check(m0_base_bb, type = "stat_grouped",
         stat = "mean", group = "assay_group") +
  ggtitle("Grouped mean prevalence by assay: beta-binomial (m0_base_bb)")

# 6. M1 Series: Single Niche Predictors ----------------------------
# Adding one niche predictor at a time to the robust baseline (m0_beta_binom).
# Structure: ~ 1 + [NicheVar] + assay_group + host_family + pathogen_family + (1|host_species)
# Family: beta_binomial

# M1prob: Probability of Occurrence
m1prob <- brm(
  number_positive | trials(number_tested) ~ 1 + prob_occur_orderNorm + 
    assay_group + host_family + pathogen_family + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v1","m1prob")
)

# M1suit: Habitat Suitability
m1suit <- brm(
  number_positive | trials(number_tested) ~ 1 + Suitability_orderNorm + 
    assay_group + host_family + pathogen_family + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v1","m1suit")
)

# M1marg: Marginality
m1marg <- brm(
  number_positive | trials(number_tested) ~ 1 + Marginality_orderNorm + 
    assay_group + host_family + pathogen_family + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v1","m1marg")
)

# M1spec: Specificity
m1spec <- brm(
  number_positive | trials(number_tested) ~ 1 + Specificity_orderNorm + 
    assay_group + host_family + pathogen_family + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v1","m1spec")
)

# M1centr: Centroid Distance
m1centr <- brm(
  number_positive | trials(number_tested) ~ 1 + Centroid_d_orderNorm + 
    assay_group + host_family + pathogen_family + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v1","m1centr")
)

# M1bound: Boundary Distance
m1bound <- brm(
  number_positive | trials(number_tested) ~ 1 + Boundary_d_orderNorm + 
    assay_group + host_family + pathogen_family + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v1","m1bound")
)

# 7. Comparison: M0 vs M1 ------------------------------------------
loo_m1prob <-loo(m1prob) #all k good
loo_m1suit <-loo(m1suit) #all k good
loo_m1marg <-loo(m1marg) #all k good
loo_m1spec <-loo(m1spec) #all k good
loo_m1centr <- loo(m1centr) #all k good
loo_m1bound <- loo(m1bound) #all k good

# loo_m1prob <- loo(m1prob, moment_match = TRUE)
# loo_m1suit <- loo(m1suit, moment_match = TRUE)
# loo_m1marg <- loo(m1marg, moment_match = TRUE)
# loo_m1spec <- loo(m1spec, moment_match = TRUE)
# loo_m1centr <- loo(m1centr, moment_match = TRUE)
# loo_m1bound <- loo(m1bound, moment_match = TRUE)

# Compare all M1 models against the best baseline (m0_beta_binom)
loo_compare(
  loo_m0_beta_binom, 
  loo_m1prob, loo_m1suit, loo_m1marg, loo_m1spec, loo_m1centr, loo_m1bound
)

# 8. M1 Series v2: Single Niche Predictors (Minimal Baseline) -------
# Adding niche predictors to the minimal beta-binomial baseline (m0_base_bb).
# Structure: ~ 1 + [NicheVar] + assay_group + (1 | host_species)
# Family: beta_binomial

# Create directory for v2 models if it doesn't exist
if (!dir.exists(here("Results", "models_v2"))) {
  dir.create(here("Results", "models_v2"), recursive = TRUE)
}

# M1prob_v2: Probability of Occurrence
m1prob_v2 <- brm(
  number_positive | trials(number_tested) ~ 1 + prob_occur_orderNorm + 
    assay_group + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v2","m1prob_v2")
)

# M1suit_v2: Habitat Suitability
m1suit_v2 <- brm(
  number_positive | trials(number_tested) ~ 1 + Suitability_orderNorm + 
    assay_group + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v2","m1suit_v2")
)

# M1marg_v2: Marginality
m1marg_v2 <- brm(
  number_positive | trials(number_tested) ~ 1 + Marginality_orderNorm + 
    assay_group + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v2","m1marg_v2")
)

# M1spec_v2: Specificity
m1spec_v2 <- brm(
  number_positive | trials(number_tested) ~ 1 + Specificity_orderNorm + 
    assay_group + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v2","m1spec_v2")
)

# M1centr_v2: Centroid Distance
m1centr_v2 <- brm(
  number_positive | trials(number_tested) ~ 1 + Centroid_d_orderNorm + 
    assay_group + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v2","m1centr_v2")
)

# M1bound_v2: Boundary Distance
m1bound_v2 <- brm(
  number_positive | trials(number_tested) ~ 1 + Boundary_d_orderNorm + 
    assay_group + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v2","m1bound_v2")
)

# 9. Comparison: M0 vs M1 (v2) -------------------------------------
loo_m1prob_v2 <- loo(m1prob_v2)
loo_m1suit_v2 <- loo(m1suit_v2)
loo_m1marg_v2 <- loo(m1marg_v2)
loo_m1spec_v2 <- loo(m1spec_v2)
loo_m1centr_v2 <- loo(m1centr_v2)
loo_m1bound_v2 <- loo(m1bound_v2)

print(loo_m1prob_v2) # all pareto k good
print(loo_m1suit_v2) # all pareto k good
print(loo_m1marg_v2) # all pareto k good
print(loo_m1spec_v2) # all pareto k good
print(loo_m1centr_v2)# all pareto k good
print(loo_m1bound_v2)# all pareto k good

# loo_m1prob_v2 <- loo(m1prob_v2, moment_match = TRUE)
# loo_m1suit_v2 <- loo(m1suit_v2, moment_match = TRUE)
# loo_m1marg_v2 <- loo(m1marg_v2, moment_match = TRUE)
# loo_m1spec_v2 <- loo(m1spec_v2, moment_match = TRUE)
# loo_m1centr_v2 <- loo(m1centr_v2, moment_match = TRUE)
# loo_m1bound_v2 <- loo(m1bound_v2, moment_match = TRUE)

# Compare all M1_v2 models against the best minimal baseline (m0_base_bb)
loo_compare(
  loo_m0_base_bb, 
  loo_m1prob_v2, loo_m1suit_v2, loo_m1marg_v2, 
  loo_m1spec_v2, loo_m1centr_v2, loo_m1bound_v2
)

# 10. M2 Series v2: Multivariate Niche Predictors (Minimal Baseline) -----
# Testing combinations of niche predictors (additive and interactions)
# Based on M1_v2 results, we focus on the top performers:
# - prob_occur (best single predictor)
# - Marginality (close second)
# - Suitability (classic hypothesis from literature)
# Structure: ~ 1 + [NicheVar combination] + assay_group + (1 | host_species)
# Family: beta_binomial

# M2a_v2: Additive (prob_occur + Marginality)
# Tests if both best M1 predictors independently improve the model
m2a_v2 <- brm(
  number_positive | trials(number_tested) ~ 1 + prob_occur_orderNorm + 
    Marginality_orderNorm + assay_group + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v2","m2a_v2")
)

# M2b_v2: Interaction (prob_occur * Marginality)
# Tests if the effect of occurrence probability depends on marginality
m2b_v2 <- brm(
  number_positive | trials(number_tested) ~ 1 + prob_occur_orderNorm * 
    Marginality_orderNorm + assay_group + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v2","m2b_v2")
)

# M2c_v2: Additive (Suitability + Marginality)
# Classic niche centrality hypothesis: both habitat quality and position matter
m2c_v2 <- brm(
  number_positive | trials(number_tested) ~ 1 + Suitability_orderNorm + 
    Marginality_orderNorm + assay_group + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v2","m2c_v2")
)

# M2d_v2: Interaction (Suitability * Marginality)
# Tests if effect of suitability depends on marginality (niche center hypothesis)
# This was the "winner" in earlier exploratory analysis (m2b_simple)
m2d_v2 <- brm(
  number_positive | trials(number_tested) ~ 1 + Suitability_orderNorm * 
    Marginality_orderNorm + assay_group + (1 | host_species),
  family = beta_binomial(link = "logit"),
  data = scaled_covars,
  prior = priors_beta_binom,
  iter = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 123,
  file = here("Results","models_v2","m2d_v2")
)

# 11. Comparison: M0 vs M1 vs M2 (v2) ----------------------------------
loo_m2a_v2 <- loo(m2a_v2)
loo_m2b_v2 <- loo(m2b_v2)
loo_m2c_v2 <- loo(m2c_v2)
loo_m2d_v2 <- loo(m2d_v2)

print(loo_m2a_v2) 
print(loo_m2b_v2) 
print(loo_m2c_v2) 
print(loo_m2d_v2) 

# Compare all M2_v2 models against baseline and best M1
loo_compare(
  loo_m0_base_bb, 
  loo_m1prob_v2,     # Best single predictor
  loo_m1marg_v2,     # Second best single predictor
  loo_m2a_v2,        # prob_occur + Marginality (additive)
  loo_m2b_v2,        # prob_occur * Marginality (interaction)
  loo_m2c_v2,        # Suitability + Marginality (additive, classic)
  loo_m2d_v2         # Suitability * Marginality (interaction, classic)
)

# 12. Visualize Best M2 Model (if interactions are supported) ------------
# Only run this section if one of the interaction models performs well

# Check if m2d_v2 (Suitability * Marginality) is among top models
# This replicates the approach from code/2. prevalence_model.R

# Conditional effects plot for interaction
cond_effects_m2d <- conditional_effects(
  m2d_v2, 
  effects = "Suitability_orderNorm:Marginality_orderNorm", 
  int_conditions = list(Marginality_orderNorm = c(-1, 0, 1))
)

p_interaction_v2 <- plot(cond_effects_m2d, plot = FALSE)[[1]] +
  labs(
    title = "Interaction: Suitability × Marginality (Minimal Model)",
    subtitle = "Beta-binomial with minimal baseline (assay + host_species RE)",
    x = "Suitability (Standardized)",
    y = "Predicted Prevalence Probability",
    fill = "Marginality",
    color = "Marginality"
  ) +
  scale_fill_viridis_d(labels = c("-1 SD (Central)", "Mean", "+1 SD (Marginal)"), alpha = 0.2) +
  scale_color_viridis_d(labels = c("-1 SD (Central)", "Mean", "+1 SD (Marginal)")) +
  theme_minimal()

print(p_interaction_v2)

# Save plot
ggsave(
  here("Results/Figures/interaction_m2d_v2.png"), 
  p_interaction_v2, 
  width = 8, 
  height = 6, 
  dpi = 300
)

# Data density check: where do we actually have data?
p_data_density_v2 <- scaled_covars %>%
  mutate(obs_prevalence = number_positive / number_tested) %>%
  ggplot(aes(x = Suitability_orderNorm, y = Marginality_orderNorm)) +
  geom_point(aes(color = obs_prevalence, size = number_tested), alpha = 0.7) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0.5,
    name = "Observed\nPrevalence"
  ) +
  scale_size_continuous(range = c(1, 5), name = "Sample Size") +
  geom_hline(yintercept = c(-1, 0, 1), linetype = "dashed", color = "gray50") +
  labs(
    title = "Data Density: Suitability vs Marginality",
    subtitle = "Check for data support in predicted interaction regions",
    x = "Suitability (Standardized)",
    y = "Marginality (Standardized)"
  ) +
  theme_minimal()

print(p_data_density_v2)
ggsave(
  here("Results/Figures/data_density_v2.png"), 
  p_data_density_v2, 
  width = 7, 
  height = 6, 
  dpi = 300
)

# Extract and visualize coefficients for M2 models
fixef(m2a_v2)
fixef(m2b_v2)
fixef(m2c_v2)
fixef(m2d_v2)
