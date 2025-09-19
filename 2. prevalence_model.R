#--------------------------------------#
# Predict zoonotic pathogen prevalence #
#--------------------------------------#

# 0. Script purpose ----
# - Load SDM estimates, & host & pathogen data
# - Clean & format data
# - Formulate Bayesian GLMs to predict prevalence given covariates
# - Visualise outputs


# 1. Load packages ----
pacman::p_load(sf, tidyverse, brms, bayesplot)


# 2. Load data ----
arha <- readRDS("Data/Project_ArHa_database_2025-09-18.rds")  # ArHA
sdm_dat <- read.csv("Data/Data_partial.csv")  # SDM estimates


# 3. Wrangle data ----

## a. ArHA ----

# Get host & pathogen data
arha_host <- arha$host 
arha_path <- arha$pathogen

# Join pathogen to host data
hostpath <- arha_host %>% 
  left_join(arha_path) %>% 
  rename(n_assayed  = number_tested,
         n_positive = number_positive,
         decimalLatitude = latitude,
         decimalLongitude = longitude) %>% 
  # Remove NAs
  filter(if_all(c(n_positive, n_assayed, host_species, pathogen_species_original, assay), ~ !is.na(.))) %>% 
  filter(tolower(assay) != "missing", tolower(coord_status) != "missing") %>% 
  # Add pathogen family name if missing & have pathogen sp name
  mutate(
    pathogen_family = case_when(
      str_detect(tolower(pathogen_species_original), "leptosp")         ~ "Leptospiraceae",
      str_detect(tolower(pathogen_species_original), "rickettsia")      ~ "Rickettsiaceae",
      str_detect(tolower(pathogen_species_original), "bartonella")      ~ "Bartonellaceae",
      str_detect(tolower(pathogen_species_original), "borrelia")        ~ "Spirochaetaceae",
      str_detect(tolower(pathogen_species_original), "poliovirus")      ~ "Picornaviridae",
      str_detect(tolower(pathogen_species_original), "yersinia pestis") ~ "Enterobacteriaceae",
      TRUE ~ pathogen_family
    )
  ) %>% 
  filter(n_assayed != 0) %>% 
  # Group data by temporal resolution
  mutate(
    date_interval = interval(start_date, end_date),
    month_diff = date_interval %/% months(1), 
    month_diff_group = case_when(
      month_diff < 1                    ~ "<1",
      month_diff >= 1 & month_diff < 3  ~ "1-3",
      month_diff >= 3 & month_diff < 6  ~ "3-6",
      month_diff >= 6 & month_diff < 12 ~ "6-12",
      month_diff >= 12                  ~ ">12"),
    month_diff_group = factor(month_diff_group, levels = c("<1", "1-3", "3-6", "6-12", ">12"))
  )

# Filter data for desired spatial & temporal resolution
num_left <- hostpath %>%
  filter(
    # month_diff_group %in% c("<1", "1-3"),
    coordinate_resolution_processed %in% c("site", "village", "town", "city", "adm3")
  ) %>%
  reframe(n = sum(n_assayed, na.rm=T))

num_left/sum(hostpath$n_assayed, na.rm=T)  # proportion of no. tests remaining


## b. SDM ----

# Wrangle SDM data
sdm_clean <- sdm_dat %>% 
  rename(decimalLatitude = bp_Y, decimalLongitude = bp_X, host_species = Species,
         prob_occur = mean) %>% 
  # Remove NAs
  filter(if_all(c(n_positive, n_assayed, host_species, prob_occur, Marginality, 
                  Specificity, Suitability, Centroid_d, Boundary_d), ~ !is.na(.))) %>%
  # Transform to integers
  mutate(n_positive = as.integer(n_positive),
         n_assayed = as.integer(n_assayed)) 


## c. Join ArHA to SDM data ----

# ** NB: need ID columns in SDM data to join properly
# Problematic ones have same lat/long, same n assayed & positive but diff pathogen/host etc [or repeats that don't allow merging the data]

joined_dat <- sdm_clean %>% 
  left_join(hostpath) %>%
  filter(!is.na(assay), assay != "Missing", !is.na(pathogen_family), !is.na(host_family)) %>% 
  mutate(assay = case_when(
    assay == "Serology" | assay == "Western Blot" ~ "Serology", 
    assay != "Serology" ~ "Culture/molecular"))


# 4. Scale covariates ----

# Function to apply unit scaling to covariates (mean of 0, SD of 1)
unitScale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Give covariates unit scaling
scaled_covars <- joined_dat %>% 
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

# 5. Formulate Bayesian GLM ----

# Fit model
mod <- brm(n_positive | trials(n_assayed) ~ 1 + 
            (1|host_family) + (1|pathogen_family) + (1|assay) +
             prob_occur + Marginality + Specificity + Suitability + Centroid_d + Boundary_d,
           family = "binomial"(link = "logit"), data = scaled_covars, 
           prior = c(prior(normal(0, 2), class = "b"),
                     prior(normal(0, 2), class = "Intercept")),
           iter = 2000, chains = 1, seed = 1)


# 6. Plot posteriors ----
niche_metrics <- c("prob_occur", "Marginality", "Specificity", "Suitability", "Centroid_d", "Boundary_d")
posterior <- as.matrix(mod)
plot_title <- ggtitle("Posterior distributions", "with medians and 80% credible intervals")
mcmc_areas(posterior, pars = paste0("b_", niche_metrics), prob = 0.8) + plot_title
