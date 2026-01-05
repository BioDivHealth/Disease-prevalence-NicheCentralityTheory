## Assay Type Balance across Host Species
## ------------------------------------------------------------
## Goal: Assess how much the estimated “assay_group” effect is driven by 
## genuine methodological differences vs. being confounded with particular species.
##
## Run from project root:
##   Rscript "code/2.6_assay_species_balance_analysis.R"

# 1. Packages -------------------------------------------------------------
pacman::p_load(tidyverse, here, patchwork)

# 2. Load data ------------------------------------------------------------
dat <- readRDS(here("Data", "dat_clean_agg.rds"))

message("Total records in aggregated data: ", nrow(dat))

# 3. Species-Assay Distribution -------------------------------------------
species_assay_balance <- dat |>
  group_by(host_species, assay_group) |>
  summarise(
    n_obs = n(),
    total_tested = sum(number_tested),
    .groups = "drop"
  ) |>
  group_by(host_species) |>
  mutate(
    total_species_obs = sum(n_obs),
    prop_obs = n_obs / total_species_obs,
    n_assay_types = n_distinct(assay_group)
  ) |>
  ungroup()

# Pivot to wide format for easier flagging
balance_wide <- species_assay_balance |>
  select(host_species, assay_group, n_obs) |>
  pivot_wider(names_from = assay_group, values_from = n_obs, values_fill = 0) |>
  rename(
    culture_molecular = `Culture/molecular`,
    serology = `Serology`
  ) |>
  mutate(
    total_obs = culture_molecular + serology,
    culture_prop = culture_molecular / total_obs,
    serology_prop = serology / total_obs,
    
    # Flags
    is_single_assay = (culture_molecular == 0 | serology == 0),
    is_strongly_imbalanced = !is_single_assay & (culture_molecular < 10 | serology < 10),
    status = case_when(
      is_single_assay ~ "Single Assay",
      is_strongly_imbalanced ~ "Strongly Imbalanced (<10)",
      TRUE ~ "Balanced"
    )
  ) |>
  arrange(desc(total_obs))

# 4. Summary Statistics ---------------------------------------------------
n_total_species <- nrow(balance_wide)
n_single_assay <- sum(balance_wide$is_single_assay)
n_imbalanced <- sum(balance_wide$is_strongly_imbalanced)
n_balanced <- sum(balance_wide$status == "Balanced")

message("\n--- Species-level Assay Balance Summary ---")
message("Total species in dataset: ", n_total_species)
message("Balanced species: ", n_balanced)
message("Single-assay species: ", n_single_assay)
message("Strongly imbalanced species (<10 obs): ", n_imbalanced)

# 5. Visualisation --------------------------------------------------------
plot_data <- species_assay_balance |>
  left_join(balance_wide |> select(host_species, status), by = "host_species") |>
  mutate(host_species = reorder(host_species, total_species_obs))

balance_plot <- ggplot(plot_data, aes(y = host_species, x = prop_obs, fill = assay_group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = n_obs), position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_fill_manual(values = c("Culture/molecular" = "#1b9e77", "Serology" = "#d95f02")) +
  labs(title = "Assay Type Distribution per Host Species",
       subtitle = "Labels indicate number of observations (rows) per assay type",
       x = "Proportion of Observations",
       y = "Host Species",
       fill = "Assay Group") +
  theme_minimal() +
  facet_grid(status ~ ., scales = "free_y", space = "free_y")

# 6. Save results ---------------------------------------------------------
results_dir <- here("Results", "assay_analysis")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

ggsave(file.path(results_dir, "assay_species_balance_plot.png"), balance_plot, width = 10, height = 12)
write.csv(balance_wide, file.path(results_dir, "species_assay_balance_details.csv"), row.names = FALSE)

# 7. Summary Report -------------------------------------------------------
sink(file.path(results_dir, "assay_species_balance_report.txt"))
cat("Assay Type Balance across Host Species Report\n")
cat("=============================================\n\n")
cat("Date: ", as.character(Sys.Date()), "\n\n")

cat("Summary Statistics:\n")
cat("------------------\n")
cat("Total species analyzed: ", n_total_species, "\n")
cat("  - Balanced (>=10 obs per type): ", n_balanced, "\n")
cat("  - Strongly Imbalanced (<10 obs): ", n_imbalanced, "\n")
cat("  - Single Assay Only: ", n_single_assay, "\n\n")

cat("Species List and Status:\n")
cat("----------------------\n")
balance_wide |> 
  select(host_species, culture_molecular, serology, status) |>
  print(n = Inf)

cat("\nInterpretation and Robustness Note:\n")
cat("----------------------------------\n")
if (n_balanced > 0) {
  cat(n_balanced, " species provide direct within-species contrasts between assays.\n")
} else {
  cat("WARNING: No species provide balanced contrasts between both assay types.\n")
}

cat("\nRobustness of Assay-Effect Estimate:\n")
if (n_balanced / n_total_species > 0.5) {
  cat("HIGH: More than 50% of species have balanced assay representation. The assay effect\n")
  cat("is likely robust and reflects methodological differences.\n")
} else if (n_balanced > 0) {
  cat("MODERATE: Only ", round(n_balanced/n_total_species*100, 1), "% of species are balanced. While a contrast exists,\n")
  cat("the global assay effect is heavily influenced by a subset of species. It should be\n")
  cat("treated primarily as a control variable.\n")
} else {
  cat("LOW: The assay effect is almost entirely confounded with species identity.\n")
  cat("Interpretation of this coefficient as 'methodological sensitivity' should be very cautious.\n")
}

sink()

message("\nAnalysis complete. Results saved to: ", results_dir)


