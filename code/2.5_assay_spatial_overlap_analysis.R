## Assay Group vs. Spatial Location Analysis
## ------------------------------------------------------------
## Goal: Determine if 'assay_group' provides information beyond what the 
## spatial Gaussian process captures by checking for overlapping assay groups
## at unique geographical locations.
##
## Run from project root:
##   Rscript "code/2.5_assay_spatial_overlap_analysis.R"

# 1. Packages -------------------------------------------------------------
pacman::p_load(tidyverse, here, sf, patchwork)

# 2. Load data ------------------------------------------------------------
dat <- readRDS(here("Data", "dat_clean_agg.rds"))

message("Total records in aggregated data: ", nrow(dat))

# 3. Define unique sites and analyze assay distribution -------------------

# We define a site by its unique longitude and latitude.
# In the models, these are used for the Gaussian Process.
sites_assay <- dat |>
  group_by(longitude, latitude) |>
  summarise(
    n_records = n(),
    n_assays = n_distinct(assay_group),
    assays_present = paste(sort(unique(assay_group)), collapse = ", "),
    n_species = n_distinct(host_species),
    species_present = paste(sort(unique(host_species)), collapse = ", "),
    total_tested = sum(number_tested),
    .groups = "drop"
  ) |>
  mutate(has_multiple_assays = n_assays > 1)

# 4. Global Statistics ----------------------------------------------------
n_total_sites <- nrow(sites_assay)
n_overlap_sites <- sum(sites_assay$has_multiple_assays)
pct_overlap_sites <- (n_overlap_sites / n_total_sites) * 100

message("\n--- Site-level Assay Summary ---")
message("Total unique locations: ", n_total_sites)
message("Locations with multiple assay groups: ", n_overlap_sites, " (", round(pct_overlap_sites, 2), "%)")

# How much of our data (tested individuals) comes from these overlapping sites?
overlap_tested <- sites_assay |> 
  filter(has_multiple_assays) |> 
  pull(total_tested) |> 
  sum()
total_tested_all <- sum(sites_assay$total_tested)
pct_overlap_tested <- (overlap_tested / total_tested_all) * 100

message("Proportion of total individuals tested at overlapping sites: ", round(pct_overlap_tested, 2), "%")

# 5. Species-level Analysis -----------------------------------------------
# Check for which species this overlap is most common.
species_overlap <- dat |>
  left_join(sites_assay |> select(longitude, latitude, has_multiple_assays), 
            by = c("longitude", "latitude")) |>
  group_by(host_species) |>
  summarise(
    n_sites = n_distinct(paste(longitude, latitude)),
    n_overlap_sites = n_distinct(paste(longitude[has_multiple_assays], latitude[has_multiple_assays])),
    pct_overlap = (n_overlap_sites / n_sites) * 100,
    n_assays_used = n_distinct(assay_group),
    .groups = "drop"
  ) |>
  arrange(desc(pct_overlap))

message("\n--- Top 10 Species by Site-level Assay Overlap ---")
print(head(species_overlap, 10))

# 6. Visualisations -------------------------------------------------------
# 6.1 Map of Overlap
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

map_plot <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray80") +
  geom_point(data = sites_assay, 
             aes(x = longitude, y = latitude, color = has_multiple_assays, size = total_tested),
             alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
                     name = "Multiple Assays at Site") +
  scale_size_continuous(name = "Total Tested") +
  labs(title = "Spatial Distribution of Assay Group Overlap",
       subtitle = "Red dots indicate sites where both Serology and Culture/Molecular were used") +
  theme_minimal() +
  coord_sf()

# 6.2 Histogram of assays per site
hist_plot <- ggplot(sites_assay, aes(x = factor(n_assays))) +
  geom_bar(fill = "steelblue", color = "black") +
  labs(title = "Number of Assay Groups per Unique Location",
       x = "Number of Assay Groups",
       y = "Count of Locations") +
  theme_minimal()

# Save results
results_dir <- here("Results", "assay_analysis")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

ggsave(file.path(results_dir, "assay_site_overlap_map.png"), map_plot, width = 12, height = 8)
ggsave(file.path(results_dir, "assay_per_site_hist.png"), hist_plot, width = 8, height = 6)
write.csv(sites_assay, file.path(results_dir, "site_assay_details.csv"), row.names = FALSE)
write.csv(species_overlap, file.path(results_dir, "species_assay_overlap.csv"), row.names = FALSE)

# 7. Summary Report -------------------------------------------------------
sink(file.path(results_dir, "assay_overlap_report.txt"))
cat("Assay Group vs. Spatial Location Overlap Report\n")
cat("==============================================\n\n")
cat("Date: ", as.character(Sys.Date()), "\n\n")
cat("Summary Statistics:\n")
cat("------------------\n")
cat("Total unique sites (lon/lat): ", n_total_sites, "\n")
cat("Sites with multiple assay groups: ", n_overlap_sites, " (", round(pct_overlap_sites, 2), "%)\n")
cat("Total individuals tested: ", total_tested_all, "\n")
cat("Individuals tested at overlapping sites: ", overlap_tested, " (", round(pct_overlap_tested, 2), "%)\n\n")

cat("Interpretation Guide:\n")
cat("--------------------\n")
cat("1. If overlap is low (<5%): Assay group is largely redundant with the spatial GP.\n")
cat("   Consider removing assay_group from fixed effects or merging it with the spatial term.\n")
cat("2. If overlap is moderate/high (>15%): Assay group captures within-site methodological\n")
cat("   variation that the GP cannot. Keep as a fixed effect.\n\n")

cat("Top Species with Assay Overlap:\n")
print(head(species_overlap, 15))
sink()

message("\nAnalysis complete. Results saved to: ", results_dir)


