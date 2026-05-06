# 2.1c_spatial_coverage
# Goal: explore spatial sampling coverage quality across hosts and host–pathogen pairs.
# Scope: compute/export coverage metrics

library(pacman)
p_load(tidyverse, here)

# Parameters ---------------------------------------------------------------
site_mode <- "exact" # "exact" or "rounded"
exact_round_digits <- 3
coord_round_digits <- 2 # ~1.1 km latitude

mode_tag <- if (site_mode == "exact") "exact" else paste0("rounded", coord_round_digits)
analysis_dir <- here("Results", "analysis_metadata", mode_tag)

# 1. Load data -------------------------------------------------------------
path_dat_temporal <- here("Data", paste0("dat_with_temporal_metadata_", mode_tag, ".rds"))
path_dat_fallback <- here("Data", "dat_clean_agg2.rds")

dat <- if (file.exists(path_dat_temporal)) {
  readRDS(path_dat_temporal)
} else {
  readRDS(path_dat_fallback)
} %>%
  select(-any_of("date_interval"))

# Define sites using rounded coordinates to avoid artificial inflation of n_sites
# from small coordinate jitter across studies.
use_site_name <- site_mode != "exact" &&
  all(c("site_name", "lon_round", "lat_round") %in% names(dat))

dat_sites <- if (use_site_name) {
  dat %>%
    mutate(
      lon_site = coalesce(lon_round, round(longitude, coord_round_digits)),
      lat_site = coalesce(lat_round, round(latitude, coord_round_digits)),
      site_id = if_else(!is.na(site_name) & nzchar(site_name), site_name,
                        paste(lon_site, lat_site, sep = "_"))
    )
} else {
  dat %>%
    mutate(
      lon_site = round(longitude, if_else(site_mode == "exact", exact_round_digits, coord_round_digits)),
      lat_site = round(latitude, if_else(site_mode == "exact", exact_round_digits, coord_round_digits)),
      site_id = paste(lon_site, lat_site, sep = "_")
    )
}

site_lut <- dat_sites %>%
  distinct(site_id, lon_site, lat_site) %>%
  arrange(lon_site, lat_site) %>%
  mutate(site_name = paste0("site_", row_number()))

dat_sites <- dat_sites %>%
  left_join(site_lut, by = c("site_id", "lon_site", "lat_site")) %>%
  mutate(site_name = coalesce(site_name.y, site_name.x)) %>%
  select(-any_of(c("site_name.x", "site_name.y")))

# 2. Spatial extent helpers ------------------------------------------------
# Rough geographic scaling for quick “extent” summaries.
# This is a screening metric, not a precise geodesic distance.
km_per_deg_lat <- 111.32

# 3. Host-species spatial coverage ----------------------------------------
spatial_coverage_species <- dat_sites %>%
  group_by(host_species) %>%
  summarise(
    n_rows = n(),
    n_sites = n_distinct(site_name),
    total_tested = sum(number_tested, na.rm = TRUE),
    total_positive = sum(number_positive, na.rm = TRUE),
    prop_zero_rows = mean(number_positive == 0, na.rm = TRUE),
    boundary_d_q10 = as.numeric(stats::quantile(Boundary_d, 0.1, na.rm = TRUE)),
    boundary_d_q50 = as.numeric(stats::quantile(Boundary_d, 0.5, na.rm = TRUE)),
    boundary_d_q90 = as.numeric(stats::quantile(Boundary_d, 0.9, na.rm = TRUE)),
    centroid_d_q10 = as.numeric(stats::quantile(Centroid_d, 0.1, na.rm = TRUE)),
    centroid_d_q50 = as.numeric(stats::quantile(Centroid_d, 0.5, na.rm = TRUE)),
    centroid_d_q90 = as.numeric(stats::quantile(Centroid_d, 0.9, na.rm = TRUE)),
    lon_min = min(longitude, na.rm = TRUE),
    lon_max = max(longitude, na.rm = TRUE),
    lat_min = min(latitude, na.rm = TRUE),
    lat_max = max(latitude, na.rm = TRUE),
    mean_lat = mean(latitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prevalence = if_else(total_tested > 0, total_positive / total_tested, NA_real_),
    lon_range_deg = lon_max - lon_min,
    lat_range_deg = lat_max - lat_min,
    lat_km = lat_range_deg * km_per_deg_lat,
    lon_km = lon_range_deg * abs(km_per_deg_lat * cos(mean_lat * pi / 180)),
    bbox_diag_km = sqrt(lat_km^2 + lon_km^2),
    bbox_area_km2_approx = lat_km * lon_km
  ) %>%
  select(-lon_min, -lon_max, -lat_min, -lat_max, -mean_lat, -lat_km, -lon_km) %>%
  arrange(desc(n_sites))

# 4. Host–pathogen spatial coverage ---------------------------------------
spatial_coverage_hp <- dat_sites %>%
  filter(!is.na(pathogen_species_cleaned)) %>%
  mutate(host_pathogen = paste(host_species, pathogen_species_cleaned, sep = " | ")) %>%
  group_by(host_pathogen, host_species, pathogen_species_cleaned) %>%
  summarise(
    assay_groups_n = n_distinct(assay_group),
    n_rows = n(),
    n_sites = n_distinct(site_name),
    total_tested = sum(number_tested, na.rm = TRUE),
    total_positive = sum(number_positive, na.rm = TRUE),
    prop_zero_rows = mean(number_positive == 0, na.rm = TRUE),
    boundary_d_q10 = as.numeric(stats::quantile(Boundary_d, 0.1, na.rm = TRUE)),
    boundary_d_q50 = as.numeric(stats::quantile(Boundary_d, 0.5, na.rm = TRUE)),
    boundary_d_q90 = as.numeric(stats::quantile(Boundary_d, 0.9, na.rm = TRUE)),
    centroid_d_q10 = as.numeric(stats::quantile(Centroid_d, 0.1, na.rm = TRUE)),
    centroid_d_q50 = as.numeric(stats::quantile(Centroid_d, 0.5, na.rm = TRUE)),
    centroid_d_q90 = as.numeric(stats::quantile(Centroid_d, 0.9, na.rm = TRUE)),
    lon_min = min(longitude, na.rm = TRUE),
    lon_max = max(longitude, na.rm = TRUE),
    lat_min = min(latitude, na.rm = TRUE),
    lat_max = max(latitude, na.rm = TRUE),
    mean_lat = mean(latitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prevalence = if_else(total_tested > 0, total_positive / total_tested, NA_real_),
    lon_range_deg = lon_max - lon_min,
    lat_range_deg = lat_max - lat_min,
    lat_km = lat_range_deg * km_per_deg_lat,
    lon_km = lon_range_deg * abs(km_per_deg_lat * cos(mean_lat * pi / 180)),
    bbox_diag_km = sqrt(lat_km^2 + lon_km^2),
    bbox_area_km2_approx = lat_km * lon_km
  ) %>%
  select(-lon_min, -lon_max, -lat_min, -lat_max, -mean_lat, -lat_km, -lon_km) %>%
  arrange(desc(n_sites))

# 5. Export ---------------------------------------------------------------
out_dir <- analysis_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  spatial_coverage_species,
  here(out_dir, "spatial_coverage_host_species.csv"),
  row.names = FALSE
)

write.csv(
  spatial_coverage_hp,
  here(out_dir, "spatial_coverage_host_pathogen.csv"),
  row.names = FALSE
)

cat("\n--- SPATIAL COVERAGE: HOST SPECIES (top 10 by n_sites) ---\n")
print(spatial_coverage_species %>%
  select(host_species, n_rows, n_sites, bbox_diag_km, prop_zero_rows) %>%
  slice_head(n = 10))

cat("\n--- SPATIAL COVERAGE: HOST × PATHOGEN (top 10 by n_sites) ---\n")
print(spatial_coverage_hp %>%
  select(host_pathogen, n_rows, n_sites, bbox_diag_km, assay_groups_n) %>%
  slice_head(n = 10))
