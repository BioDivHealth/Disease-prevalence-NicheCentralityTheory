# 2.1a_temporal_data
# Goal: add temporal covariates + identify temporal candidate species / host–pathogen pairs.
# Scope: data engineering + candidate filtering only (plots live in 2.1b).

# Load libraries
library(pacman)
# Note: The package name is case-sensitive: MODISTools
p_load(tidyverse, lubridate, here, MODISTools)

# Parameters ---------------------------------------------------------------
max_period_default <- 15
max_period_sensitivity <- c(7, 30, 60, 120)
max_period_grid <- sort(unique(c(max_period_default, max_period_sensitivity)))

# Spatial mode: exact vs rounded coordinates
site_mode <- "exact" # "exact" or "rounded"
exact_round_digits <- 3
coord_round_digits <- 2

mode_tag <- if (site_mode == "exact") "exact" else paste0("rounded", coord_round_digits)
analysis_dir <- here("Results", "analysis_metadata", mode_tag)

# MODIS phenology extraction ------------------------------------------------
# Keep opt-in: avoids re-downloading and allows temporal audit without MODIS.
use_modis <- TRUE
run_modis_download <- FALSE

# Helpers ------------------------------------------------------------------
gap_vec <- function(dates) {
  dates <- sort(unique(dates[!is.na(dates)]))
  if (length(dates) < 2) return(numeric(0))
  as.numeric(diff(dates), units = "days")
}

add_gap_summaries <- function(df, gaps_col = "gaps"){
  gaps <- df[[gaps_col]]
  df %>%
    mutate(
      avg_gap_days = map_dbl(gaps, ~ if (length(.x) > 0) mean(.x) else NA_real_),
      median_gap_days = map_dbl(gaps, ~ if (length(.x) > 0) median(.x) else NA_real_),
      min_gap_days = map_dbl(gaps, ~ if (length(.x) > 0) min(.x) else NA_real_),
      max_gap_days = map_dbl(gaps, ~ if (length(.x) > 0) max(.x) else NA_real_),
      gap_q75_days = map_dbl(gaps, ~ if (length(.x) > 0) as.numeric(stats::quantile(.x, 0.75)) else NA_real_)
    )
}

make_temporal_frame <- function(data, max_period) {
  data %>%
    mutate(
      max_period = max_period,
      short_window = days_diff <= max_period,
      # Use midpoint_date (Date) for assigning an event date for short windows
      event_date = if_else(short_window, midpoint_date, as.Date(NA)),
      coverage_days = if_else(short_window, as.numeric(days_diff) + 1, 0)
    )
}

# 1. Load & pre-filter data ------------------------------------------------
dat <- readRDS(here("Data", "dat_clean_agg2.rds"))

# Keep records with <= 1 year sampling window; keep hosts with enough records
# (Used for temporal audits + MODIS join)
dat_year <- dat %>%
  filter(days_diff < 370) %>%
  group_by(host_species) %>%
  filter(n() >= 15) %>%
  ungroup()

cat("Number of unique host species in dat_year: ", length(unique(dat_year$host_species)))
cat("\nNumber of rows in dat_year: ", nrow(dat_year), "\n")

# 2. Temporal variables ----------------------------------------------------
# The date_interval column is a lubridate interval object
# Use midpoint for assigning seasonal/yearly context
dat_year <- dat_year %>%
  mutate(
    start_date = as.Date(int_start(date_interval)),
    end_date = as.Date(int_end(date_interval)),
    midpoint_date = as.Date(start_date + (end_date - start_date) / 2),
    month = month(midpoint_date),
    quarter = quarter(midpoint_date),
    year = year(midpoint_date),
    doy = yday(midpoint_date),
    doy_sin = sin(2 * pi * doy / 365.25),
    doy_cos = cos(2 * pi * doy / 365.25),
    sampling_duration_days = as.numeric(days_diff),
    lon_site = round(longitude, if_else(site_mode == "exact", exact_round_digits, coord_round_digits)),
    lat_site = round(latitude, if_else(site_mode == "exact", exact_round_digits, coord_round_digits))
  ) %>%
  # Avoid carrying Interval columns into dplyr summarise workflows
  select(-date_interval)

# 3. Site definitions (mode-specific coordinates) ---------------------------
# Used for: (a) reducing MODIS API calls, (b) measuring repeat sampling at same locality
sites <- dat_year %>%
  distinct(lat_site, lon_site) %>%
  arrange(lon_site, lat_site) %>%
  mutate(
    site_id = paste(lon_site, lat_site, sep = "_"),
    site_name = paste0("site_", row_number())
  )

# Join site_name back to main data for later merging
dat_year <- dat_year %>%
  left_join(sites, by = c("lat_site", "lon_site"))

# 4. MODIS phenology metadata (MCD12Q2: MidGreenup) ------------------------
# Notes:
# - Product: MCD12Q2 (Land Cover Dynamics)
# - Band: MidGreenup.Num_Modes_01
# - MODIS availability: 2000-02-18 onwards

years_needed <- sort(unique(dat_year$year))
start_year <- max(2000, min(years_needed, na.rm = TRUE))
end_year <- max(years_needed, na.rm = TRUE)

message(paste("MODIS availability: 2000-present. Data includes years back to", min(years_needed)))
message(paste(
  "Ready to extract MODIS data for",
  nrow(sites),
  "unique", mode_tag, "sites for",
  start_year, "to", end_year
))

# MODIS cache paths
modis_path_round <- here("Data", paste0("modis_mcd12q2_raw_round", coord_round_digits, ".rds"))
modis_path_exact <- here("Data", "modis_mcd12q2_raw.rds")

modis_path_mode <- if (site_mode == "rounded") modis_path_round else modis_path_exact
modis_source <- if (file.exists(modis_path_mode)) site_mode else "missing"

if (isTRUE(use_modis)) {
  message(
    "MODIS input source: ", modis_source,
    if (modis_source != "missing") paste0(" (", basename(modis_path_mode), ")") else ""
  )
}

if (isTRUE(use_modis) && isTRUE(run_modis_download)) {
  modis_res <- mt_batch_subset(
    df = sites %>% rename(lat = lat_site, lon = lon_site),
    product = "MCD12Q2",
    band = "MidGreenup.Num_Modes_01",
    start = paste0(start_year, "-01-01"),
    end = paste0(end_year, "-12-31")
  )
  saveRDS(modis_res, modis_path_mode)
  modis_source <- site_mode
}

# Default: build dat_final even if MODIS is missing
# (keeps downstream code stable and avoids failing temporal audit runs).
dat_final <- dat_year %>%
  mutate(
    greenup_date = as.Date(NA),
    mean_doy = as.numeric(NA),
    pheno_source = "Missing",
    effective_greenup = as.Date(NA),
    days_since_midgreenup = as.numeric(NA)
  )

if (isTRUE(use_modis) && modis_source != "missing") {
  if (!exists("modis_res")) {
    modis_res <- readRDS(modis_path_mode)
  }

  # Standardize coordinate column names for exact-mode files.
  if (modis_source == "exact") {
    if ("latitude" %in% names(modis_res) && !("lat" %in% names(modis_res))) {
      modis_res <- dplyr::rename(modis_res, lat = latitude)
    }
    if ("longitude" %in% names(modis_res) && !("lon" %in% names(modis_res))) {
      modis_res <- dplyr::rename(modis_res, lon = longitude)
    }
    if (!("lat" %in% names(modis_res)) || !("lon" %in% names(modis_res))) {
      stop(
        "MODIS file exists but lacks lat/lon columns: ", modis_path_mode,
        "\nExpected columns 'lat'/'lon' (or 'latitude'/'longitude')."
      )
    }
  }

  if (!("calendar_date" %in% names(modis_res)) || !("value" %in% names(modis_res))) {
    stop("Unexpected MODIS file format: missing calendar_date/value in ", modis_path_mode)
  }

  if (modis_source == "rounded") {
    modis_clean <- modis_res %>%
      filter(value != 32767) %>%
      mutate(
        greenup_date = as.Date(value, origin = "1970-01-01"),
        modis_year = year(as.Date(calendar_date))
      ) %>%
      transmute(site_name = site, modis_year, greenup_date) %>%
      group_by(site_name, modis_year) %>%
      summarise(
        greenup_date = as.Date(round(mean(as.numeric(greenup_date), na.rm = TRUE)), origin = "1970-01-01"),
        .groups = "drop"
      )
  } else {
    modis_clean <- modis_res %>%
      filter(value != 32767) %>%
      mutate(
        greenup_date = as.Date(value, origin = "1970-01-01"),
        modis_year = year(as.Date(calendar_date)),
        lat_site = round(lat, exact_round_digits),
        lon_site = round(lon, exact_round_digits)
      ) %>%
      left_join(sites, by = c("lat_site", "lon_site")) %>%
      filter(!is.na(site_name)) %>%
      group_by(site_name, modis_year) %>%
      summarise(
        greenup_date = as.Date(round(mean(as.numeric(greenup_date), na.rm = TRUE)), origin = "1970-01-01"),
        .groups = "drop"
      )
  }

  site_climatology <- modis_clean %>%
    mutate(doy = yday(greenup_date)) %>%
    group_by(site_name) %>%
    summarise(mean_doy = mean(doy, na.rm = TRUE), .groups = "drop")

  dat_final <- dat_year %>%
    left_join(modis_clean, by = c("site_name", "year" = "modis_year")) %>%
    left_join(site_climatology, by = "site_name") %>%
    mutate(
      pheno_source = case_when(
        !is.na(greenup_date) ~ "Year-specific",
        !is.na(mean_doy) ~ "Climatology (Site Avg)",
        TRUE ~ "Missing"
      ),
      effective_greenup = case_when(
        pheno_source == "Year-specific" ~ greenup_date,
        !is.na(mean_doy) ~ as.Date(paste0(year, "-01-01")) + days(round(mean_doy) - 1),
        TRUE ~ as.Date(NA)
      ),
      days_since_midgreenup = as.numeric(difftime(midpoint_date, effective_greenup, units = "days"))
    )

  cat("\nPhenology Data Match Summary:\n")
  print(table(dat_final$pheno_source))
} else if (isTRUE(use_modis)) {
  warning(
    "MODIS phenology file not found for mode '", site_mode, "'. Tried:\n- ",
    modis_path_mode,
    "\nContinuing without MODIS metadata."
  )
}

# Save temporal-only metadata product (always)
data_temporal_rds <- here("Data", paste0("dat_with_temporal_metadata_", mode_tag, ".rds"))
data_temporal_csv <- here("Data", paste0("dat_with_temporal_metadata_", mode_tag, ".csv"))

data_modis_rds <- here("Data", paste0("dat_with_modis_metadata_", mode_tag, ".rds"))
data_modis_csv <- here("Data", paste0("dat_with_modis_metadata_", mode_tag, ".csv"))

saveRDS(dat_year, data_temporal_rds)
write.csv(dat_year, data_temporal_csv, row.names = FALSE)

# Save MODIS-enriched metadata only when MODIS is available
# (prevents overwriting a previously-built MODIS file with NA placeholders).
if (isTRUE(use_modis) && modis_source != "missing") {
  saveRDS(dat_final, data_modis_rds)
  write.csv(dat_final, data_modis_csv, row.names = FALSE)
} else {
  message("Skipping write of ", basename(data_modis_rds), " (MODIS not available).")
}

# 5. General temporal audit (candidates) -----------------------------------
# Treat records as:
# - short_window: days_diff <= max_period
# - long_window:  days_diff >  max_period
# For short windows, event_date = midpoint_date.

# Default period (used by 2.1b plots)
dat_year_temporal <- make_temporal_frame(dat_year, max_period_default)

# 5.1 Candidate hosts (species level)
# NOTE: This aggregates across pathogens; it is mainly a data quality summary.
species_site_rows <- dat_year_temporal %>%
  filter(short_window) %>%
  count(host_species, site_name, name = "n_short_rows_site")

species_top_site <- species_site_rows %>%
  group_by(host_species) %>%
  summarise(
    n_sites_short = n_distinct(site_name),
    top_site_short_rows = max(n_short_rows_site, na.rm = TRUE),
    top_site_short_rows_share = top_site_short_rows / sum(n_short_rows_site),
    .groups = "drop"
  )

temporal_candidates_species_all <- dat_year_temporal %>%
  group_by(host_species) %>%
  summarise(
    max_period = first(max_period),
    n_total_records = n(),
    n_sites_total = n_distinct(site_name),
    total_tested = sum(number_tested, na.rm = TRUE),
    total_positive = sum(number_positive, na.rm = TRUE),
    prevalence = if_else(total_tested > 0, total_positive / total_tested, NA_real_),
    prop_zero_rows = mean(number_positive == 0, na.rm = TRUE),
    n_short_windows = sum(short_window, na.rm = TRUE),
    n_long_windows = sum(!short_window, na.rm = TRUE),
    n_event_days = n_distinct(event_date, na.rm = TRUE),
    event_span_days = {
      d <- event_date[!is.na(event_date)]
      if (length(d) > 1) as.numeric(difftime(max(d), min(d), units = "days")) else NA_real_
    },
    n_months = n_distinct(month(event_date), na.rm = TRUE),
    n_quarters = n_distinct(quarter(event_date), na.rm = TRUE),
    n_years = n_distinct(year(event_date), na.rm = TRUE),
    gaps = list(gap_vec(event_date)),
    .groups = "drop"
  ) %>%
  add_gap_summaries("gaps") %>%
  left_join(species_top_site, by = "host_species") %>%
  select(-gaps)

# Filtered candidate subset for convenience
temporal_candidates_species <- temporal_candidates_species_all %>%
  filter(n_event_days >= 30, !is.na(avg_gap_days)) %>%
  arrange(avg_gap_days)

cat("\n--- GENERAL TEMPORAL AUDIT (Host species) ---\n")
cat("Top 10 host species with most consistent reporting:\n")
print(
  temporal_candidates_species %>%
    select(host_species, n_event_days, avg_gap_days, max_gap_days, n_sites_total, n_sites_short, top_site_short_rows_share) %>%
    slice_head(n = 10)
)

# 5.2 Candidate host x pathogen
# This is closer to the modelling unit and should be preferred for pilot selection.
hp_site_rows <- dat_year_temporal %>%
  filter(!is.na(pathogen_species_cleaned), short_window) %>%
  mutate(host_pathogen = paste(host_species, pathogen_species_cleaned, sep = " | ")) %>%
  count(host_pathogen, site_name, name = "n_short_rows_site")

hp_top_site <- hp_site_rows %>%
  group_by(host_pathogen) %>%
  summarise(
    n_sites_short = n_distinct(site_name),
    top_site_short_rows = max(n_short_rows_site, na.rm = TRUE),
    top_site_short_rows_share = top_site_short_rows / sum(n_short_rows_site),
    .groups = "drop"
  )

temporal_candidates_hp_all <- dat_year_temporal %>%
  filter(!is.na(pathogen_species_cleaned)) %>%
  mutate(host_pathogen = paste(host_species, pathogen_species_cleaned, sep = " | ")) %>%
  group_by(host_pathogen, host_species, pathogen_species_cleaned) %>%
  summarise(
    max_period = first(max_period),
    assay_groups_n = n_distinct(assay_group),
    n_total_records = n(),
    n_sites_total = n_distinct(site_name),
    total_tested = sum(number_tested, na.rm = TRUE),
    total_positive = sum(number_positive, na.rm = TRUE),
    prevalence = if_else(total_tested > 0, total_positive / total_tested, NA_real_),
    prop_zero_rows = mean(number_positive == 0, na.rm = TRUE),
    n_short_windows = sum(short_window, na.rm = TRUE),
    n_long_windows = sum(!short_window, na.rm = TRUE),
    n_event_days = n_distinct(event_date, na.rm = TRUE),
    event_span_days = {
      d <- event_date[!is.na(event_date)]
      if (length(d) > 1) as.numeric(difftime(max(d), min(d), units = "days")) else NA_real_
    },
    n_months = n_distinct(month(event_date), na.rm = TRUE),
    n_quarters = n_distinct(quarter(event_date), na.rm = TRUE),
    n_years = n_distinct(year(event_date), na.rm = TRUE),
    gaps = list(gap_vec(event_date)),
    .groups = "drop"
  ) %>%
  add_gap_summaries("gaps") %>%
  left_join(hp_top_site, by = "host_pathogen") %>%
  select(-gaps)

temporal_candidates_hp <- temporal_candidates_hp_all %>%
  filter(n_event_days >= 20, !is.na(avg_gap_days)) %>%
  arrange(avg_gap_days)

cat("\n--- TEMPORAL AUDIT (Host × Pathogen) ---\n")
cat("Top 10 combos with most consistent reporting:\n")
print(
  temporal_candidates_hp %>%
    select(host_pathogen, n_event_days, avg_gap_days, max_gap_days, n_sites_total, n_sites_short, top_site_short_rows_share, assay_groups_n) %>%
    slice_head(n = 10)
)

# 5.3 Candidate host x pathogen x site (within-locality time series)
temporal_candidates_hp_site_all <- dat_year_temporal %>%
  filter(!is.na(pathogen_species_cleaned)) %>%
  group_by(host_species, pathogen_species_cleaned, site_name) %>%
  summarise(
    max_period = first(max_period),
    assay_groups_n = n_distinct(assay_group),
    n_total_records = n(),
    total_tested = sum(number_tested, na.rm = TRUE),
    total_positive = sum(number_positive, na.rm = TRUE),
    prevalence = if_else(total_tested > 0, total_positive / total_tested, NA_real_),
    prop_zero_rows = mean(number_positive == 0, na.rm = TRUE),
    n_short_windows = sum(short_window, na.rm = TRUE),
    n_long_windows = sum(!short_window, na.rm = TRUE),
    n_event_days = n_distinct(event_date, na.rm = TRUE),
    event_span_days = {
      d <- event_date[!is.na(event_date)]
      if (length(d) > 1) as.numeric(difftime(max(d), min(d), units = "days")) else NA_real_
    },
    n_months = n_distinct(month(event_date), na.rm = TRUE),
    n_quarters = n_distinct(quarter(event_date), na.rm = TRUE),
    n_years = n_distinct(year(event_date), na.rm = TRUE),
    gaps = list(gap_vec(event_date)),
    site_lon = first(lon_site),
    site_lat = first(lat_site),
    .groups = "drop"
  ) %>%
  add_gap_summaries("gaps") %>%
  select(-gaps) %>%
  arrange(desc(n_event_days), avg_gap_days)

# A very strict “pilot shortlist” to explore within-locality time series with enough events and multi-season coverage.
pilot_hp_site <- temporal_candidates_hp_site_all %>%
  filter(n_event_days >= 30, n_quarters >= 3, !is.na(avg_gap_days))

# 6. Exports ---------------------------------------------------------------
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)

# Default-period objects used by 2.1b plotting
saveRDS(dat_year_temporal, here(analysis_dir, "dat_year_temporal.rds"))
saveRDS(temporal_candidates_species, here(analysis_dir, "temporal_candidates_species.rds"))
write.csv(
  temporal_candidates_species,
  here(analysis_dir, "temporal_candidates_species.csv"),
  row.names = FALSE
)

# Host-pathogen default
saveRDS(
  dat_year_temporal %>%
    filter(!is.na(pathogen_species_cleaned)) %>%
    mutate(host_pathogen = paste(host_species, pathogen_species_cleaned, sep = " | ")),
  here(analysis_dir, "dat_year_temporal_host_pathogen.rds")
)
saveRDS(temporal_candidates_hp, here(analysis_dir, "temporal_candidates_host_pathogen.rds"))
write.csv(
  temporal_candidates_hp,
  here(analysis_dir, "temporal_candidates_host_pathogen.csv"),
  row.names = FALSE
)

# Also export full (unfiltered) summaries for auditing
saveRDS(temporal_candidates_species_all, here(analysis_dir, "temporal_candidates_species_all.rds"))
write.csv(
  temporal_candidates_species_all,
  here(analysis_dir, "temporal_candidates_species_all.csv"),
  row.names = FALSE
)

saveRDS(temporal_candidates_hp_all, here(analysis_dir, "temporal_candidates_host_pathogen_all.rds"))
write.csv(
  temporal_candidates_hp_all,
  here(analysis_dir, "temporal_candidates_host_pathogen_all.csv"),
  row.names = FALSE
)

# Host-pathogen-site time series
saveRDS(temporal_candidates_hp_site_all, here(analysis_dir, "temporal_candidates_host_pathogen_site_all.rds"))
write.csv(
  temporal_candidates_hp_site_all,
  here(analysis_dir, "temporal_candidates_host_pathogen_site_all.csv"),
  row.names = FALSE
)

saveRDS(pilot_hp_site, here(analysis_dir, "pilot_temporal_hp_site_candidates.rds"))
write.csv(
  pilot_hp_site,
  here(analysis_dir, "pilot_temporal_hp_site_candidates.csv"),
  row.names = FALSE
)

# Gap distributions (for thinking about temporal independence thresholds)
gaps_hp_site <- dat_year_temporal %>%
  filter(!is.na(pathogen_species_cleaned), short_window, !is.na(event_date)) %>%
  group_by(host_species, pathogen_species_cleaned, site_name) %>%
  summarise(event_dates = list(sort(unique(event_date))), .groups = "drop") %>%
  mutate(gap_days = map(event_dates, gap_vec)) %>%
  select(-event_dates) %>%
  unnest(gap_days)

write.csv(
  gaps_hp_site,
  here(analysis_dir, "temporal_gaps_host_pathogen_site_max15.csv"),
  row.names = FALSE
)

# Assay-stratified host-pathogen candidates (confounding check)
assay_levels <- sort(unique(dat_year_temporal$assay_group))
for (a in assay_levels) {
  df_a <- dat_year_temporal %>% filter(assay_group == a)
  cand_a <- df_a %>%
    filter(!is.na(pathogen_species_cleaned)) %>%
    mutate(host_pathogen = paste(host_species, pathogen_species_cleaned, sep = " | ")) %>%
    group_by(host_pathogen, host_species, pathogen_species_cleaned) %>%
    summarise(
      max_period = first(max_period),
      assay_group = first(assay_group),
      n_total_records = n(),
      n_sites_total = n_distinct(site_name),
      total_tested = sum(number_tested, na.rm = TRUE),
      total_positive = sum(number_positive, na.rm = TRUE),
      prevalence = if_else(total_tested > 0, total_positive / total_tested, NA_real_),
      n_short_windows = sum(short_window, na.rm = TRUE),
      n_event_days = n_distinct(event_date, na.rm = TRUE),
      gaps = list(gap_vec(event_date)),
      .groups = "drop"
    ) %>%
    add_gap_summaries("gaps") %>%
    select(-gaps) %>%
    filter(n_event_days >= 10, !is.na(avg_gap_days)) %>%
    arrange(avg_gap_days)

  a_safe <- stringr::str_replace_all(a, "[^A-Za-z0-9]+", "_")
  out_a <- here(analysis_dir, paste0("temporal_candidates_host_pathogen_assay_", a_safe, ".csv"))
  write.csv(cand_a, out_a, row.names = FALSE)
}

# Sensitivity tables for other max_period values (filtered candidates only)
sens_dir <- here(analysis_dir, "temporal_sensitivity")
dir.create(sens_dir, showWarnings = FALSE, recursive = TRUE)

for (mp in max_period_grid) {
  df_mp <- make_temporal_frame(dat_year, mp)

  cand_sp_mp <- df_mp %>%
    group_by(host_species) %>%
    summarise(
      max_period = first(max_period),
      n_total_records = n(),
      n_sites_total = n_distinct(site_name),
      n_short_windows = sum(short_window, na.rm = TRUE),
      n_event_days = n_distinct(event_date, na.rm = TRUE),
      gaps = list(gap_vec(event_date)),
      .groups = "drop"
    ) %>%
    add_gap_summaries("gaps") %>%
    select(-gaps) %>%
    filter(n_event_days >= 30, !is.na(avg_gap_days)) %>%
    arrange(avg_gap_days)

  cand_hp_mp <- df_mp %>%
    filter(!is.na(pathogen_species_cleaned)) %>%
    mutate(host_pathogen = paste(host_species, pathogen_species_cleaned, sep = " | ")) %>%
    group_by(host_pathogen) %>%
    summarise(
      max_period = first(max_period),
      assay_groups_n = n_distinct(assay_group),
      n_total_records = n(),
      n_sites_total = n_distinct(site_name),
      n_short_windows = sum(short_window, na.rm = TRUE),
      n_event_days = n_distinct(event_date, na.rm = TRUE),
      gaps = list(gap_vec(event_date)),
      .groups = "drop"
    ) %>%
    add_gap_summaries("gaps") %>%
    select(-gaps) %>%
    filter(n_event_days >= 20, !is.na(avg_gap_days)) %>%
    arrange(avg_gap_days)

  mp_tag <- stringr::str_pad(mp, width = 3, side = "left", pad = "0")
  write.csv(cand_sp_mp, here(sens_dir, paste0("temporal_candidates_species_max", mp_tag, ".csv")), row.names = FALSE)
  write.csv(cand_hp_mp, here(sens_dir, paste0("temporal_candidates_host_pathogen_max", mp_tag, ".csv")), row.names = FALSE)
}

# Helpful QA
cat("\nSummary of sampling window lengths (days):\n")
print(summary(dat_year$sampling_duration_days))

invisible(list(
  dat_year = dat_year,
  dat_final = dat_final,
  temporal_candidates_species = temporal_candidates_species,
  temporal_candidates_hp = temporal_candidates_hp
))
