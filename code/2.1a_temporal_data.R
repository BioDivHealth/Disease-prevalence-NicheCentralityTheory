# 2.1a_temporal_data
# Goal: add temporal covariates and optional MODIS MidGreenup metadata.
# Scope: row-level data engineering only. Temporal candidate/audit tables live in
# code/2.1a_temporal_candidates.R.

# Load libraries
library(pacman)
# Note: The package name is case-sensitive: MODISTools
p_load(tidyverse, lubridate, here, MODISTools)

# Parameters ---------------------------------------------------------------
# Spatial mode: exact vs rounded coordinates
site_mode <- "exact" # "exact" or "rounded"
exact_round_digits <- 3
coord_round_digits <- 2

mode_tag <- if (site_mode == "exact") "exact" else paste0("rounded", coord_round_digits)

# MODIS phenology extraction ------------------------------------------------
# Keep opt-in: avoids re-downloading and allows temporal metadata runs without MODIS.
use_modis <- TRUE
run_modis_download <- FALSE

# 1. Load & pre-filter data ------------------------------------------------
path_dat_clean <- here("Data", "dat_clean_agg2.rds")
dat_clean_download_url <- "https://github.com/BioDivHealth/Disease-prevalence-NicheCentralityTheory/blob/artur/Data/dat_clean_agg2.rds"

if (!file.exists(path_dat_clean)) {
  stop(
    "Missing required input file: ", path_dat_clean,
    "\nDownload or restore it from: ", dat_clean_download_url,
    "\nThen rerun code/2.1a_temporal_data.R.",
    call. = FALSE
  )
}

dat <- readRDS(path_dat_clean)

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

# Helpful QA
cat("\nSummary of sampling window lengths (days):\n")
print(summary(dat_year$sampling_duration_days))

invisible(list(
  dat_year = dat_year,
  dat_final = dat_final
))
