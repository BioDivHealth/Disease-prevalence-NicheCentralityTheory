# 2.1d_spatial_coverage_plots
# Goal: plot spatial coverage for top host species and host–pathogen combos.

library(pacman)
p_load(tidyverse, lubridate, here, maps, patchwork)

# Parameters ---------------------------------------------------------------
site_mode <- "exact" # "exact" or "rounded"
exact_round_digits <- 3
coord_round_digits <- 2
max_top <- 10
inset_pos <- list(left = 0.55, bottom = 0.75, right = 1, top = 1)

mode_tag <- if (site_mode == "exact") "exact" else paste0("rounded", coord_round_digits)
analysis_dir <- here("Results", "analysis_metadata", mode_tag)
plots_dir <- here("Results", "plots", mode_tag)

# Inputs ------------------------------------------------------------------
path_cov_species <- here(analysis_dir, "spatial_coverage_host_species.csv")
path_cov_hp <- here(analysis_dir, "spatial_coverage_host_pathogen.csv")
path_dat_temporal <- here("Data", paste0("dat_with_temporal_metadata_", mode_tag, ".rds"))
path_dat_fallback <- here("Data", "dat_clean_agg2.rds")

stopifnot(file.exists(path_cov_species))
stopifnot(file.exists(path_cov_hp))
stopifnot(file.exists(path_dat_fallback))

coverage_species <- read.csv(path_cov_species, stringsAsFactors = FALSE)
coverage_hp <- read.csv(path_cov_hp, stringsAsFactors = FALSE)

dat <- if (file.exists(path_dat_temporal)) {
  readRDS(path_dat_temporal)
} else {
  readRDS(path_dat_fallback)
}

dat <- dat %>%
  select(-any_of("date_interval"))

# Add/fill a year column if temporal metadata exists
date_fields <- intersect(c("event_date", "midpoint_date", "start_date", "end_date"), names(dat))
if (length(date_fields) > 0) {
  for (col in date_fields) {
    dat[[col]] <- as.Date(dat[[col]])
  }
}

year_vals <- if ("year" %in% names(dat)) dat$year else rep(NA_integer_, nrow(dat))
if (length(date_fields) > 0) {
  for (col in date_fields) {
    year_vals <- coalesce(year_vals, lubridate::year(dat[[col]]))
  }
}

dat$year <- year_vals

# Define sites using mode-specific coordinates to match 2.1c
sites <- dat %>%
  mutate(
    lon_site = round(longitude, if_else(site_mode == "exact", exact_round_digits, coord_round_digits)),
    lat_site = round(latitude, if_else(site_mode == "exact", exact_round_digits, coord_round_digits)),
    site_id = paste(lon_site, lat_site, sep = "_")
  )

# Ranking helpers ----------------------------------------------------------
rank_by_coverage <- function(df) {
  df %>%
    mutate(
      n_sites_scaled = as.numeric(scale(n_sites)),
      bbox_scaled = as.numeric(scale(bbox_diag_km)),
      rank_score = n_sites_scaled + bbox_scaled
    ) %>%
    arrange(desc(rank_score))
}

coverage_species_ranked <- rank_by_coverage(coverage_species)
coverage_hp_ranked <- rank_by_coverage(coverage_hp)

selected_species <- coverage_species_ranked %>%
  slice_head(n = min(max_top, nrow(coverage_species_ranked)))

selected_hp <- coverage_hp_ranked %>%
  slice_head(n = min(max_top, nrow(coverage_hp_ranked)))

# Plot helpers -------------------------------------------------------------
world <- map_data("world")

create_inset_map <- function(plot_data, title, subtitle, inset_rect, inset_pos) {
  year_limits <- range(plot_data$year, na.rm = TRUE)
  if (!all(is.finite(year_limits))) {
    year_limits <- NULL
  }

  year_breaks <- NULL
  if (any(!is.na(plot_data$year))) {
    year_breaks <- scales::pretty_breaks(n = 4)(plot_data$year)
  }

  span_ratio <- diff(inset_rect$xlims) / diff(inset_rect$ylims)
  barwidth_cm <- max(2.5, min(8, 8 * span_ratio))

  p_main <- ggplot() +
    geom_polygon(
      data = world,
      aes(x = long, y = lat, group = group),
      fill = "grey92",
      color = "white",
      linewidth = 0.1
    ) +
    geom_point(
      data = plot_data,
      aes(x = longitude, y = latitude, color = year),
      size = 1.6,
      alpha = 0.45,
      position = position_jitter(width = 0.1, height = 0.1)
    ) +
    coord_quickmap(xlim = inset_rect$xlims, ylim = inset_rect$ylims) +
    scale_color_viridis_c(
      option = "D",
      limits = year_limits,
      breaks = year_breaks,
      na.value = "grey70",
      guide = guide_colorbar(
        title.position = "top",
        barwidth = unit(barwidth_cm, "cm"),
        barheight = unit(0.4, "cm")
      )
    ) +
    labs(title = title, subtitle = subtitle, color = "Year") +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "aliceblue", color = "grey80"),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  p_inset <- ggplot() +
    geom_polygon(
      data = world,
      aes(x = long, y = lat, group = group),
      fill = "grey75",
      color = "white",
      linewidth = 0.05
    ) +
    geom_rect(
      aes(
        xmin = inset_rect$xlims[1],
        xmax = inset_rect$xlims[2],
        ymin = inset_rect$ylims[1],
        ymax = inset_rect$ylims[2]
      ),
      color = "red",
      fill = "red",
      alpha = 0.2,
      linewidth = 0.2
    ) +
    coord_quickmap() +
    theme_void() +
    theme(panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3))

  p_main + inset_element(
    p_inset,
    left = inset_pos$left,
    bottom = inset_pos$bottom,
    right = inset_pos$right,
    top = inset_pos$top,
    align_to = "panel"
  )
}

compute_extent <- function(plot_data) {
  lon_range <- range(plot_data$longitude, na.rm = TRUE)
  lat_range <- range(plot_data$latitude, na.rm = TRUE)

  lon_span <- max(diff(lon_range), 0.1)
  lat_span <- max(diff(lat_range), 0.1)

  min_ratio <- 0.5
  min_span_deg <- 2

  target_lon_span <- max(lon_span, lat_span * min_ratio, min_span_deg)
  target_lat_span <- max(lat_span, lon_span * min_ratio, min_span_deg)

  lon_center <- mean(lon_range)
  lat_center <- mean(lat_range)

  lon_span <- target_lon_span
  lat_span <- target_lat_span

  lon_buffer <- max(lon_span * 0.3, 1.5)
  lat_buffer <- max(lat_span * 0.3, 1.5)

  xlims <- c(lon_center - lon_span / 2 - lon_buffer, lon_center + lon_span / 2 + lon_buffer)
  ylims <- c(lat_center - lat_span / 2 - lat_buffer, lat_center + lat_span / 2 + lat_buffer)

  xlims <- c(max(-180, xlims[1]), min(180, xlims[2]))
  ylims <- c(max(-90, ylims[1]), min(90, ylims[2]))

  list(xlims = xlims, ylims = ylims)
}

format_subtitle <- function(n_sites, bbox_diag_km, n_rows) {
  paste0(
    "Sites: ", n_sites,
    " | Rows: ", n_rows,
    " | BBox diag (km): ", round(bbox_diag_km, 1)
  )
}

safe_name <- function(x) {
  str_replace_all(x, "[^A-Za-z0-9]+", "_")
}

# Output directories -------------------------------------------------------
output_root <- here(plots_dir, "spatial_coverage")
output_species <- here(output_root, "top_species")
output_hp <- here(output_root, "top_host_pathogen")

dir.create(output_species, showWarnings = FALSE, recursive = TRUE)
dir.create(output_hp, showWarnings = FALSE, recursive = TRUE)

# 1. Species maps ----------------------------------------------------------
for (i in seq_len(nrow(selected_species))) {
  species_name <- selected_species$host_species[i]

  plot_data <- sites %>%
    filter(host_species == species_name)

  if (nrow(plot_data) == 0) {
    next
  }

  extent <- compute_extent(plot_data)
  subtitle <- format_subtitle(
    n_sites = selected_species$n_sites[i],
    bbox_diag_km = selected_species$bbox_diag_km[i],
    n_rows = selected_species$n_rows[i]
  )

  p <- create_inset_map(
    plot_data = plot_data,
    title = str_wrap(species_name, width = 40),
    subtitle = subtitle,
    inset_rect = extent,
    inset_pos = inset_pos
  )

  out_path <- here(output_species, paste0("map_", safe_name(species_name), ".png"))
  ggsave(out_path, p, width = 11, height = 8, dpi = 300)

}

# 2. Host–pathogen maps ----------------------------------------------------
if (nrow(selected_hp) > 0) {
  sites_hp <- sites %>%
    filter(!is.na(pathogen_species_cleaned)) %>%
    mutate(host_pathogen = paste(host_species, pathogen_species_cleaned, sep = " | "))

  for (i in seq_len(nrow(selected_hp))) {
    hp_name <- selected_hp$host_pathogen[i]

    plot_data <- sites_hp %>%
      filter(host_pathogen == hp_name)

    if (nrow(plot_data) == 0) {
      next
    }

    extent <- compute_extent(plot_data)
    subtitle <- format_subtitle(
      n_sites = selected_hp$n_sites[i],
      bbox_diag_km = selected_hp$bbox_diag_km[i],
      n_rows = selected_hp$n_rows[i]
    )

    p <- create_inset_map(
      plot_data = plot_data,
      title = str_wrap(hp_name, width = 45),
      subtitle = subtitle,
      inset_rect = extent,
      inset_pos = inset_pos
    )

    out_path <- here(output_hp, paste0("map_", safe_name(hp_name), ".png"))
    ggsave(out_path, p, width = 11, height = 8, dpi = 300)

  }
}

cat("\nSpatial coverage maps saved in:", output_root, "\n")
cat("Top species mapped:", nrow(selected_species), "\n")
cat("Top host-pathogen mapped:", nrow(selected_hp), "\n")
