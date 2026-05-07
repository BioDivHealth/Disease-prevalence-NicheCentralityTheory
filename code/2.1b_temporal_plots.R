# 2.1b_temporal_plots
# Goal: plotting + figure export for temporal audit outputs from
# 2.1a_temporal_candidates.R.

library(pacman)
p_load(tidyverse, lubridate, here, maps, patchwork)

# Inputs written by 2.1a_temporal_candidates.R
site_mode <- "exact" # "exact" or "rounded"
coord_round_digits <- 2
mode_tag <- if (site_mode == "exact") "exact" else paste0("rounded", coord_round_digits)
analysis_dir <- here("Results", "analysis_metadata", mode_tag)
plots_dir <- here("Results", "plots", mode_tag)

path_dat_final <- here("Data", paste0("dat_with_modis_metadata_", mode_tag, ".rds"))
path_dat_year_temporal <- here(analysis_dir, "dat_year_temporal.rds")
path_candidates_species <- here(analysis_dir, "temporal_candidates_species.rds")
path_dat_year_temporal_hp <- here(analysis_dir, "dat_year_temporal_host_pathogen.rds")
path_candidates_hp <- here(analysis_dir, "temporal_candidates_host_pathogen.rds")

stopifnot(file.exists(path_dat_year_temporal))
stopifnot(file.exists(path_candidates_species))
stopifnot(file.exists(path_dat_year_temporal_hp))
stopifnot(file.exists(path_candidates_hp))

# Load
# dat_year_temporal contains the filtered (<= 1 year) dataset + event_date/short_window
# dat_final contains MODIS-derived metadata (optional; required for phenology QA plot)
dat_year_temporal <- readRDS(path_dat_year_temporal)
temporal_candidates <- readRDS(path_candidates_species)
dat_year_temporal_hp <- readRDS(path_dat_year_temporal_hp)
temporal_candidates_hp <- readRDS(path_candidates_hp)

# Output folder
out_dir <- plots_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Histogram: time between start and end dates ----
p_zoom_count <- ggplot(dat_year_temporal, aes(x = days_diff)) +
  geom_histogram(
    binwidth = 3, boundary = 0,
    colour = "white", linewidth = 0.3, fill = "coral3"
  ) +
  coord_cartesian(xlim = c(0, 370)) +
  scale_x_continuous(
    breaks = seq(0, 370, by = 30),
    minor_breaks = seq(0, 370, by = 3),
    expand = c(0, 0)
  ) +
  labs(
    title = "Time between start and end dates (0–365 days)",
    subtitle = "Histogram (3-day bins)",
    x = "Days (end − start)",
    y = "Count"
  ) +
  theme_classic(base_size = 13)

ggsave(
  filename = here(out_dir, "temporal_window_histogram.png"),
  plot = p_zoom_count,
  width = 8,
  height = 4.5,
  dpi = 300
)

# 2. Phenology QA plot (requires dat_with_modis_metadata.rds) ----
if (file.exists(path_dat_final)) {
  dat_final <- readRDS(path_dat_final)

  p_pheno_qa <- ggplot(
    dat_final %>% filter(pheno_source != "Missing"),
    aes(x = days_since_midgreenup, fill = pheno_source)
  ) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = "Distribution of Sampling Relative to Mid-Greenup",
      subtitle = "0 = Peak greening rate; positive = post-greenup sampling",
      x = "Days since Mid-Greenup",
      y = "Count"
    ) +
    theme_minimal()

  ggsave(
    filename = here(out_dir, "phenology_days_since_midgreenup.png"),
    plot = p_pheno_qa,
    width = 8,
    height = 4.5,
    dpi = 300
  )
}

# 3. Temporal coverage plot (species) ----
max_period <- 15

top_10_species <- temporal_candidates$host_species[seq_len(min(10, nrow(temporal_candidates)))]

p_temporal_coverage <- dat_year_temporal %>%
  filter(host_species %in% top_10_species, short_window) %>%
  ggplot(aes(x = event_date, y = host_species)) +
  geom_jitter(alpha = 0.4, size = 0.6, color = "darkred", width = 0, height = 0.15) +
  labs(
    title = "General Temporal Sampling Density: Top 10 Species",
    subtitle = paste0("Short windows (<=", max_period, "d) shown using midpoint/event dates"),
    x = "Date",
    y = "Host Species"
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme_minimal()

p_temporal_coverage

ggsave(
  filename = here(out_dir, "temporal_coverage_top10_species.png"),
  plot = p_temporal_coverage,
  width = 10,
  height = 4.8,
  dpi = 300
)

# 4. Temporal coverage plot (host x pathogen) ----
top_10_hp <- temporal_candidates_hp$host_pathogen[seq_len(min(10, nrow(temporal_candidates_hp)))]

p_temporal_coverage_hp <- dat_year_temporal_hp %>%
  filter(host_pathogen %in% top_10_hp, short_window) %>%
  ggplot(aes(x = event_date, y = host_pathogen)) +
  geom_jitter(alpha = 0.4, size = 0.6, color = "darkgreen", width = 0, height = 0.18) +
  labs(
    title = "Temporal Sampling Density: Top 10 Host × Pathogen Combos",
    subtitle = paste0("Short windows (<=", max_period, "d) shown using midpoint/event dates"),
    x = "Date",
    y = "Host | Pathogen"
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))

ggsave(
  filename = here(out_dir, "temporal_coverage_top10_host_pathogen.png"),
  plot = p_temporal_coverage_hp,
  width = 10,
  height = 5.5,
  dpi = 300
)

# 5. Faceted temporal audit by pathogen (top combos) ----
top_10_data <- dat_year_temporal_hp %>%
  filter(host_pathogen %in% top_10_hp, short_window)

adaptive_breaks <- function(lims) {
  if (any(is.na(lims)) || length(lims) < 2) return(NULL)
  span_days <- as.numeric(diff(lims))
  span_years <- span_days / 365.25
  if (is.na(span_years) || span_years <= 0) return(NULL)

  if (span_days < 30) {
    seq(from = lims[1], to = lims[2], by = "3 days")
  } else if (span_years < 1.5) {
    start <- lubridate::floor_date(lims[1], "month")
    end <- lubridate::ceiling_date(lims[2], "month")
    unique(seq(start, end, by = "1 month"))
  } else if (span_years <= 4) {
    start <- lubridate::floor_date(lims[1], "year")
    end <- lubridate::ceiling_date(lims[2], "year")
    unique(seq(start, end, by = "1 year"))
  } else {
    start_year <- lubridate::year(lims[1])
    end_year <- lubridate::year(lims[2])
    years_seq <- seq(start_year, end_year + 1, by = 2)
    lubridate::make_date(year = years_seq, month = 1, day = 1)
  }
}

adaptive_labels <- function(dates) {
  if (length(dates) == 0) return(character(0))
  if (all(is.na(dates))) return(rep(NA_character_, length(dates)))

  valid_dates <- dates[!is.na(dates)]
  if (length(valid_dates) < 2) {
    return(format(dates, "%b %d"))
  }

  span_days <- as.numeric(max(valid_dates) - min(valid_dates))
  span_years <- span_days / 365.25

  if (span_days < 30) {
    format(dates, "%b %d")
  } else if (span_years < 1.5) {
    format(dates, "%b %Y")
  } else {
    format(dates, "%Y")
  }
}

p_temporal_pathogen_facets <- ggplot(top_10_data, aes(x = event_date, y = host_species)) +
  geom_jitter(alpha = 0.5, size = 0.8, color = "darkblue", width = 0, height = 0.2) +
  facet_wrap(~pathogen_species_cleaned, scales = "free", ncol = 1, axes = "all_x") +
  labs(
    title = "Temporal Sampling Density Faceted by Pathogen (Free Scales)",
    subtitle = paste0("Showing host species for top-sampled pathogens (Short windows <=", max_period, "d)"),
    x = "Date",
    y = "Host Species"
  ) +
  scale_x_date(breaks = adaptive_breaks, labels = adaptive_labels) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 8),
    panel.spacing = unit(1, "lines")
  )

ggsave(
  filename = here(out_dir, "temporal_facets_by_pathogen_top10.png"),
  plot = p_temporal_pathogen_facets,
  width = 10,
  height = 10,
  dpi = 300
)

# 6. Explicit date ranges for top combos ----
hp_ranges <- top_10_data %>%
  group_by(host_pathogen, host_species, pathogen_species_cleaned) %>%
  summarise(
    min_date = min(event_date, na.rm = TRUE),
    max_date = max(event_date, na.rm = TRUE),
    .groups = "drop"
  )

p_hp_ranges <- ggplot(hp_ranges, aes(y = host_species)) +
  geom_segment(
    aes(x = min_date, xend = max_date, y = host_species, yend = host_species),
    color = "grey70",
    linewidth = 1.5,
    alpha = 0.5
  ) +
  geom_jitter(
    data = top_10_data,
    aes(x = event_date, y = host_species),
    alpha = 0.6,
    size = 0.8,
    color = "darkblue",
    width = 0,
    height = 0.1
  ) +
  facet_wrap(~pathogen_species_cleaned, scales = "free", ncol = 1, axes = "all_x") +
  labs(
    title = "Sampling Date Ranges by Pathogen (Free Scales)",
    subtitle = "Lines show the full span from first to last sample; dots show individual events",
    x = "Date",
    y = "Host Species"
  ) +
  scale_x_date(breaks = adaptive_breaks, labels = adaptive_labels) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 8),
    panel.spacing = unit(1, "lines")
  )

ggsave(
  filename = here(out_dir, "temporal_date_ranges_by_pathogen_top10.png"),
  plot = p_hp_ranges,
  width = 10,
  height = 10,
  dpi = 300
)

# 7. Spatial plots (optional) ----
# These are intentionally opt-in because they can be slow and/or produce many files.
make_spatial_maps <- FALSE
make_all_hp_maps <- FALSE

if (make_spatial_maps) {
  p_load(grid)

  world <- map_data("world")

  create_inset_map <- function(hp_name, data, world_data) {
    hp_data <- data %>% filter(host_pathogen == hp_name)
    year_limits <- range(hp_data$year, na.rm = TRUE)

    lon_range <- range(hp_data$longitude, na.rm = TRUE)
    lat_range <- range(hp_data$latitude, na.rm = TRUE)

    lon_buffer <- max(diff(lon_range) * 0.3, 1.5)
    lat_buffer <- max(diff(lat_range) * 0.3, 1.5)

    xlims <- c(max(-180, lon_range[1] - lon_buffer), min(180, lon_range[2] + lon_buffer))
    ylims <- c(max(-90, lat_range[1] - lat_buffer), min(90, lat_range[2] + lat_buffer))

    p_main <- ggplot() +
      geom_polygon(
        data = world_data,
        aes(x = long, y = lat, group = group),
        fill = "grey92",
        color = "white",
        linewidth = 0.1
      ) +
      geom_point(
        data = hp_data,
        aes(x = longitude, y = latitude, color = year),
        size = 1.8,
        alpha = 0.4,
        position = position_jitter(width = 0.1, height = 0.1)
      ) +
      coord_quickmap(xlim = xlims, ylim = ylims) +
      labs(title = str_wrap(hp_name, width = 40)) +
      scale_color_viridis_c(
        option = "D",
        limits = year_limits,
        guide = guide_colorbar(direction = "horizontal", barwidth = grid::unit(7, "cm"))
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.background = element_rect(fill = "aliceblue", color = "grey80"),
        plot.title = element_text(face = "bold", size = 9, margin = margin(b = 2)),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 8),
        panel.grid.major = element_line(color = "white", linewidth = 0.2)
      )

    p_inset <- ggplot() +
      geom_polygon(
        data = world_data,
        aes(x = long, y = lat, group = group),
        fill = "grey70",
        color = "white",
        linewidth = 0.05
      ) +
      geom_rect(
        aes(xmin = xlims[1], xmax = xlims[2], ymin = ylims[1], ymax = ylims[2]),
        color = "red",
        fill = "red",
        alpha = 0.2,
        linewidth = 0.2
      ) +
      coord_quickmap() +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        plot.margin = margin(0, 0, 0, 0)
      )

    p_main + inset_element(p_inset, left = 0.5, bottom = 0.8, right = 1, top = 1, align_to = "panel")
  }

  map_list <- lapply(top_10_hp, create_inset_map, data = top_10_data, world_data = world)

  p_spatial_dist_inset <- wrap_plots(map_list, ncol = 3) +
    plot_annotation(
      title = "Spatial Distribution Host-Pathogen Combos",
      subtitle = "Main map: Zoomed to occurrences | Top-right: Global context",
      theme = theme(
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 10, color = "grey30")
      )
    )

  ggsave(
    here(out_dir, "top_hp_spatial_distribution.png"),
    p_spatial_dist_inset,
    width = 12,
    height = 12,
    dpi = 300,
    scale = 0.7
  )
}

if (make_all_hp_maps) {
  p_load(grid)

  world <- map_data("world")
  all_hp_data <- dat_year_temporal_hp
  all_hp_combos <- sort(unique(all_hp_data$host_pathogen))

  out_all <- here(out_dir, "all_hp_maps")
  dir.create(out_all, showWarnings = FALSE, recursive = TRUE)

  create_inset_map <- function(hp_name, data, world_data) {
    hp_data <- data %>% filter(host_pathogen == hp_name)
    year_limits <- range(hp_data$year, na.rm = TRUE)

    lon_range <- range(hp_data$longitude, na.rm = TRUE)
    lat_range <- range(hp_data$latitude, na.rm = TRUE)

    lon_buffer <- max(diff(lon_range) * 0.3, 1.5)
    lat_buffer <- max(diff(lat_range) * 0.3, 1.5)

    xlims <- c(max(-180, lon_range[1] - lon_buffer), min(180, lon_range[2] + lon_buffer))
    ylims <- c(max(-90, lat_range[1] - lat_buffer), min(90, lat_range[2] + lat_buffer))

    ggplot() +
      geom_polygon(
        data = world_data,
        aes(x = long, y = lat, group = group),
        fill = "grey92",
        color = "white",
        linewidth = 0.1
      ) +
      geom_point(
        data = hp_data,
        aes(x = longitude, y = latitude, color = year),
        size = 1.2,
        alpha = 0.4,
        position = position_jitter(width = 0.1, height = 0.1)
      ) +
      coord_quickmap(xlim = xlims, ylim = ylims) +
      labs(title = str_wrap(hp_name, width = 40)) +
      scale_color_viridis_c(option = "D", limits = year_limits) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "none",
        axis.title = element_blank()
      )
  }

  for (i in seq_along(all_hp_combos)) {
    hp_name <- all_hp_combos[i]
    hp_safe <- str_replace_all(hp_name, "[^A-Za-z0-9]+", "_")
    out_path <- here(out_dir, "all_hp_maps", paste0("map_", hp_safe, ".png"))

    tryCatch(
      {
        p <- create_inset_map(hp_name, data = all_hp_data, world_data = world)
        ggsave(out_path, p, width = 8, height = 7, dpi = 150)
        if (i %% 10 == 0) cat("Processed", i, "/", length(all_hp_combos), "maps...\n")
      },
      error = function(e) {
        message(paste("Error generating map for", hp_name, ":", e$message))
      }
    )
  }

  cat("Batch generation complete. Files located in Results/plots/all_hp_maps/\n")
}
