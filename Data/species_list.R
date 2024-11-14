library(dplyr)
library(tidyr)
library(here)

write_rds(host_path_wide, file="./data/Species_list/host_path_wide.rds")

host_path_wide = readRDS("./data/Species_list/host_path_wide.rds")

# Treat coordinates as character
host_path_wide$decimalLongitude <- as.character(host_path_wide$decimalLongitude)
host_path_wide$decimalLatitude <- as.character(host_path_wide$decimalLatitude)

# Create unique coordinate pairs
host_path_wide <- host_path_wide %>%
  mutate(coordinate_pair = paste(decimalLongitude, decimalLatitude, sep = ", "))

# Summarize data for each host species including the count of unique coordinate pairs
host_summary <- host_path_wide %>%
  group_by(host_name) %>%
  summarise(
    tested = sum(as.numeric(tested), na.rm = TRUE),
    positives = sum(as.numeric(positive), na.rm = TRUE),
    unique_coordinates = n_distinct(coordinate_pair),
    all_coordinates = list(unique(coordinate_pair)) 
  )

# Expand the list of coordinates into individual columns
host_summary_coordinates <- host_summary %>%
  unnest_wider(all_coordinates, names_sep = "_coord_")
