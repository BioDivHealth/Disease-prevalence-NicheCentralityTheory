library(dplyr)
library(tidyr)
library(here)
library(ggplot2)

host_path_wide = readRDS("./data/host_path_wide_SDM.rds")

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

# Remove rows with "spp" in host_names
species_list <- host_summary_coordinates[!grepl("spp", host_summary_coordinates$host_name, ignore.case = TRUE), ]

# Add prevalence column
species_list$prevalence <- species_list$positives / species_list$tested

# Order by higher number of positives and by unique_coordinates
# Exclude rows with unique_coordinates = 1, then order and filter top 20
top_20 <- species_list %>%
  filter(unique_coordinates > 1) %>%
  arrange(desc(positives), desc(unique_coordinates)) %>%
  slice_head(n = 20) 

write.csv(top_20, "data/top_20.csv", row.names = FALSE)

# Barplot
ggplot(top_20, aes(x = reorder(host_name, -unique_coordinates), y = unique_coordinates)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = sprintf("%.2f", prevalence), y = unique_coordinates + 0.5), 
            color = "blue", size = 3, hjust = 0) +
  coord_flip() +
  labs(
    title = "Species list by unique coordinates and prevalence",
    x = "Species",
    y = "Number of Unique Coordinates"
  ) +
  theme_minimal()

