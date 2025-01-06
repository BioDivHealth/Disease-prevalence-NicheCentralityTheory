library(dplyr)
library(ggplot2)
library(maps)

df <- readRDS("./Data/clean_sdm_data.rds")

# plot mus musculus sites in uruguay
host_name_group <- "Mus musculus"
filtered_df <- unique_sites_df %>%
  filter(host_name == host_name_group)

ggplot(data = filtered_df, aes(x = decimalLongitude, y = decimalLatitude)) +
  borders("world", region = "Uruguay", colour = "gray85", fill = "gray80") +  # Add map borders for Uruguay
  geom_point(color = "red", size = 3) +  # Plot the points
  labs(title = paste("Locations for", host_name_group, "in Uruguay"),
       x = "Longitude",
       y = "Latitude") +
  coord_fixed(xlim = c(-59, -53), ylim = c(-35, -30)) +  # Zoom into Uruguay's coordinates
  theme_minimal()

# plot mus musculus sites in uruguay
host_name_group <- "Baiomys taylori"
filtered_df <- unique_sites_df %>%
  filter(host_name == host_name_group)

# Plot the coordinates on a map of Texas
ggplot(data = filtered_df, aes(x = decimalLongitude, y = decimalLatitude)) +
  borders("state", region = "texas", colour = "gray85", fill = "gray80") +  # Add map borders for Texas
  geom_point(color = "red", size = 3) +  # Plot the points
  labs(title = paste("Locations for Host:", host_name_group, "in Texas"),
       x = "Longitude",
       y = "Latitude") +
  coord_fixed(xlim = c(-100, -99), ylim = c(28, 28.5)) +  # Zoom into Texas's coordinates
  theme_minimal()
