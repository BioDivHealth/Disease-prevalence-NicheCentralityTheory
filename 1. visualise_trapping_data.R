library(dplyr)
library(ggplot2)
library(maps)
library(sf)

top_20_df <- readRDS("./Data/clean_site_data.rds")

host_site_summary <- top_20_df %>%
  group_by(host_name) %>%
  summarise(unique_coordinates = n_distinct(decimalLatitude, decimalLongitude))

#########
# Barplot of number of locations per species

ggplot(host_site_summary, aes(x = reorder(host_name, unique_coordinates), y = unique_coordinates)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(
    title = "Species list by unique coordinates and prevalence",
    x = "Species",
    y = "Number of Unique Coordinates"
  ) +
  theme_minimal()

#########
# Maps of trapped rodents
###

# Filter the data for a specific host_name
host_name_group <- "Apodemus flavicolus"  # Change to your desired host_name
filtered_df <- top_20_df %>%
  filter(host_name == host_name_group)

# Plot the coordinates on a map of Volyn Oblast, Ukraine
ggplot(data = filtered_df, aes(x = decimalLongitude, y = decimalLatitude)) +
  borders("world", region = "Ukraine", colour = "gray85", fill = "gray80") +  # Add map borders for Ukraine
  geom_point(color = "red", size = 3) +  # Plot the points
  labs(title = paste("Locations for Host:", host_name_group, "in Volyn Oblast, Ukraine"),
       x = "Longitude",
       y = "Latitude") +
  coord_fixed(xlim = c(23.96, 24.13), ylim = c(51.175, 51.3)) +  # Zoom into Volyn Oblast's coordinates
  theme_minimal()



# Filter the data for a specific host_name
host_name_group <- "Suncus murinus"  # Change to your desired host_name
filtered_df <- top_20_df %>%
  filter(host_name == host_name_group)

# Plot the coordinates on a map of Asia
ggplot(data = filtered_df, aes(x = decimalLongitude, y = decimalLatitude)) +
  borders("world", colour = "gray85", fill = "gray80") +  # Add map borders for Asia
  geom_point(color = "red", size = 3) +  # Plot the points
  labs(title = paste("Locations for Host:", host_name_group, "in Asia"),
       x = "Longitude",
       y = "Latitude") +
  coord_fixed(xlim = c(60, 135), ylim = c(-10, 38)) +  # Adjust limits for Asia
  theme_minimal()





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

#### with shapefile of species range

# Load the Shapefile (replace 'path_to_shapefile' with the actual path to your .shp file)
shapefile_path <- "./Data/data_0.shp"  # Adjust this path
polygon_data <- st_read(shapefile_path)



# Plot the map of Uruguay with the Shapefile overlay
ggplot() +
  borders("world", region = "Uruguay", colour = "gray85", fill = "gray80") +  # Base map of Uruguay
  geom_sf(data = polygon_data, fill = NA, color = "blue", size = 0.5) +  # Overlay the Shapefile polygon
  geom_point(data = filtered_df, aes(x = decimalLongitude, y = decimalLatitude), color = "red", size = 3) +  # Plot points
  labs(title = paste("Locations for Host:", host_name_group, "with Shapefile Overlay"),
       x = "Longitude",
       y = "Latitude") +
  coord_sf(xlim = c(-59, -53), ylim = c(-35, -30)) +  # Zoom into Uruguay's coordinates
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
