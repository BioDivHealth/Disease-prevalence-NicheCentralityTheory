library(dplyr)
library(ggplot2)
library(maps)
library(sf)
library(tidyr)
library(scatterpie)
library(rnaturalearth)
library(ggspatial)
setwd("./trapping_data_scripts/")
top_20_df <- readRDS("../Data/clean_site_data.rds")

host_site_summary <- top_20_df %>%
  group_by(host_name) %>%
  summarise(unique_coordinates = n_distinct(decimalLatitude, decimalLongitude))

#########
# Barplot of number of locations per species

ggplot(host_site_summary, aes(x = reorder(host_name, unique_coordinates), y = unique_coordinates)) +
  geom_bar(stat = "identity", fill = "orangered", alpha = 0.7) +
  coord_flip() +
  labs(
    title = "Species list by unique coordinates and prevalence",
    x = "Species",
    y = "Number of Unique Coordinates"
  ) +
  theme(
    legend.position = "bottom",          # Adjust legend position
    plot.background = element_rect(fill = "white"),# Set panel background to white
  )  # Adjust legend position)


ggsave("../Results/Figures/n_unique_sites.png", dpi=500)

#########
# Maps of trapped rodents
###

### ALL ON ONE MAP

# Use the filtered unique_top_20_df directly
# Plot the coordinates on a world map with color representing different host_name groups
ggplot(data = top_20_df, aes(x = decimalLongitude, y = decimalLatitude, color = host_name)) +
  borders("world", colour = "gray85", fill = "gray80") +  # Add world map borders
  geom_point(size = 3) +  # Plot the points
  labs(title = "Locations by species",
       x = "Longitude",
       y = "Latitude",
       color = "Host Name") +
  theme_minimal() +
  theme(
    legend.position = "bottom",          # Adjust legend position
    plot.background = element_rect(fill = "white"),# Set panel background to white
    panel.grid = element_blank()         # Remove grid lines
  )  # Adjust legend position

ggsave("../Results/Figures/world_sites_map.png", dpi=500)


### SCATTERPIE WORLD MAP

# Load continent data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Function to find continent for a given coordinate
get_continent <- function(latitude, longitude) {
  point <- st_point(c(longitude, latitude)) %>%
    st_sfc(crs = st_crs(world))  # Create a spatial point with the same CRS as the continent data
  
  continent <- world %>%
    st_contains(point, sparse = FALSE) %>%
    apply(1, any) %>%
    which() %>%
    {if (length(.) > 0) world$continent[.] else NA}
  
  return(continent)
}
a

# Add continent information to the data
top_20_df <- top_20_df %>%
  rowwise() %>%
  mutate(continent = get_continent(decimalLatitude, decimalLongitude)) %>%
  ungroup()

# Summarize the data by continent and host_name
continent_summary <- top_20_df %>%
  group_by(continent, host_name) %>%
  summarize(tested = n(), .groups = 'drop') %>%
  pivot_wider(names_from = host_name, values_from = tested, values_fill = 0)

# Create a dataframe with the centroids of continents for plotting pie charts
continent_centroids <- world %>%
  group_by(continent) %>%
  summarize(geometry = st_union(geometry), .groups = 'drop') %>%
  st_centroid() %>%
  filter(continent %in% unique(continent_summary$continent))

# Merge centroid coordinates with the summary data
continent_pie_data <- left_join(continent_centroids, continent_summary, by = "continent")

# Plot the map with pie charts
ggplot() +
  geom_sf(data = world, fill = "orange", color = "gray85", size = 0.5, alpha = 0.5) +  # World map
  geom_scatterpie(data = continent_pie_data,
                  aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2]),
                  cols = names(continent_summary)[-1], pie_scale = 0.1) +  # Add pie charts
  coord_sf() +
  labs(title = "Proportion of Hosts by Continent",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

##### UKRAINE

# Filter the data for Apodemus flavicolus
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
  coord_fixed(xlim = c(23.99, 24.10), ylim = c(51.195, 51.275)) +  # Zoom into Volyn Oblast's coordinates
  theme_minimal()


###just one transect

# Plot the coordinates on a map of Volyn Oblast, Ukraine
coordinates_sf <- st_as_sf(filtered_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

ggplot() +
  geom_sf(data = coordinates_sf, color = "blue", size = 3) +  # Plot the points
  annotation_scale(location = "bl", width_hint = 0.5, style = "bar", unit_category = "metric") +  # Add scale bar in km
  coord_sf(xlim = c(24.08, 24.095), ylim = c(51.198, 51.208), expand = FALSE) +  # Zoom in on a specific region
  labs(title = paste("Locations for", host_name_group, "in Volyn Oblast, Ukraine"),
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(       
    plot.background = element_rect(fill = "lightgrey")
  )

ggsave("../Results/Figures/Apodemus_flavicolus_Ukraine.png", dpi=500)

#### ASIA

# Filter the data for Suncus murinus
host_name_group <- "Suncus murinus"  # Change to your desired host_name
filtered_df <- top_20_df %>%
  filter(host_name == host_name_group)

# Plot the coordinates on a map of Asia
ggplot(data = filtered_df, aes(x = decimalLongitude, y = decimalLatitude)) +
  borders("world", colour = "grey30", fill = "lightgrey") +  # Add map borders for Asia
  geom_point(color = "red", size = 3) +  # Plot the points
  labs(title = paste("Locations for Host:", host_name_group, "in Asia"),
       x = "Longitude",
       y = "Latitude") +
  coord_fixed(xlim = c(25, 135), ylim = c(-25, 35)) +  # Adjust limits for Asia
  theme_minimal()

### now overlay shapefile

# Load the Shapefile (replace 'path_to_shapefile' with the actual path to your .shp file)
shapefile_path <- "../Data//iucn_data/iucn_data_suncus_murinus/data_0.shp" 
polygon_data <- st_read(shapefile_path)

#### URUGUAY

# Plot the map of Uruguay with the Shapefile overlay
ggplot() +
  borders("world", colour = "grey30", fill = "lightgrey") +  # Add map borders for Asia
  geom_sf(data = polygon_data, color = "darkorange", fill = "orange", alpha = 0.5, size = 0.5) +  # Overlay the Shapefile polygon
  geom_point(data = filtered_df, aes(x = decimalLongitude, y = decimalLatitude), color = "red", size = 3, shape = 18) +  # Plot points
  labs(title = paste(host_name_group, "trapping locations and species range"),
       x = "Longitude",
       y = "Latitude") +
  coord_sf(xlim = c(25, 135), ylim = c(-25, 35)) +  # Adjust limits for Asia
  theme_minimal() +
  theme(
    legend.position = "bottom",          # Adjust legend position
    plot.background = element_rect(fill = "white"),# Set panel background to white
    panel.grid = element_blank()         # Remove grid lines
  )

ggsave("../Results/Figures/Suncus_murinus_map.png", dpi=500)

# plot mus musculus sites in uruguay
host_name_group <- "Mus musculus"
filtered_df <- top_20_df %>%
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
shapefile_path <- "../Data//iucn_data/iucn_data_mus_musculus/data_0.shp"   # Adjust this path
polygon_data <- st_read(shapefile_path)


# Plot the map of Uruguay with the Shapefile overlay
ggplot() +
  borders("world", colour = "grey30", fill = "lightgrey", size = 1) +  # Base map of Uruguay
  geom_sf(data = polygon_data, color = "darkgreen", fill = "springgreen3", alpha = 0.5, linewidth = 1) +  # Overlay the Shapefile polygon
  geom_point(data = filtered_df, aes(x = decimalLongitude, y = decimalLatitude), color = "red", size = 5, shape = 18) +  # Plot points
  labs(title = paste("Trapping locations for Mus musculus"),
       x = "Longitude",
       y = "Latitude") +
  coord_sf(xlim = c(-64, -50), ylim = c(-38, -27)) +  # Zoom into Uruguay's coordinates
  theme_minimal() +
  theme(
    legend.position = "bottom",          # Adjust legend position
    plot.background = element_rect(fill = "white"),# Set panel background to white
    panel.grid = element_blank()         # Remove grid lines
  )

ggsave("../Results/Figures/Mus_musculus_uruguay.png", dpi=500)

#### TEXAS

# plot Baiomys taylori sites in texas
host_name_group <- "Baiomys taylori"
filtered_df <- top_20_df %>%
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
