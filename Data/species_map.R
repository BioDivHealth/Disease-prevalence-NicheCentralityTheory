# Load necessary libraries
library(leaflet)
library(dplyr)  # For data manipulation
library(tidyr)  # For reshaping the data

# Step 1: Reshape the data to gather coordinate columns into a single column
top_20_long <- top_20 %>%
  pivot_longer(
    cols = starts_with("all_coordinates_coord_"),  # Selects coordinate columns
    names_to = "coord_number", 
    values_to = "coords"
  ) 

# Step 2: Remove rows where 'coords' is NA
top_20_long <- top_20_long %>%
  filter(!is.na(coords)) 

# Step 3: Split 'coords' into longitude and latitude
top_20_long <- top_20_long %>%
  mutate(
    longitude = as.numeric(sapply(strsplit(coords, ","), `[`, 1)),
    latitude = as.numeric(sapply(strsplit(coords, ","), `[`, 2))
  )

# Step 4: Plot with a natural color map (using Esri World Imagery or CartoDB Positron)
m <- leaflet() %>% 
  addTiles() %>%  # Default tiles for zoom levels
  setView(lng = mean(top_20_long$longitude, na.rm = TRUE), 
          lat = mean(top_20_long$latitude, na.rm = TRUE), 
          zoom = 5) %>% 
  addProviderTiles("Esri.WorldImagery") %>%  # Change to green/brown-themed map
  addCircleMarkers(
    lng = top_20_long$longitude, 
    lat = top_20_long$latitude, 
    popup = paste("Host:", top_20_long$host_name, "<br>Tested:", top_20_long$tested),
    radius = log(top_20_long$tested + 1) * 0.5,  # Adjust dot size
    color = "red", 
    fillOpacity = 0.8
  )

# View the map
m
