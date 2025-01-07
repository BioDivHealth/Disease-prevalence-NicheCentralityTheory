library(dplyr)
library(ggplot2)


df <- readRDS("./Data/trapping_data_NC.rds")

df <- df %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))

# Summarize the 'count' for each unique combination of coordinates
summed <- df %>%
  group_by(locality, decimalLatitude, decimalLongitude) %>%
  summarise(total_tested = sum(as.numeric(tested)),
            total_positive = sum(as.numeric(positive)), 
            .groups = "drop")

# Get unique combinations and add the corresponding coordinate pairs and summed 'count' for each locality
unique_sites_df <- df %>%
  distinct(host_name, locality) %>%  # Keep unique combinations of host_name, start_date, locality
  left_join(summed, by = c("locality"))   # Join with summarized 'count' data

### Keep only the top 20 species

# Filter to keep only the 20 most common host_name values
top_20_host_names <- unique_df %>%
  count(host_name, sort = TRUE) %>%       # Count occurrences of each host_name and sort by frequency
  slice(1:20) %>%                         # Keep only the top 20 most frequent host_name values
  select(host_name)                       # Select only the host_name column

# Filter the original dataframe to include only rows with the top 20 host_names
unique_sites_top_20 <- unique_sites_df %>%
  semi_join(top_20_host_names, by = "host_name")  # Keep only rows with host_name in top 20


saveRDS(unique_sites_top_20, file = "./Data/clean_site_data.rds")
