library(dplyr)
library(ggplot2)
library(maps)


df <- readRDS("./Data/host_path_wide_SDM.rds")

unique_sites_df <- df %>%
  count(host_name, sort = TRUE) %>%        # Count occurrences of each host_name and sort by frequency
  slice(1:20) %>%                          # Keep only the top 20 most frequent host_name values
  select(host_name) %>%                    # Select only the host_name column
  inner_join(df, by = "host_name") %>%     # Join back with the original dataframe to filter rows
  distinct(host_name, start_date, locality) %>%  # Keep unique combinations of host_name, start_date, locality
  arrange(host_name) %>%                   # Sort by host_name in ascending order
  left_join(df %>%                         # Join with coordinate data
              distinct(locality, decimalLatitude, decimalLongitude), # Keep first pair for each locality
            by = "locality")

saveRDS(unique_sites_df, file = "./Data/clean_sdm_data.rds")
