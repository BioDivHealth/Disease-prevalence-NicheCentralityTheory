coordinates <- readRDS("github/rodent-zoonoses/data/host_path_wide.rds") #Add host_path to github

coordinates <- coordinates %>%
  mutate(associatedTaxa = case_when(
    associatedTaxa == "Akodon azare" ~ "Akodon azarae",
    associatedTaxa == "brush mice" ~ "Peromyscus boylii",
    associatedTaxa == "Clethrionomys glareolus" ~ "Myodes glareolus",
    associatedTaxa == "Clethryonomys glareolus" ~ "Myodes glareolus",
    associatedTaxa == "Oligoryzomys favescens" ~ "Oligoryzomys flavescens",
    associatedTaxa == "P. californicus" ~ "Peromyscus californicus",
    associatedTaxa == "dusky-footed" ~ "Neotoma fuscipes",
    TRUE ~ associatedTaxa  # Retain original name if no match
  ))


# Create spatial object from Data.Frame
data<-sf::st_as_sf(coordinates,coords=c("decimalLongitude","decimalLatitude"))
plot(data %>% st_geometry())

# Look up for duplicates
data<-data %>% mutate(data %>% sf::st_coordinates())

if((data$rodent_record_id %>% unique() %>% length()) == nrow(data)){
  print("No duplicated records (rodent_record_id)")
  
}else{
  print("there are duplicated records")
  index_dupli<-data$rodent_record_id[duplicated(data$rodent_record_id,fromLast=TRUE)]
  duplicated<-data[data$rodent_record_id %in% index_dupli,]
  
  
}

pt_sp <- paste(data$rodent_record_id)