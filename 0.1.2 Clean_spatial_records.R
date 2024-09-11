rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored
td<-tempdir()

#
#
#/\/\/\/\/\/\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/////\/\/\/\/\/\//\/\/\/>//\/\///\
##                                                                          \//\/|/      
#### Spatial data gathering and preparation for the Range-Disease analysis   
##                                                                          \//\/\/\
#/\/\/\/\/\/\/\/\\/\/\/\/\/\/\//\/\/\/\/\/\/////\/\/\/\/\/\///\/\\//\/\/\///\
#
# 0. load the needed libraries----
list.of.packages<-c("tidyverse","doParallel","foreach","rstudioapi","colorspace","readxl",
                    "data.table","parallel","rredlist")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# project paramteres
crs_p <- "EPSG:4326"


# Load the needed functions
functions<-"./Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 1. Load the species information----
data_route <- "./Data/Species_list" ; data_route %>% list.files(pattern=".csv$",full.names = TRUE)

sp_list <- data_route %>% list.files(pattern="Species_analysis.csv$",full.names = TRUE) %>% read.csv()
IUCN_info <-lapply("./Data/IUCN_info/sp_info" %>% list.files(recursive=TRUE,full.names=TRUE),read.csv) %>% rbindlist()

# species points----
sp_points <- "./Data/Sp_info/raw_records" %>% list.files(recursive=TRUE,full.names=TRUE)

# 1.b Species ranges from IUCN-RedList----
# We have downloaded the list of range polygons for the mammals (https://www.iucnredlist.org/resources/spatial-data-download)
paste(td,"iucn_pols",sep="/") %>% dir.create(showWarnings = FALSE,recursive = TRUE)
range_route <-"./Data/Species_Ranges" %>% list.files(pattern=".zip$",full.names = TRUE)

range_route %>% unzip(exdir=paste(td,"iucn_pols",sep="/")) # Unzip the spatial information into the temporal folder
range_sp <- paste(td,"iucn_pols",sep="/") %>% list.files(pattern = ".dbf$",recursive = TRUE,full.names = TRUE) %>% st_read()

range_sp <- range_sp %>% filter(sci_name %in% sp_list$IUCN_name)
range_sp %>% st_write(paste("./Data/Species_ranges","Range_analysis.gpkg",sep="/"))

range_sp <-paste("./Data/Species_ranges","Range_analysis.gpkg",sep="/") %>% st_read()
paste(td,"iucn_pols",sep="/") %>% unlink(recursive = T)

# 2. Clean the spatial records---- 


# 1.d. Check the IUCN Red List Spatial information ----
# Look for the polygon

# collect and save
# RM the unzipped data

# 2. Clean species records ----




# 3. Collect the spatial information
species_occ_list
species_polygons_list









# 1.c Check for spatial data for the species----
pol_route <- "D:/Data/Spatial information/IUCN spatial data/All polygons"
range_pols <- data.frame(route=pol_route %>% list.files(pattern = ".shp$",recursive=TRUE,full.names = TRUE),
                         species=pol_route %>% list.files(pattern = ".shp$",recursive=TRUE,full.names = TRUE) %>% 
                           basename() %>% gsub(pattern=".shp$",replacement=""))

# 2. Download the spatial information for the species ----
for(i in 417:length(species)){
  
  # if(species[i] %in% range_pols$species){
  #                     y_range<-range_pols %>% filter(species == species[i]) %>% 
  #                                   dplyr::select("route") %>% sf::read_sf() %>% 
  #                                                 sf::st_as_sfc() %>% sf::st_convex_hull() %>% sf::st_as_text()
  #                     
  #                 }else{
  #               y_range<-NULL
  #                   }
  
  try(Spatial_spp(sci_sp = species[i], start_date = 2019),silent=FALSE)# range_sp=y_range)
}





unlink(td,recursive=TRUE)
