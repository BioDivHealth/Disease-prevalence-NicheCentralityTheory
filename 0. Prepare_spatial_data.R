rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored
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

# Load the needed functions
functions<-"./Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 1. Load the species information----
data_route <- "./Data/Species_list" ; data_route %>% list.files(pattern=".xlsx$")
sp_list <- read_xlsx(data_route %>% list.files(pattern=".xlsx$",full.names = TRUE),sheet=1)

# Check the IUCN API and version
# IUCN API token
options(iucn_redlist_key="eb704359f6ea22d50235efebf7f4a2f5f843a66fa951cf6a376df71cf7268986")
RL.version<-rredlist::rl_version() ; print(paste("RedList Verion",RL.version))

# 1.b Download the IUCN species information----
export_route <- "./Data/IUCN_info" ; export_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)# output route for the IUCN information

# Check species names
sp_names <- lapply(sp_list$Species,function(x) retrieve_syns(spp_name=x)$TaxDat) %>% rbindlist()
sp_analysis <- sp_names %>% filter(!is.na(IUCN_name)) %>% dplyr::select(c("Or_name","IUCN_name")) # Species with different_names under the IUCN Red_list

# Get the Run in parallel 
  lapply(sp_analysis$IUCN_name %>% unlist(),function(y) IUCN_red_List(x=y,export=T,exit_route = export_route))

# 1.c Download the spatial information from Gbif----
  
  points_route <- paste("./Data/Sp_info/raw_records") ; points_route %>% dir.create(recursive=TRUE,showWarnings = FALSE)
  
  for(i in 1:length(sp_analysis$IUCN_name)){
    try(Spatial_spp(sci_sp = sp_analysis$IUCN_name[i],
                    p.route = points_route,
                      start_date = 2019),
                        silent=FALSE)
    }

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