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

sp_names %>% write.csv(paste("./Data/Species_list","Species_analysis.csv",sep="/"),row.names = F)  

# Get the Run in parallel 
  lapply(sp_analysis$IUCN_name %>% unlist(),function(y) IUCN_red_List(x=y,export=T,exit_route = export_route))

# 1.c Download the spatial information from Gbif (this takes time)----
  
  points_route <- paste("./Data/Sp_info/raw_records") ; points_route %>% dir.create(recursive=TRUE,showWarnings = FALSE)
  
  for(i in 1:length(sp_analysis$IUCN_name)){
    try(Spatial_spp(#sci_sp = sp_analysis$IUCN_name[i],
                    sci_sp="Sorex minutus",
                      p.route = points_route,
                      start_date = 2000),
                        silent=FALSE)
    }

#
# This script takes time to download the spatial information
# End of the script
#  