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
                    "data.table","parallel","rredlist","sf")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# Load the needed functions
functions<-"./Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 1. Load the species information----
# Ana's list
data_route <- "./Data/Species_list" ; data_route %>% list.files(pattern=".xlsx$")
sp_list <- read_xlsx(data_route %>% list.files(pattern=".xlsx$",full.names = TRUE),sheet=1)

# List of captured species within the sampling localities
sp_list_site_data <- readRDS("./Data/clean_site_data.rds")
(sp_list_site_data$host_name %>% unique()) %>% length() # number of unique species

((sp_list_site_data$host_name %>% unique()) %in% sp_list$Species) %>% sum() # Species in both datasets

# Create a table with the list of species
sp_combined <- data.frame(Ana=sp_list,
                       Harry=c(unique(sp_list_site_data$host_name),
                               rep(NA,times=nrow(sp_list)-length(unique(sp_list_site_data$host_name))))
                        )

write.csv(sp_combined,"./Data/Species_list/Combined_list.csv")

# Create a full list of host/vector species
sp_list <- c(sp_list$Species,unique(sp_list_site_data$host_name))
sp_list <- unique(sp_list)

# Check the IUCN API and version
# IUCN API token
options(iucn_redlist_key="eb704359f6ea22d50235efebf7f4a2f5f843a66fa951cf6a376df71cf7268986")
RL.version<-rredlist::rl_version() ; print(paste("RedList Verion",RL.version))

# 1.b Download the IUCN species information----
export_route <- "./Data/IUCN_info" ; export_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)# output route for the IUCN information

# Check species names
sp_names <- lapply(sp_list,function(x) retrieve_syns(spp_name=x)$TaxDat) %>% rbindlist()
sp_analysis <- sp_names %>% filter(!is.na(IUCN_name)) %>% dplyr::select(c("Or_name","IUCN_name")) # Species with different_names under the IUCN Red_list

sp_names %>% write.csv(paste("./Data/Species_list","Species_analysis.csv",sep="/"),row.names = F)  

# Get the Run in parallel 
  lapply(sp_analysis$IUCN_name %>% unlist(),function(y) IUCN_red_List(x=y,export=T,exit_route = export_route))

# 1.c Download the spatial information from Gbif----
  points_route <- paste("./Data/Sp_info/raw_records") ; points_route %>% dir.create(recursive=TRUE,showWarnings = FALSE)
  
  for(i in 1:length(sp_analysis$IUCN_name)){
    try(Spatial_spp(sci_sp = sp_analysis$IUCN_name[i],
                    p.route = points_route,
                      start_date = 2000),
                        silent=FALSE)
    }

# 1.d. Check the IUCN Red List Spatial information ----
# 1.d.1 Look for the polygon
route_to_polygons <- "D:/Data/Spatial information/IUCN spatial data/All polygons"
  IUCN_pols <- lapply(sp_analysis$IUCN_name,function(w){ p <- route_to_polygons %>% 
                                                             list.files(pattern=paste0(w,".shp"),recursive = TRUE,
                                                             full.names = TRUE)
                                                              ifelse(length(p)==0,return(NA),return(p))})

  sp_analysis <- cbind(sp_analysis,polygons=unlist(IUCN_pols))  

# 1.d.2 Collect the spatial data and save them----
  polygons <- lapply(sp_analysis$polygons[!is.na(sp_analysis$polygons)],sf::st_read)
  polygons_species <- do.call("rbind",polygons)
  
  # Check the spatial data  
  if(FALSE %in% c(polygons_species %>% st_is_valid())){
    xp <- polygons_species %>% st_is_valid()
    polygons_species <- polygons_species[xp,]
    
  }
  
  polygons_species %>% st_geometry() %>% plot(col=viridis::viridis(nrow(polygons_species))%>% 
                                                adjustcolor(alpha.f = 0.15))

# 2. Clean species records ----
records_route <- data.frame(Species=points_route %>% list.files(pattern = ".csv") %>% basename() %>% gsub(pattern=".csv",replacement=""),
                            route=points_route %>% list.files(pattern = ".csv",full.names = TRUE))   
records_sp <- list()  

for(i in 1:nrow(records_route)){  
  
  range_x <- polygons_species %>% filter(BINOMIAL==records_route$Species[i])
  
  if(nrow(range_x)==0) range_x <- NULL
  
  # load the records  
  records_x <- records_route %>% filter(Species==records_route$Species[i]) %>% dplyr::select("route") %>% unlist()
  
  skip_to_next <- FALSE
  
  # Note that print(b) fails since b doesn't exist
  tryCatch(records_x <- records_x %>% read.csv(), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next){ 
    print(paste("no records for",records_route$Species[i]))
    next
    }  
  
  if(ncol(records_x)<10){
    print(paste("no records for",records_route$Species[i]))
    next
  }
  
    # Check coordinates for missing values
  obs_index <- cbind(records_x$decimalLatitude %>% is.na(),records_x$decimalLongitude %>% is.na()) %>% rowSums()
  records_x <- records_x[obs_index==0,]
  
  if(nrow(records_x)==0){
    print(paste("no records for",records_route$Species[i]))
    next
  }
  
  # Filter the points  
  records_sp[[i]] <- Prepare_points(points_sp=records_x, range_sp=range_x)
  names(records_sp)[i] <- records_route$Species[i]
  print(records_route$Species[i])
  
  }
  
# 3. Export the species spatial data----
# Summary of the data
  summary_data <- lapply(records_sp,nrow)
  summary_data <- do.call("rbind",summary_data) %>% as.data.frame()
  
  summary_data <- cbind(species=row.names(summary_data),summary_data)
  summary_data <- cbind(sp_names[match(summary_data$species,sp_names$IUCN_name),],summary_data)
  
  "./Data/Summary_data" %>% dir.create(recursive = TRUE,showWarnings = FALSE)
  write.csv(summary_data,"./Data/Summary_data/Summary_sp_records.csv")

# Export the species records
"./Data/Sp_records" %>% dir.create(recursive = TRUE,showWarnings = FALSE)

for(i in 1:length(records_sp)){
  w  <- records_sp[[i]]
  if(is.null(w)) next()
  if(nrow(w)<1) next()
  sp <- names(records_sp)[i]
  
  write.csv(w,paste("./Data/Sp_records",paste0(sp,".csv"),sep="/"))
  print(sp)
  }

# Exports the range data
st_write(polygons_species,paste("./Data/Species_Ranges",paste0("SpeciesRanges",".shp"),sep="/"), append=FALSE)

#
# ~~~~The species data is ready for the analysis~~~~
# End of the script
#