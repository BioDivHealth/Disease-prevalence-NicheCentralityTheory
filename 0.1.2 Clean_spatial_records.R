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
range_spp <- paste(td,"iucn_pols",sep="/") %>% list.files(pattern = ".dbf$",recursive = TRUE,full.names = TRUE) %>% st_read()

range_spp <- range_spp %>% filter(sci_name %in% sp_list$IUCN_name)
range_spp %>% st_write(paste("./Data/Species_ranges","Range_analysis.gpkg",sep="/"),append=FALSE)

range_spp <-paste("./Data/Species_ranges","Range_analysis.gpkg",sep="/") %>% st_read()
paste(td,"iucn_pols",sep="/") %>% unlink(recursive = T)

# 2. Clean the spatial records----
sp_info <- list()
r_route <- "./Data/Sp_info/clean_records" ; r_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)

for(p in 1:length(sp_points)){
  # a Load the point data----
  p.x <- sp_points[p] %>% read.csv()
  sp <- sp_points[p] %>% basename() %>% gsub(pattern=".csv",replacement="")
  
  init.rec <- nrow(p.x)
  
  p.x <- p.x[!c(is.na(p.x$decimalLatitude)|is.na(p.x$decimalLongitude)),]
  
  if(nrow(p.x)<1){
    print(paste("No spatial informatoin for",sp))
  }
  
  # b. Check the range data----
  sp <- sp_points[p] %>% basename() %>% gsub(pattern=".csv",replacement="")
  IUCN.sp <-sp_list[sp_list$Or_name==sp,"IUCN_name"]
  
  if(length(IUCN.sp)<1){
   r.x <- NULL
   IUCN.sp  <- NA
    }else{
  r.x <- range_spp %>% filter(sci_name == IUCN.sp)
  if(nrow(r.x)<1){
    r.x <-NULL
    }
  }
  # c. Clean the Gbif points----
  p.y <- Prepare_points(points_sp=p.x,
                        range_sp=r.x,
                        xy.c=c("decimalLongitude","decimalLatitude"), # Variables describing the points coordinates (check GBIF variable naming to correctly assign thsi)
                        b.width=1, 
                        crs.r=crs_p)
  
  # d. Additional filters (this can be easily added if needed)----
  # Coordinate uncertainty
  # Incidents
  
  p.y %>% write.csv(paste(r_route,paste0("C_",sp,".csv"),sep="/"))
  
  # Display the data
  # plot(p.x$decimalLongitude,p.x$decimalLatitude,pch=19,col="grey88")
  # points(x=p.y$decimalLongitude,y=p.y$decimalLatitude,pch=19,col="firebrick")
  # plot(r.x %>% st_geometry(),add=TRUE)
  
  end.rec <- nrow(p.y)
  
  sp_info[[p]]<- data.frame(Species=sp,IUCN_name=IUCN.sp,
                            Initial.Records=init.rec,FinalRecords=end.rec,
                            IUCN_range=ifelse(is.null(r.x),FALSE,TRUE),
                            route.points=paste(r_route,paste0("C_",sp,".csv"),sep="/"))
  
 rm(p.x,sp,IUCN.sp,r.x,end.rec,init.rec,p.y)
 gc()
}

sp_info<-sp_info %>% rbindlist() ; sp_info %>% write.csv(paste("./Data/Sp_info/Records_summary.csv"))

#
# End of the script
#
unlink(td,recursive=TRUE)
