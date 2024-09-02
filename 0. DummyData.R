rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored
td<-tempdir()
dir.create(td,showWarnings = FALSE)
#
# PREPARE SOME DUMMY DATA FOR THE ANALYSIS
#
# 0. Load the packages----
list.of.packages<-c("sf","zip","tidyverse")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0. polygon_range data
range <- "./Data" %>% list.files(pattern=".gpkg$",full.names = TRUE) %>% st_read()
range <- range[sample(1:nrow(range),1),]

range %>% st_geometry() %>% plot(col="grey80")
mtext(side=3,adj=0,"Distribution of X species",col="black",font=2,line=-1)

# 1. Species points (we wanto some diversity)
rbox <- st_as_sfc(range %>% st_bbox()) # create a bounding box around the area
plot(rbox,border="tomato2",add=TRUE) ; mtext(side=3,adj=0,"Sampling area",col="tomato4",font=2,line=-3)

pX<-sf::st_sample(rbox,size=100) # create 200 random locations
pBk<-sf::st_sample(rbox %>% st_difference(range),size=100)

pX %>% plot(add=TRUE,col="purple3",pch=19) ; pBk %>% plot(add=TRUE,col="skyblue",pch=19) 
mtext(side=3,adj=0,c("Presence_points","Absence_points"),col=c("purple3","skyblue"),font=2,line=c(-4,-5))

# Export the data into gpkg format (it support larger formats than shp)
"./Data/Sp_info" %>% dir.create(recursive = TRUE,showWarnings = FALSE)

range %>% st_write("./Data/Sp_info/range.gpkg",append = FALSE)
rbox %>% st_write("./Data/Sp_info/study_area.gpkg",append = FALSE)

pX %>% st_write("./Data/Sp_info/presence.gpkg",append = FALSE)
pBk %>% st_write("./Data/Sp_info/absence.gpkg",append = FALSE)

# End of the script