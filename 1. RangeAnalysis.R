rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored
td<-tempdir()
dir.create(td,showWarnings = FALSE)
#
#
###//\/\/\/\/\/\/\////\/\/\/\/\/\/\///\\\\\\//\/\/\/\///////////////////////////////##-#
##                      Overlap and range distance analysis                         ##-#  
###///\/\/\/\/\/\/\////\/\/\//\/\/\/\/\/\///\\\\\\\\\..\.\\\\.\\\\\\\\\\\\\\\><\\\\\##-#
#'
#' In this case we are not going to replicate the 14 different combinations of bii/bii_diff/NCPs/WWF/WDPA that we have calculated for the local
#' scenarios. Instead we are going to calculate the values for bii and bii_diff for each country and how these values are distributed across WDPA and WWF biomes
#'
# 0. Load the packages----
list.of.packages<-c("sf","terra","tidyverse","performance","dismo","car")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0.1 Load the needed functions----
functions<-"./Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 0.2 Data for the analysis----
# Species range
rX <- st_read("./Data/Sp_info" %>% list.files(pattern = "range.gpkg",full.names = TRUE)) 
study_area <-st_read("./Data/Sp_info" %>% list.files(pattern = "study_area.gpkg",full.names = TRUE))

# Species points
pX <- st_read("./Data/Sp_info" %>% list.files(pattern = "presence.gpkg",full.names = TRUE))
pBk <- st_read("./Data/Sp_info" %>% list.files(pattern = "absence.gpkg",full.names = TRUE))

pX$presence<-1 ; pBk$presence<-0 # add a column to the dataset to mark presence and absence/Background data
points <- rbind(pX,pBk) # Combine the data

# display the data
plot(rX %>% st_geometry(),col="grey90")

p.col=ifelse(points$presence==1,"firebrick3","cyan4")
plot(points %>% st_geometry(),pch=19,col=p.col,add=TRUE)
legend("bottom",legend=c("Range","Presence","Absence"),col=c("grey90","firebrick3","cyan4"),
       pt.cex=1.2,pch=c(15,19,19),horiz=T,bty="o",xpd=TRUE)

# 1. Range distance analysis----
# Testing the function with our data
route_figs <- "./Results/Figures"; route_figs %>% dir.create(recursive=TRUE,showWarnings = FALSE)
png(paste(route_figs,"Distance_range.png",sep="/"),res=600,units="cm",height=14,width=18)

distance_ranges(range_sp = rX ,points = points,plot=TRUE,
                full=TRUE,units_d="km")

dev.off()

# 2.b Simple SDM (using glm)----
# Get the environmental data
env_route <-"./Data/Env_vars"; env_route %>% dir.create(recursive=TRUE,showWarnings = FALSE)

if(c("./Data/Env_vars/Processed" %>% list.files(recursive=FALSE) %>% length())==0){
  
  if(list.files(env_route,pattern = ".tif$",recursive = TRUE)%>%length()==0){
    # Prepare the routes to the environmental data
    paste0(env_route,"/land_cover") %>% dir.create()
    paste0(env_route,"/bioclim") %>% dir.create()
    
    # Dowload and export the data
    geodata::worldclim_global(var="bio",res=2.5,path=paste(env_route,"bioclim",sep="/"))
    geodata::landcover(var="trees",path=paste(env_route,"land_cover",sep="/"))
    geodata::landcover(var="grassland",path=paste(env_route,"land_cover",sep="/"))
    geodata::landcover(var="cropland",path=paste(env_route,"land_cover",sep="/"))
  }else{
    # Get the spatial information
    env_list<-list.files(env_route,recursive = TRUE,full.names = TRUE,pattern=".tif$")
  }
  
  # Environmental data harmonization---- 
  # We only need the data for the extension of our data
  st_crs(study_area)$units
  env_info<-resample.rast(y=env_list,paralell = TRUE,
                          ex=st_bbox(rX),
                          results.r="./Data/Env_vars/Processed",
                          mask_r=study_area %>% sf::st_buffer(dist=100000) # we wanto to preserve the data from the study area with a 100k buffer
                          )
  env <- "./Data/Env_vars/Processed" %>% list.files(full.names = TRUE) %>% rast()
  
}else{
  env <- "./Data/Env_vars/Processed" %>% list.files(full.names = TRUE) %>% rast() 
  }

# Check the variables (We need to reproject the points)
plot(env[[1]],alpha=0.3,bg="white",box="n",mar=c(2,3,5,5),legend=F) ; plot(env[[1]] %>% mask(rX %>% st_geometry() %>% st_transform(crs(env)) %>% vect()),add=TRUE,legend=T)
#plot(study_area %>% st_transform(crs(env)) %>% st_geometry(),add=TRUE,border="white",lwd=2)
plot(rX %>% st_transform(crs(env)) %>% st_geometry(),add=TRUE,border="black",lwd=1)
plot(points %>% st_transform(crs(env)) %>% st_geometry(),col=ifelse(points$presence==1,"skyblue2","grey33"),add=TRUE,pch=19)
legend("topleft",legend=c("Presence","Background"),pch=19,col=c("skyblue2","grey33"),bty="n",xpd=TRUE)
mtext(side=3,adj=0,"Distribution of observations",line=2,font=2)

# 2.c Extract the data from the rasters (terra::rast object) using the species points
Values_env <- terra::extract(x=env,y=points %>% st_transform(crs(env)) %>% vect())
points <- cbind(points,Values_env) ; class(points)

set.seed(14)
points<-points[sample(1:nrow(points)),] # randomize the order of the observations

index_t <- sample(1:nrow(points),size=c(nrow(points)*0.3) %>% round(digits=0))

train_d <- points[!c(1:nrow(points)) %in% index_t,]
test_d <- points[index_t,]

mod.full <- glm(presence~.,data=train_d[,c("presence",names(env))] %>% st_drop_geometry(),
                family=binomial) #simple glm using all environmental variables

# Very basic checks
summary(mod.full) ; performance::model_performance(mod.full)
car::Anova(mod.full) # check parameters influence

test.m<-dismo::evaluate(p=test_d[test_d$presence==1,] %>% st_drop_geometry(),
                        a=test_d[test_d$presence==0,] %>% st_drop_geometry(),
                        model=mod.full,tr=seq(0,1,by=0.015))

max(test.m@kappa) # Which is the maximun value of Kappa?
treshold_kappa<-dismo::threshold(test.m)$kappa # At which threshold we obtain the maximun value of Kappa?

# 2.d Use the model to make predictions over the data----
m.pr <- terra::predict(env,mod.full,type="response") # Raw predictions with present data

# 3. Prepare the polygons and calculate the minimun distances----
# 3.1 Extract the polygons with the best areas for the species
p.model <- Pred_to_polygons(x=m.pr,
                 pol.x= rX %>% st_transform(crs(m.pr)),
                 t_value=treshold_kappa,
                 plot.r=TRUE,
                 export="./Results/Figures",
                 name.mod="dummy")
    
# 3.2 Calculate the mininum distance of the points to the polygons----
dist.list<-distance_p_pols(points.d = points,
                polygons.d = p.model$pol_intersects,
                full = TRUE,
                id_field="ID")

# 4. Plot the results----
# Model predictions:
# Configure the plotting area
route_figs <- "./Results/Figures"; route_figs %>% dir.create(recursive=TRUE,showWarnings = FALSE)
png(paste(route_figs,"Polygons_AMPO.png",sep="/"),res=600,units="cm",height=14,width=18)
lt<-layout(matrix(c(rep(c(1,1,1,1,2,2),2),rep(1,30)),ncol=6,nrow=7,byrow = TRUE))
layout.show(lt)
colfun_p<-colorRampPalette(c("#e77c71ff","#e7d771ff","#71e7b2ff","#71b7e7ff")%>%rev())
# colfun_p<-colorRampPalette(c("#E40303ff","#FF8C00ff","#FFED00ff","#008026ff","#004CFFff","#732982ff")%>%rev())
#colfun_p<-colorRampPalette(c("black","grey50","white"))

# 1 mod predictions
plot(m.pr,axes=T,legend=F,bg=NA,mar=c(2,2,3,4),alpha=0.35,box="n",col=colfun_p(200))
plot(m.pr %>% mask(rX %>% st_transform(crs(m.pr))),legend=TRUE,col=colfun_p(200),add=TRUE)
plot(p.model$pol_intersects %>% st_geometry(),add=TRUE,border="grey36",lwd=0.5)
mtext(side=3,adj=0,"Model Predictions",font=2)

# Add the points
plot(points %>% st_transform(crs(m.pr)) %>% st_geometry(),col=ifelse(points$presence==1,"firebrick","grey50" %>% adjustcolor(alpha.f = 0.5)),
     cex=ifelse(points$presence==1,1.5,0.75),pch=19,add=TRUE)
plot(dist.list$links %>% st_geometry(),add=TRUE,lwd=0.5,col="grey32",lty=3)

# add the distance frequency
par(mar=c(5,2,5,5))
p<-hist(dist.list$data$Boundary_d,plot=F)
hist(dist.list$data$Boundary_d,main="Distances of points\nto AMPO polygons",
     xlab="Distance (km)",ylab="Frequency",col=colfun_p(p$breaks %>% length())%>%rev())
dev.off()

# End of the script