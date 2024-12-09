rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored
#
#
#/\/\/\/\/\/\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/////\/\/\/\/\/\//\/\/\/>//\/\///\
##                                                                          \//\/|/      
#### Script/Function to perform a Environmental Niche Factor Analysis (ENFA)   
##                                                                          \//\/\/\
#/\/\/\/\/\/\/\/\\/\/\/\/\/\/\//\/\/\/\/\/\/////\/\/\/\/\/\///\/\\//\/\/\///\
#
#
list.of.packages<-c("tidyr","terra","sf","raster","data.table","dplyr","bestNormalize")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# Load the functions----
# 0.1 Load the needed functions----
functions<-"./Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 1. Prepare the spatial information----
# a. Load the spatial data----
env<-rast("./ResultsRast/Resample_rast.tif")
points<-read.csv("./Data/Sp_info/raw_records/Peromyscus californicus.csv")
range<-sf::st_read("./Data/Species_ranges/IUCN_range.gpkg")

plot(env[[c(1:8)]])

# b. Select the environmental variables of interest----
#env<-env[[c(1:8)]]
points <- points %>% 
          sf::st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=crs(env)) %>% 
          st_transform(crs=crs(env))

# c. Crop the environmental information to mach our study area
sty_area <- st_bbox(points) 
sty_area <- sty_area + c(-1,-1,1,1) # Add a degree to the bounding box
sty_area <- sty_area %>% poly_from_ext(crs_p=crs(points)) # Create the study area polygon

env2 <- env %>% crop(sty_area %>%vect()) # crop the raster-stack

plot(env2[[1]],axes=F,bty="n")
sty_area %>% plot(add=TRUE,lwd=2,border="tomato")
plot(points %>% sf::st_geometry(),add=TRUE,pch=19)
mtext(side=3,adj=0,text="Study area & species records")

# 2. Prepare the information for the analysis----
# a. Extract the values from the raster stack and retain the na and index information ----
env_values<-rast_to_vect(env2)

# b. Extract the position of the presence points----
cells_points <- env2 %>% terra::cellFromXY(xy=st_coordinates(points))
cells_points <-xtabs(~cells_points)

# c. Presence absence index:----
# For the ENFA analysis we need a vector of values that reflect the number of presence records
# for each "row" or cell of the raster
p_index <- rep(0,times=ncell(env2))
p_index[cells_points %>% names() %>% as.numeric()]<-cells_points
p_index <- p_index[-c(env_values$index_missin)]

# d. Check the data (it needs to be normalize and standarize for the analysis)----
t_env <- list()

for(i in 2:ncol(env_values$tab)){
  x<-env_values$tab[,i]
  t_env[[i-1]]<-best_normalization(x)
  }
names(t_env)<-names(env_values$tab)[-1]

#apply(env_values$tab[,-1],2,function(x) best_normalization(x)) # normalization of the data
# Extract the normalize data
x<-data.frame(cells=env_values$tab$cell)

for(i in 1:length(t_env)){
  x<-cbind(x,t_env[[i]]$t.values)
  names(x)[i+1]<-names(t_env)[i]
  if(i==1){
    t.method<-t_env[[i]]$method
  }else{
    t.method <- c(t.method,t_env[[i]]$method)
  }
}

# Compare the transformed and not transformed data
par(mfrow=c(5,5)) # check this before plotting
apply(env_values$tab[,!colnames(env_values$tab)%in% "cell"],2,hist,col="skyblue")
par(mfrow=c(5,5)) # check this before plotting
apply(x[,-1],2,hist,col="firebrick") # Compare transformed and raw data

# d.2 Standardize the data (x-mean)/sd (since the data has been transformed using the Ordered Quantile Normalization there is no need for this)
# x<-apply(x[,-1],2,function(p) y=(p-mean(p,na.rm=TRUE))/sd(p,na.rm = TRUE))

# 3. Run the ENFA analysis----
head(x)
enfa.1<-ENFA_function(data=x[,-1], # exclude the index
                        presence_index = p_index,
                        n_speciation_axes=1)

enfa.1$coordinates_axis

# 3.a Plot the results----
plot_enfa(mar=enfa.1$marginality_specificity_vals$Marginality,
          spc=enfa.1$marginality_specificity_vals$Specialization1,
          m=enfa.1$niche_centroid_coordinates,
          sp_rec=enfa.1$presence_index,pts=T)

arrows(x0=0,y0=0,x1=enfa.1$coordinates_axis$Marginality,
       y1=enfa.1$coordinates_axis$Specialization1,length = 0)

text(x=enfa.1$coordinates_axis$Marginality,
     y=enfa.1$coordinates_axis$Specialization1,
     labels=row.names(enfa.1$coordinates_axis),cex=0.5)

# 3.b Get the suitability projections----
# Configure the data to follow the same structure as the original raster
y <- cbind(x,enfa.1$marginality_specificity_vals,enfa.1$prediction)
y.na <-data.frame(cells=env_values$index_missin) ; y.na[,names(y)[-1]]<-NA
                  
suit.1 <- rbind(y,y.na) ; suit.1<-suit.1[order(suit.1$cells),]

# Create the raster stack with the new data (need to check if the results make sense)
data_analysis<-list()

for(i in 2:length(suit.1)){ # the first column is the cell index of the original raster
  data_analysis[[i-1]] <- rast(x=matrix(suit.1[,i],
                                nrow=env_values$dim["rows"],
                                ncol=env_values$dim["colums"],byrow = TRUE),
                                  crs=env_values$crs,
                                    extent=env_values$entent)
  }

data_analysis<-rast(data_analysis) ; names(data_analysis)<-names(suit.1)[-1]
plot(data_analysis)
