rm(list=ls())
gc()
#.rs.restartR()
options(java.parameters = "-Xmx15g") # increase the memory space for jave before loading any package
options("rgdal_show_exportToProj4_warnings"="none") # Silence packages updates warnings

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
crs.r = "EPSG:4326"

# a. Route for the analysis----
records_r <- "./Data/Sp_records"
range_r <- "./Data/Species_Ranges"
env_r <- "./Data/Env_vars/Process"

# b. Species records----
sp_records <- records_r %>% list.files(pattern=".csv$",full.names = TRUE)

# Information about the records
sp_list <- "./Data/Summary_data/Summary_sp_records.csv" %>% read.csv()
#sp_list <- sp_list %>% filter(V1>100)

# c. Species Sampling locations----
sp_site_data <- readRDS("./Data/clean_site_data.rds")
((sp_site_data$host_name %>% unique()) %in% sp_list$species) %>% sum()

# d. Environmental variables---
env_data <- env_r %>% list.files(pattern=".tif",full.names = TRUE) %>% rast()
names(env_data)

# e. Range species----
range_sp <- range_r %>% list.files(pattern=".shp$",full.names = TRUE) %>% st_read()

# 1. Run the environmental suitability and ENFA analysis----
for(i in 1:length(sp_records)){
  
  # 1.1 Get the species information----
    sp <- sp_records[i]
    species <- sp %>% basename() %>% gsub(pattern=".csv$",replacement="")  
    Sp_dist <- sp %>% read.csv() %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=crs.r)
    
  # Sp trapping locations----
    if(species %in% (sp_site_data$host_name %>% unique())){
      "Combining sampling locations and Gbif data" %>% print()
      Sampling_locations <- sp_site_data %>% filter(host_name == species)
      
      # get unique sampling locations
      Sampling_locations <- Sampling_locations %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=crs.r)
      # Sampling_locations <- Sampling_locations[!duplicated(Sampling_locations$geometry),] # Should we remove duplicated localities?
      
      # Reduce the number of variables
      names(Sampling_locations)[1] <- "species"
      Sp_dist <- rbind(Sp_dist %>% dplyr::select(c("species","geometry")) %>% mutate(Origin="Gbif",.after="species"),
                       Sampling_locations %>% dplyr::select(c("species","geometry")) %>% mutate(Origin="Sampling",.after="species")
                       )
      
      Sp_dist <- Sp_dist %>% mutate(.before=1,id=paste(1:nrow(Sp_dist),abbreviate(species,2),sep="."))
      
    }else{
      Sp_dist <- Sp_dist %>% dplyr::select(c("species","geometry"))
      Sp_dist <- Sp_dist %>% mutate(Origin="Gbif",.after="species")
      Sp_dist <- Sp_dist %>% mutate(.before=1,id=paste(1:nrow(Sp_dist),abbreviate(species,2),sep="."))
    }
   
  # Number of observations enough
  if(nrow(Sp_dist)<20){
    print(paste("Not enough information for",species))
    next
      }  
    
  # 1.2 Prepare the spatial information----
    rX <- range_sp %>% filter(BINOMIAL == species)
  
  # If there is no range data
    if(nrow(rX)<1){
      # Create mcp from points
        p.index <- chull(Sp_dist %>% st_coordinates())
        xy.hull <- st_coordinates(Sp_dist)[c(p.index,p.index[1]),]
        
        sp.pol <- st_sf(data.frame(ID=1,geom=st_sfc(st_polygon(list(xy.hull)))),crs=crs.r)
        rX <- sp.pol %>% st_cast("POLYGON") %>% st_geometry()
        
  # # Get the study area (bounding box)
  #     rX <- sp.pol %>% st_bbox() %>% poly_from_ext(crs_p = NULL)
  #     st_crs(rX) <- crs.r ; rX <- rX %>% st_transform(crs.r)
  #     st_crs(Sp_dist) 
  #     
  #     rX <- rX %>% st_buffer(dist=25*1000)
  }
  
  # 1.3 Distance of the distribution points to the centroid and boundary of the species distribution
    dist_range <- distance_ranges(range_sp = rX ,points = Sp_dist,plot=TRUE,
                                  full=TRUE,units_d="km")
  
  # Repeat this step for the points with the disease prevalence
    if(species %in% (sp_site_data$host_name %>% unique())){
      samp_points <- sp_site_data %>% filter(host_name == species)
      samp_points <- samp_points %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=crs.r)
      
      
      dist_sampling <- distance_ranges(range_sp = rX ,points = samp_points,plot=TRUE,
                                       full=TRUE,units_d="km")
      }
    
    gc() ; gc()
    
  # 1.4 Species Presence Probabilities or SDM----
    env_model <- env_data %>% terra::crop(rX %>% vect)
    
    r.maxent <- Auto_maxent(presence_dat=Sp_dist, 
                            predictors=env_model, 
                            rm.dp = F,
                            crs.r = crs.r,
                            name.mod = species, 
                            type_bk = "Random", #[Random,BwData,BwData_inv,EnvBK]
                            world_pol = NULL, 
                            select_var = "NUMERICAL", 
                            sp_range = rX,
                            random_features = F,
                            beta.val = c(7),
                            n_bk = 10000,
                            Test_n = 20,
                            # Model Selection
                            mod.select = T, n.mods = 10)
    
    gc() ; gc()
    
    # 1.4.2 Transform the probabilities into ranges
    # Extract the polygons with the best areas for the species
    p.model <- Pred_to_polygons(x=r.maxent$avr.preds,
                                pol.x= rX,
                                t_value=r.maxent$params$TSS.mean.TEST %>% mean(na.rm=T),
                                plot.r=FALSE,
                                export=NULL,
                                name.mod="dummy")
    
    gc() ; gc()
    
    # 3.2 Calculate the mininum distance of the points to the polygons----
    dist.list<-distance_p_pols(points.d = Sp_dist,
                               polygons.d = p.model$pol_intersects,
                               full = TRUE,
                               id_field="id")
    
    gc() ; gc()
    
    # 3.3 Run the ENFA analysis----
    ind_x <- env_data %>% crop(rX %>% vect)
    
    r.dat <- rast_to_vect(ind_x)
    n.row.dat <- prod(r.dat[["dim"]])
    
    pres_index<-rep(0,times=n.row.dat)
    
    obs <- terra::extract(x=ind_x, y=Sp_dist %>% vect(), cells=T)$cell
    pres_index[obs]<-1
    
    obs_index <- pres_index[-r.dat$index_missin]
    
    ENFA.r <- ENFA_function(data = r.dat$tab[,!colnames(r.dat$tab) %in% "cell"], # Data.frame containing the environmental information with no NAs
                            presence_index = obs_index)
    
    # Transform the ENFA_results into rasters for the export
    empty_rast <- ind_x[[1]] ; empty_rast[-is.na(empty_rast)] <- NA
    
    # Get the different ENFA values
    maha <- empty_rast ; maha[r.dat$tab$cell] <- ENFA.r$prediction
    Marginality <- empty_rast ; Marginality[r.dat$tab$cell] <- ENFA.r$marginality_specificity_vals$Marginality
    Specialization <- empty_rast ; Specialization[r.dat$tab$cell] <- ENFA.r$marginality_specificity_vals$Specialization1
    
    ENFA_rast <- c(maha,Marginality,Specialization) ; names(ENFA_rast) <- c("Mahalanobis_dist","Marginality","Specificity")
    
    # 3.3.b Get the rest of the ENFA parameters
    ENFA_extra<- plot_enfa(mar=ENFA.r$marginality_specificity_vals$Marginality, # Marginality vector
                            spc=ENFA.r$marginality_specificity_vals$Specialization1, # Specialization vector
                              m=ENFA.r$niche_centroid_coordinates, # Niche centroid
                                sp_rec=obs_index, # Species records index
                                  plot_sp=F, # should we plot the results
                                    pts=FALSE)
                  
    # 4. Plot the results----
    # Output format
      exit_route <- paste("./Results/Distance_metrics",species,sep="/")
      exit_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)
    
        # Model predictions:
        # Configure the plotting area
        route_figs <- exit_route ; route_figs %>% dir.create(recursive=TRUE,showWarnings = FALSE)
        
        png(paste(route_figs,"Polygons_AMPO.png",sep="/"),res=600,units="cm",height=14,width=18)
        lt<-layout(matrix(c(rep(c(1,1,1,1,2,2),2),rep(1,30)),ncol=6,nrow=7,byrow = TRUE))
        
        #layout.show(lt)
        colfun_p<-colorRampPalette(c("#e77c71ff","#e7d771ff","#71e7b2ff","#71b7e7ff")%>%rev())
        # colfun_p<-colorRampPalette(c("#E40303ff","#FF8C00ff","#FFED00ff","#008026ff","#004CFFff","#732982ff")%>%rev())
        # colfun_p<-colorRampPalette(c("black","grey50","white"))
    
        # 1 mod predictions
        plot(r.maxent$avr.preds,axes=T,legend=F,bg=NA,mar=c(5,2,2,4),alpha=0.35,box="n",col=colfun_p(200),
             #plg = list(loc = "bottom", size=c(0.5), title = "Probability of Occurrence"))
        )
        plot(r.maxent$avr.preds %>% mask(rX %>% st_transform(crs(ENFA_rast$Marginality))),legend=TRUE,col=colfun_p(200),add=TRUE,
             plg = list(loc = "bottom", size=c(0.5), title="Probability of Occurrence",adj=0))
        
        plot(p.model$pol_intersects %>% st_geometry(),add=TRUE,border="grey36",lwd=0.5)
        mtext(side=3,adj=0,"MaxEnt Predictions",font=2)
        
        # Add the points
        plot(Sp_dist %>% st_geometry(),col="firebrick" %>% adjustcolor(alpha.f = 0.5),
             cex=1.5,pch=19,add=TRUE)
        
        plot(dist.list$links %>% st_geometry(),add=TRUE,lwd=0.5,col="grey32",lty=3)
    
        # add the distance frequency
        par(mar=c(5,2,5,5))
        p<-hist(dist.list$data$Boundary_d,plot=F)
        hist(dist.list$data$Boundary_d,main="Distances of points\nto AMPO polygons",
             xlab="Distance (km)",ylab="Frequency",col=colfun_p(p$breaks %>% length())%>%rev())
        
        dev.off()
      
      # ENFA analysis:
        png(paste(route_figs,"ENFA_metrics.png",sep="/"),res=600,units="cm",height=22,width=22)
        #lt<-layout(matrix(c(rep(c(1,1,1,1,2,2),2),rep(1,30)),ncol=6,nrow=7,byrow = TRUE))
        par(mfrow=c(2,2))
        #layout.show(lt)
        colors.r<-list(colfun_maha<-colorRampPalette(c("#e77c71ff","#e7d771ff","#71e7b2ff","#71b7e7ff")%>%rev()),
              colfun_mar<-colorRampPalette(c("#E40303ff","#FF8C00ff","#FFED00ff","#732982ff")%>%rev()),
                colfun_spe<-colorRampPalette(c("black","grey50","white"))
                )
        
        # ENFA results
        for(w in 1:nlyr(ENFA_rast)){
        # Configure the color
          col.rast<-do.call(colors.r[[w]],list(200))
        
        # Configure the layer  
        if(names(ENFA_rast)[w]=="Mahalanobis_dist"){
          x <- ENFA_rast[[w]] %>% log1p()
            }else{
          x <- ENFA_rast[[w]]
          }
        
        plot(x,axes=T,legend=F,bg=NA,mar=c(3,2,2,4),alpha=0.35,box="n",col=col.rast,
             #plg = list(loc = "bottom", size=c(0.5), title = "Probability of Occurrence"))
        )
        plot(x %>% mask(rX),axes=T,legend=T,bg=NA,mar=c(5,2,2,4),box="n",col=col.rast,add=TRUE,
             plg = list(loc = "bottom", size=c(0.5), title=names(ENFA_rast)[w],adj=0))
        
        plot(rX %>% st_geometry(),add=TRUE,border="grey36",lwd=0.5)
        mtext(side=3,adj=0,ifelse(w==1,"Distance to niche centroid\nlog(Mahalanobis+1)",names(ENFA_rast)[w]),font=2,cex=0.8)
        
        # Add the points
        plot(Sp_dist %>% st_geometry(),col="firebrick" %>% adjustcolor(alpha.f = 0.5),
             cex=1.5,pch=19,add=TRUE)
        }
        
        # Plot the Marginality and Specificity Dimensions
        plot_enfa(mar=ENFA.r$marginality_specificity_vals$Marginality, # Marginality vector
                  spc=ENFA.r$marginality_specificity_vals$Specialization1, # Specialization vector
                  m=ENFA.r$niche_centroid_coordinates, # Niche centroid
                  sp_rec=obs_index, # Species records index
                  plot_sp=T, # should we plot the results
                  pts=FALSE)
       
      dev.off()
  
   # 5. Export the results ----
        # Generate the final table----
        dist_metrics <- cbind(Sp_dist,terra::extract(c(r.maxent$avr.preds,ENFA_rast),Sp_dist %>% terra::vect()))
        dist_metrics <- merge(dist_metrics %>% st_drop_geometry(),dist.list$data %>% st_drop_geometry(),by="id")
        
        dist_metrics %>% write.csv(paste(exit_route,"Distance_metrics.csv",sep="/"))# Export the dataframe
        
        # Numerical results
        write_rds(list(MaxEnt=r.maxent[c(1:7)],ENFA=ENFA.r),paste(exit_route,paste0("Dist_num",".rds"),sep="/"))
        
        # Combine and export the needed raster objects----
        rast_res<-c(r.maxent$avr.preds,p.model$trim_mod,ENFA_rast)
        writeRaster(rast_res,paste(exit_route,paste0("Distance_metrics",".tif"),sep="/"),overwrite=T)
        
        # Get the polygons for the distance calculations----
        write_rds(list(rX,p.model$pol_mod,p.model$pol_intersects),paste(exit_route,paste0("Dist_polygons",".rds"),sep="/"))
     
    if(exists("Sampling_locations")){
      rm("Sampling_locations")
        } # reset some important parameters
   
    rm(p,r.maxent,ENFA.r,ENFA_extra,ENFA_rast,rast_res,exit_route,r.maxent,w,Sp_dist,rX)
    print(paste("ALL metrics calculated for",species))
    gc() ; gc()
}

#
# End of the script
#