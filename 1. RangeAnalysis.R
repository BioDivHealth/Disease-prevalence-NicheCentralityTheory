rm(list=ls())
<<<<<<< HEAD
gc()
#.rs.restartR()
options(java.parameters = "-Xmx15g") # increase the memory space for jave before loading any package
options("rgdal_show_exportToProj4_warnings"="none") # Silence packages updates warnings

=======
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m")) # increase the memory for Java 
>>>>>>> b8acd646b26c63d02cee91d8c054094dbe869975
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
<<<<<<< HEAD
crs.r = "EPSG:4326"

# a. Route for the analysis----
records_r <- "./Data/Sp_records"
range_r <- "./Data/Species_Ranges"
env_r <- "./Data/Env_vars/Process"

# b. Species records----
sp_records <- records_r %>% list.files(pattern=".csv$",full.names = TRUE)

# c. Environmental variables---
env_data <- env_r %>% list.files(pattern=".tif",full.names = TRUE) %>% rast()

# d. Range species----
range_sp <- range_r %>% list.files(pattern=".shp$",full.names = TRUE) %>% st_read()

# 1. Run the environmental suitability and ENFA analysis----
for(i in 1:length(sp_records)){
  # 1.1 Get the species information----
    sp <- sp_records[i]
    species <- sp %>% basename() %>% gsub(pattern=".csv$",replacement="")  
    Sp_dist <- sp %>% read.csv() %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=crs.r)
    
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
    # Repeat this step for the points with the disease prevalence!
    dist_range <- distance_ranges(range_sp = rX ,points = Sp_dist,plot=TRUE,
                                  full=TRUE,units_d="km")
  
    gc() ; gc()
    
  # 1.4 Species Presence Probabilities or SDM----
    r.maxent <- Auto_maxent(presence_dat=Sp_dist, 
                            predictors=env_data %>% crop(rX %>% vect), 
                            rm.dp = TRUE,
                            crs.r = crs.r,
                            name.mod = species, 
                            type_bk = "Random", #[Random,BwData,BwData_inv,EnvBK]
                            world_pol = NULL, 
                            select_var = F, 
                            sp_range = rX,
                            random_features = TRUE,
                            n.m=1,
                            beta.val = c(1:15),
                            n_bk = 10000,
                            Test_n = 20,
                            # Model Selection
                            mod.select = F, n.mods = 10)
    
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
                               id_field="X")
    
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
                                  plot_sp=TRUE, # should we plot the results
                                    pts=FALSE)
                  
    # 4. Export the results ----
    exit_route <- paste("./Results/Distance_metrics/",species,sep="/")
    exit_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)
    
    # Numerical results
    write_rds(list(r.maxent[c(1:7)],ENFA.r),paste(exit_route,paste0("Dist_num",".rds"),sep="/"))
    
    # Combine and export the needed raster objects----
    rast_res<-c(r.maxent$avr.preds,p.model$trim_mod,ENFA_rast)
    writeRaster(rast_res,paste(exit_route,paste0("Distance_metrics",".tif"),sep="/"),overwrite=T)
    
    # Get the polygons for the distance calculations----
    write_rds(list(rX,p.model$pol_mod,p.model$pol_intersects),paste(exit_route,paste0("Dist_polygons",".rds"),sep="/"))
    
        # # Plot the results
    # # 4. Plot the results----
    # # Model predictions:
    # # Configure the plotting area
    # route_figs <- "./Results/Figures"; route_figs %>% dir.create(recursive=TRUE,showWarnings = FALSE)
    # png(paste(route_figs,"Polygons_AMPO.png",sep="/"),res=600,units="cm",height=14,width=18)
    # lt<-layout(matrix(c(rep(c(1,1,1,1,2,2),2),rep(1,30)),ncol=6,nrow=7,byrow = TRUE))
    # layout.show(lt)
    # colfun_p<-colorRampPalette(c("#e77c71ff","#e7d771ff","#71e7b2ff","#71b7e7ff")%>%rev())
    # # colfun_p<-colorRampPalette(c("#E40303ff","#FF8C00ff","#FFED00ff","#008026ff","#004CFFff","#732982ff")%>%rev())
    # #colfun_p<-colorRampPalette(c("black","grey50","white"))
    # 
    # # 1 mod predictions
    # plot(m.pr,axes=T,legend=F,bg=NA,mar=c(2,2,3,4),alpha=0.35,box="n",col=colfun_p(200))
    # plot(m.pr %>% mask(rX %>% st_transform(crs(m.pr))),legend=TRUE,col=colfun_p(200),add=TRUE)
    # plot(p.model$pol_intersects %>% st_geometry(),add=TRUE,border="grey36",lwd=0.5)
    # mtext(side=3,adj=0,"Model Predictions",font=2)
    # 
    # # Add the points
    # plot(points %>% st_transform(crs(m.pr)) %>% st_geometry(),col=ifelse(points$presence==1,"firebrick","grey50" %>% adjustcolor(alpha.f = 0.5)),
    #      cex=ifelse(points$presence==1,1.5,0.75),pch=19,add=TRUE)
    # plot(dist.list$links %>% st_geometry(),add=TRUE,lwd=0.5,col="grey32",lty=3)
    # 
    # # add the distance frequency
    # par(mar=c(5,2,5,5))
    # p<-hist(dist.list$data$Boundary_d,plot=F)
    # hist(dist.list$data$Boundary_d,main="Distances of points\nto AMPO polygons",
    #      xlab="Distance (km)",ylab="Frequency",col=colfun_p(p$breaks %>% length())%>%rev())
    # dev.off()
   print(paste("ALL metrics calculated for",species))
    }
=======
# Species range and global polygon
rX <- st_read("./Data/Species_Ranges" %>% list.files(pattern = "Range_analysis.gpkg",full.names = TRUE)) 
wrld_pol <- geodata::world(resolution=5,level=0,path=td) ; wrld_pol <- wrld_pol %>% st_as_sf()

# Check wrld_pol geometry
if(FALSE %in% c(wrld_pol %>% st_is_valid())){
  wrld_pol <- wrld_pol %>% st_make_valid()
  }

# wrld_pol <- wrld_pol %>% st_union() # gets problematic with some geometries, so we are going to stick with the separated polygons
# plot(wrld_pol %>% st_geometry())
# plot(rX %>% st_geometry(),col="orange",add=TRUE)

# Species points
pX <- "./Data/Sp_info/clean_records" %>% list.files(pattern = ".csv$",full.names = TRUE)

# 1. Compute the range analysis for each species----
# Analysis parameters:
  # CRS and spatial standarization
  crs_p <- "EPSG:4326"
  
  # Environmental data
  Env_variables <- "./ResultsRast/Resample_rast.tif" %>% rast()
  env_crs <- Env_variables %>% crs(describe=T) ; env_crs <- paste(env_crs$authority,env_crs$code,sep=":")

  if(env_crs != crs_p){
    env_crs <- env_crs %>% terra::project(crs=crs_p)
  }
  
  # Global Polygon data
  world_crs <- wrld_pol %>% crs(describe=T) ; world_crs <- paste(wrld_pol$authority,wrld_pol$code,sep=":")
  length(world_crs)<1
  
  
  if(length(world_crs)>0){
    if(world_crs != crs_p){ 
    wrld_pol <- wrld_pol %>% st_transform(crs=crs_p)
      }
    }
  
rm(world_crs,env_crs)  
  
# Exit routes:
analysis_r <- "./Data/Data_analysis" ; analysis_r %>% dir.create(recursive = TRUE,showWarnings = FALSE)
results_r <- "./Results/Niche_Distance" ; results_r %>% dir.create(recursive = TRUE,showWarnings = FALSE)

# 1.a Run the analysis for the different species----
  for(i in pX){
    # A. Load the point distribution data----
    points <- i %>% read.csv()
    sp <- i %>% basename() %>% gsub(pattern=".csv$",replacement = "") %>% gsub(pattern="C_",replacement="")
    
    points <- points %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=crs_p)
    
    # B. Load the spatial data and prepare the variables---- 
    range.sp <- rX %>% filter(sci_name==sp, presence %in% c(1,2,3))
    
      # b.1 Check the range crs and projection----
      range_crs <- range.sp %>% crs(describe=T) ; range_crs <- paste(range_crs$authority,range_crs$code,sep=":")
      
      if(range_crs != crs_p){
        range.sp <- range.sp %>% st_transform(crs=crs_p)
      }
      
      rm(range_crs)
      
      # b.2 Create random pseudo-absences----
      # Adjust the number of background points to be equal to 60% of the cells available in the study area
      num_bk <- Env_variables[[1]] %>% crop(range.sp %>% vect()) %>% mask(range.sp %>% vect()) %>% values()
      num_bk <- num_bk[!is.na(num_bk)] %>% length() ; num_bk <- (num_bk*0.6) %>% round(digits=0)
      
      if(num_bk>10000){
        num_bk <- 10000
      }
      
      points_bk<-backgroundPOINTS(presence=points,background_n=num_bk, # we select the default method for maxent
                                  TrainTest=c(0.7),range_samp=range.sp %>% st_bbox() %>% poly_from_ext(crs=crs(Env_variables)),
                                  weights.p="Random",buffer.dist=2,
                                  cut_area = wrld_pol) # we are going to use a global polygon to delimit the sampling area and avoid sampling points in the ocean/sea
      
    # b.2 Normalize the environmental variables----
      # d.1.2 We are going to normalize and standardize the data----
      # t_env <- list()
      # 
      # for(w in 2:nlyr(Env_variables)){ # the first column includes the cell ID so we are going to exclude it
      #   x<-Env_variables[[w]] %>% values()
      #   t_env[[w-1]]<-best_normalization(x,allow.norm=F)
      # }
      # 
      # names(t_env)<-names(env_ENFA$tab)[-1]
      # 
      # # Extract the normalize data
      # x.id<-data.frame(cells=env_ENFA$tab$cell)
      # 
      # for(w in 1:length(t_env)){
      #   x.id<-cbind(x.id,t_env[[w]]$t.values)
      #   names(x.id)[w+1]<-names(t_env)[w]
      #   
      #   if(w==1){
      #     t.method<-t_env[[w]]$method  %>% class()
      #   }else{
      #     t.method <- c(t.method,t_env[[w]]$method %>% class())
      #   }
      # }
      # 
      # # d.1.3 Standarize the data----
      # x.id[,-1]<-apply(x.id[,-1],2,function(x) (x-mean(x))/(sd(x)))
      # 
      
    # b.3 Select environmental variables (Correlation and VIF)----
    VIF.threshold <- 5 # Maximun VIF value allowed
    preserve_vars=FALSE # Do we want to preserve any specific variable for the analysis
    
        # Extract the environmental information
        env_data <- terra::extract(x=Env_variables,y=rbind(points_bk$Train,points_bk$Test) %>% vect)
        
        # d.1 Check Variables VIF values (if any transformation needs to be applied to the data do it before this step)----
        select.vars <- VIF_vars(yp=env_data,
                                vars=names(env_data)[!colnames(env_data) %in% "ID"],
                                threshold=VIF.threshold,
                                return_inf = "select")
        
        # d.2. correlation of the selected variables----
        c.int <- env_data[,colnames(env_data) %in% select.vars] %>% cor(use="complete.obs") ; c.int[upper.tri(c.int)]<-NA ; diag(c.int)<-NA
        c.var <- which(abs(c.int) > 0.7,arr.ind=TRUE) 
        
        if(nrow(c.var)==0){
          print("All correlations are lower than 0.7")
        }else{
          pairs.cor <- apply(c.var,1,function(x) data.frame(x=colnames(c.int)[x[1]],y=rownames(c.int)[x[2]])) %>% rbindlist()
              if(preserve_vars==FALSE){
                  select.vars<-names(env_data)[!colnames(env_data) %in% c(pairs.cor[,1]%>%unlist())]
          }else{
            select.vars <- names(env_data)[!colnames(env_data) %in% c(pairs.cor[,1]%>%unlist())]
            select.vars <- c(select.vars,preserve_vars)
          }
        }
    
        # d.3 Get the environmental variables for the analysis----
        env_study<-Env_variables[[names(Env_variables) %in% select.vars]] %>% crop(range.sp %>% st_buffer(dist=5) %>% vect())
    
    # C. Run the MaxEnt----
        MaxEnt_route  <- paste("./Results/Niche_Distance",sp,"MaxEnt",sep="/") ; MaxEnt_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)
        
        
      # D. ENFA analysis----
        ENFA_route  <- paste("./Results/Niche_Distance",sp,"ENFA",sep="/") ; ENFA_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)
        
        # d.1 Adapt the environmental data and Extract the position of the presence points----
        env_ENFA<-rast_to_vect(env_study)
        
        cells_points_train <- env_study %>% terra::cellFromXY(xy=st_coordinates(points))
        #cells_points_test <- env_study %>% terra::cellFromXY(xy=st_coordinates(points.test))
        
        cells_points_train <-xtabs(~cells_points_train)
        #cells_points_test <-xtabs(~cells_points_test)
        
        # Train data
        p_index_train <- rep(0,times=ncell(env_study))
        p_index_train[cells_points_train %>% names() %>% as.numeric()]<-cells_points_train
        p_index_train <- p_index_train[-c(env_ENFA$index_missin)]
        
        # Test data
        # p_index_test <- rep(0,times=ncell(env_study))
        # p_index_test[cells_points_test %>% names() %>% as.numeric()]<-cells_points_test
        # p_index_test <- p_index_test[-c(env_ENFA$index_missin)]
       
          # d.1.2 We are going to normalize and standardize the data----
          t_env <- list()

          for(w in 2:ncol(env_ENFA$tab)){ # the first column includes the cell ID so we are going to exclude it
            x<-env_ENFA$tab[,w]
            t_env[[w-1]]<-best_normalization(x,allow.norm=F)
          }

          names(t_env)<-names(env_ENFA$tab)[-1]

           # Extract the normalize data
          x.id<-data.frame(cells=env_ENFA$tab$cell)

          for(w in 1:length(t_env)){
            x.id<-cbind(x.id,t_env[[w]]$t.values)
            names(x.id)[w+1]<-names(t_env)[w]
>>>>>>> b8acd646b26c63d02cee91d8c054094dbe869975

            if(w==1){
              t.method<-t_env[[w]]$method  %>% class() %>% paste(collapse="-")
            }else{
              t.method <- c(t.method,t_env[[w]]$method %>% class() %>% paste(collapse="-"))
            }
          }

          # d.1.3 Standarize the data----
          x.id[,-1]<-apply(x.id[,-1],2,function(x) (x-mean(x))/(sd(x)))

      # d.2 Run the ENFA function----        
       enfa.1 <- ENFA_function(data=x.id[,-1],
                      presence_index = p_index_train)
                      
       print(paste("Niche centroid for",sp,"is",
                   paste(enfa.1$niche_centroid_coordinates,collapse=":"),
                   "and the marginality value for the species is of",
                   enfa.1$marginality %>% round(digits=4)))
       
       enfa.1$marginality
       enfa.1$coordinates_axis
      
          # d.2.1 Get the suitability projections----
          # Configure the data to follow the same structure as the original raster
             Suit <- cbind(x.id,enfa.1$marginality_specificity_vals,enfa.1$prediction)
             Suit.na <-data.frame(cells=env_ENFA$index_missin) ; Suit.na[,names(Suit)[-1]]<-NA
             
             Suit.1 <- rbind(Suit,Suit.na) ; Suit.1<-Suit.1[order(Suit.1$cells),]
             
             # Create the raster stack with the new data (need to check if the results make sense)
             data_ENFA<-list()
             
             for(w in 2:length(Suit.1)){ # the first column is the cell index of the original raster
               data_ENFA[[w-1]] <- rast(x=matrix(Suit.1[,w],
                                            nrow=env_ENFA$dim["rows"],
                                            ncol=env_ENFA$dim["colums"],byrow = TRUE),
                                            crs=env_ENFA$crs,
                                            extent=env_ENFA$entent)}
             
             data_ENFA<-rast(data_ENFA) ; names(data_ENFA)<-names(Suit.1)[-1]
             
             # Create a suitability map
             suitability <- data_ENFA[["Mahalanobis.Dist"]] %>% rast_01(na.rm=TRUE)
             suitability <- abs(suitability-1) # %>% plot(breaks=seq(0,1,by=0.25),col=colfun_p(4))
              
             data_ENFA$suitability <- suitability
             
             # Create the suitability classes
             suitability<-suitability %>% terra::classify(rcl=matrix(c(0,0.15,0,
                            0.15,0.25,1,0.25,0.5,2,0.5,0.75,3,0.75,1,4),byrow=TRUE,ncol=3)) #%>% plot()
                
             class_rast<-data.frame(id=0:4,suitability=c("Very-Low","Low","Intermediate","Hight","Very-hight"))
             levels(suitability)<-class_rast
             data_ENFA$Suitability_class<-suitability
             
             # Export the results
             enfa.1$data<-data_ENFA
             save(enfa.1,file = paste(ENFA_route,paste0(sp,".Rdata"),sep="/"))
             
      # E. Calculate other centrality metrics (too slow for large rasters)----    
         # # e.1 Geographical centroid----
         #   centroid_dist<-distance(env_study[[1]],range.sp %>% st_centroid() %>% vect(),unit="km") %>% mask(env_study[[1]]) 
         #     
         # # e.2 Distance to perimeter centroid----
         #   per_dist <- distance(env_study[[1]],range.sp %>% st_boundary() %>% vect(),unit="km") %>% mask(env_study[[1]]) 
         #     
         # # e.3 Distance to centroid and perimeter----
         #    x_rast<- env_study[[1]] %>% mask(range.sp %>% st_boundary() %>% vect(),updatevalue=NA,touches=TRUE)
         #    x_rast[env_study[[1]] %>% terra::cellFromXY(range.sp %>% st_centroid() %>% st_coordinates())]<-1
         #    x_rast[!is.na(x_rast)]<-1
         #    
         #    CenPer_dist <- distance(x_rast)
            
      # F. Export the results and prepare a summary report----        
        # f.1 Summary report----
            summary.1 <- data.frame(species=sp,
                                  n=nrow(points),
                                  total.modesl=length(MaxEnt_md),
                                  selected.models=length(mods),
                                  maxent_auc=lapply(MaxEnt_md,function(x) x$Test_evaluate@auc %>% round(digits = 3)) %>% unlist() %>% paste(collapse=";"),
                                  maxent_kappa=lapply(MaxEnt_md,function(x) x$Threshold_kappa) %>% unlist() %>% paste(collapse=";"),
                                  ENFA_marginality=enfa.1$marginality
                                  )
            
            summary.2 <- data.frame(vars=names(env_study),Transformation=t.method)
          
            summary.1 %>% write.csv(paste("./Results/Niche_Distance",sp,"Mods_summary.csv",sep="/"))
            summary.2 %>% write.csv(paste("./Results/Niche_Distance",sp,"Var_transformation.csv",sep="/"))
            
          # f.2 Group the result layers into a raster-stack----
            dist_r <- c(mean_maxent,sd_maxent,
                        data_ENFA$Marginality,
                        data_ENFA$Specialization1,
                        data_ENFA$Mahalanobis.Dist,
                        data_ENFA$suitability,
                        data_ENFA$Suitability_class)
             
            vect_r<-list(Range=range.sp,
                         Range_perim= range.sp %>% st_boundary(),
                         Range_centroid= range.sp %>% st_centroid(),
                         MaxEnt=pols_pred,
                         MaxEnt_cut=pols_pred_cut)
          
          # f.3 Export the results----
            terra::writeRaster(dist_r,paste("./Results/Niche_Distance",sp,"Dist_metrics.tiff",sep="/"),overwrite=TRUE)
            
            for(w in 1:length(vect_r)){
              vect_r[[w]] %>% st_write(paste("./Results/Niche_Distance",sp,paste0(names(vect_r)[w],".shp"),sep="/"),append=FALSE)
              }
          # f.4 Display the results of the anaylsis and export the results----
          # png(paste("./Results/Niche_Distance",sp,paste0("Distance_metrics_",sp,".png"),sep="/"),
          #     height = nrow(mean_maxent)*3,width = ncol(mean_maxent),units="px",res=600)
          #   
          #   lt<-layout(matrix(c(rep(1,4),2,2,4,4,3,3,4,4),byrow=TRUE,ncol=4))
          #   layout.show(lt)
          #   
          #    colfun_p<-colorRampPalette(c("#e77c71ff","#e7d771ff","#71e7b2ff","#71b9e1ff")%>%rev())
          #    colfun_p2<-colorRampPalette(c("#e77c71ff","#e7d771ff","white"))
          #    
          #    # Display the maxent results----
          #    terra::plot(mean_maxent,col=colfun_p(200),main="",axes=F) ; plot(range.sp %>% st_geometry(),add=TRUE,border="grey25") 
          #    mtext(side=3,adj=0,"MaxEnt Model Average")
          #    
          #    
          #   # ENFA results ----
          #    terra::plot(dist_r$Mahalanobis.Dist,col=colfun_p2(1000),main="Distance to niche centroid") ; plot(range.sp %>% st_geometry(),add=TRUE,border="grey25") 
          #    terra::plot(dist_r$Suitability_class,col=colfun_p2(1000) %>% rev(),main="Suitability map (1-0 Scalled)") ; plot(range.sp %>% st_geometry(),add=TRUE,border="grey25") 
          #    
          #    par(col.axis="black")
          #    plot_enfa(mar=enfa.1$marginality_specificity_vals$Marginality,
          #              spc=enfa.1$marginality_specificity_vals$Specialization1,
          #              m=enfa.1$niche_centroid_coordinates,
          #              sp_rec=enfa.1$presence_index,pts=F)
          #    axis(1); axis(2)
          #    mtext(side=3,adj=0,"ENFA results")
          #    # arrows(x0=0,y0=0,x1=enfa.1$coordinates_axis$Marginality,col="black",
          #    #        y1=enfa.1$coordinates_axis$Specialization1,length = 0)
          #    # 
          #    # text(x=enfa.1$coordinates_axis$Marginality,
          #    #      y=enfa.1$coordinates_axis$Specialization1,
          #    #      labels=row.names(enfa.1$coordinates_axis),cex=0.5,col="black")
          #    
          #    dev.off()
             
            # End of the loop
            print(paste0(sp,paste(rep("---",times=25),collapse="")))
       }
#
# End of the script
#