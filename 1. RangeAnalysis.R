rm(list=ls())
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m")) # increase the memory for Java 
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
                                  TrainTest=c(0.7),range_samp=range.sp,
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
        
        # c.1 shuffle the presence data----
        points <- points[sample(1:nrow(points),size=nrow(points)),]
        
        points.t <- points[sample(1:nrow(points),size=round(nrow(points)*0.7,digits=0)),]
        points.test <- points[!points$ID_p %in% points.t$ID_p,]
        
        print(paste("the sample size for training MaxEnt is of",nrow(points.t),"!"))
        dismo::maxent() # check if you can r
        
        # par(mfrow=c(2,2))
        MaxEnt_md<-run_maxent(train.p = points.t,test.p = points.test,
                              bk.train = points_bk$Train,bk.test = points_bk$Test,
                              predictors = env_study,
                              random_features = FALSE,
                              dp.t = FALSE,
                              n.m=10, # since random_features is equal to FALSE, all possible models will be run
                              results.maxent = MaxEnt_route,
                              mod.name = sp,
                              plot.maxent=FALSE,
                              jackknife = TRUE,
                              return = TRUE)
        
        # c.2 Select the 10 best performing models----
        # Get an index with performance metrics
        mod_selection<-data.frame(index=1:length(MaxEnt_md),
                                  AUC=lapply(MaxEnt_md,function(x) x$Test_evaluate@auc) %>% unlist(),
                                  TypeI=lapply(MaxEnt_md,function(x) x$Test_evaluate@FPR[x$Test_evaluate@t %in% 
                                                                                             x$Threshold_kappa]) %>% unlist())
        
        # Select the models with the highest AUC and lower TYPEI error
        mod_selection <- mod_selection[!duplicated(paste(mod_selection$AUC,mod_selection$TypeI)),] # remove the duplicated models (We can ensure no duplication by trying all possible combinaitons of autofeatures)
        mod_selection <- mod_selection %>% filter(AUC>0.7,TypeI < 0.1)
        
        if(nrow(mod_selection)>10){
          mod_selection <- mod_selection[order(mod_selection$AUC,decreasing=TRUE),]
          mods<-MaxEnt_md[c(mod_selection$index[1:10])]
          
        }else{
          mods <- MaxEnt_md[c(mod_selection$index[1:nrow(mod_selection)])]
        }      
        
        # c.3 Export the selected models----
        save(mods,file = paste(MaxEnt_route,paste0(sp,".Rdata"),sep="/"))
        
        # c.4 Run the model predictions and save them----  
        pred_maxend <- lapply(mods, function(x) terra::predict(model=x$Train_model,env_study,na.rm=TRUE)) %>% rast()# %>% mean(na.rm=TRUE)
        
        mean_maxent <-  pred_maxend %>% mean(na.rm=TRUE)
        sd_maxent <- pred_maxend %>% stdev(na.rm=TRUE)
        
        # c.5 Transform the predictions of the models into polygons using the kappa parameter----
        pols_pred <- list()
        
        for(w in 1:nlyr(pred_maxend)){
          pols_pred[[w]] <-  Pred_to_polygons(pred_maxend[[w]],pol.x=range.sp,
                           t_value = mods[[w]]$Threshold_kappa,
                           plot.r = F)$pol_mod #%>% st_geometry() %>% plot(add=ifelse(w==1,FALSE,TRUE),col=w)
          }
        
        pols_pred <- do.call("rbind",pols_pred)
        pols_pred <- pols_pred %>% st_union(by_feature = F) #%>% st_cast("POLYGON") %>% sf::st_as_sf()# merge the polygons that toucht
        pols_pred <- pols_pred[st_geometry_type(pols_pred) %in% c("MULTIPOLYGON","POLYGON")] %>% st_cast("POLYGON") %>% sf::st_as_sf()# merge the polygons that toucht
        
        pols_pred$Id <- 1:nrow(pols_pred)  
        
        if(FALSE %in% c(pols_pred %>% st_is_valid())){
          pols_pred <- pols_pred %>% st_make_valid()
          }
        
        pols_pred$area <- st_area(pols_pred) ; units(pols_pred$area)<-"km^2"
        pols_pred$intersects <- pols_pred %>% st_intersects(range.sp) %>% as.numeric()
        pols_pred$intersects <- ifelse(pols_pred$intersects==1,TRUE,FALSE)
        
        pols_pred_cut <- pols_pred %>% st_intersection(y=range.sp)
        
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

            if(w==1){
              t.method<-t_env[[w]]$method  %>% class()
            }else{
              t.method <- c(t.method,t_env[[w]]$method %>% class())
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
            dist_r <- c(mean_maxend,sd_maxent,
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
            terra::writeRaster(dist_r,paste("./Results/Niche_Distance",sp,"Dist_metrics.tiff",sep="/"))
            
            for(w in 1:length(vect_r)){
              vect_r[[w]] %>% st_write(paste("./Results/Niche_Distance",sp,paste0(names(vect_r)[w],".shp"),sep="/"))
              }
          # f.4 Display the results of the anaylsis and export the results----
          png(paste("./Results/Niche_Distance",sp,paste0("Distance_metrics_",sp,".png"),sep="/"),
              height = nrow(mean_maxent)*3,width = ncol(mean_maxent),units="px",res=600)
            
            lt<-layout(matrix(c(rep(1,4),2,2,4,4,3,3,4,4),byrow=TRUE,ncol=4))
            layout.show(lt)
            
             colfun_p<-colorRampPalette(c("#e77c71ff","#e7d771ff","#71e7b2ff","#71b9e1ff")%>%rev())
             colfun_p2<-colorRampPalette(c("#e77c71ff","#e7d771ff","white"))
             
             # Display the maxent results----
             terra::plot(mean_maxent,col=colfun_p(200),main="",axes=F) ; plot(range.sp %>% st_geometry(),add=TRUE,border="grey25") 
             mtext(side=3,adj=0,"MaxEnt Model Average")
             
             
            # ENFA results ----
             terra::plot(dist_r$Mahalanobis.Dist,col=colfun_p2(1000),main="Distance to niche centroid") ; plot(range.sp %>% st_geometry(),add=TRUE,border="grey25") 
             terra::plot(dist_r$Suitability_class,col=colfun_p2(1000) %>% rev(),main="Suitability map (1-0 Scalled)") ; plot(range.sp %>% st_geometry(),add=TRUE,border="grey25") 
             
             par(col.axis="black")
             plot_enfa(mar=enfa.1$marginality_specificity_vals$Marginality,
                       spc=enfa.1$marginality_specificity_vals$Specialization1,
                       m=enfa.1$niche_centroid_coordinates,
                       sp_rec=enfa.1$presence_index,pts=F)
             axis(1); axis(2)
             mtext(side=3,adj=0,"ENFA results")
             # arrows(x0=0,y0=0,x1=enfa.1$coordinates_axis$Marginality,col="black",
             #        y1=enfa.1$coordinates_axis$Specialization1,length = 0)
             # 
             # text(x=enfa.1$coordinates_axis$Marginality,
             #      y=enfa.1$coordinates_axis$Specialization1,
             #      labels=row.names(enfa.1$coordinates_axis),cex=0.5,col="black")
             
             dev.off()
             
            # End of the loop
            print(paste0(sp,paste(rep("---",times=25),collapse="")))
       }
#
# End of the script
#