##########################.
# Complementary functions
##########################.
#
# a. Transform rast into matrix and save the paramters of the raster----
rast_to_vect<-function(x){
  #
  list.of.packages<-c("tidyr","terra","sf","data.table","dplyr")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  #
  dim_r <- c(rows=nrow(x),colums=ncol(x))
  crs_r <- crs(x)
  extend_r <- x %>% terra::ext()
  y <- x %>% terra::as.data.frame(cells=TRUE,na.rm=FALSE)
  
  data <- y[y %>% complete.cases(),]
  NA_index <- y[!y %>% complete.cases(),"cell"]
  
  return(list(dim=dim_r,
              crs=crs_r,
              entent=extend_r,
              tab=data,
              index_missin=NA_index))
}

# b. transform the extension of a raster or polygon into another polygon
poly_from_ext<-function(x,crs_p){
  
  # Load the needed packages
  list.of.packages<-c("tidyr","sf","data.table","dplyr")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # main function
  x2<-x[c("xmin", "ymin", "xmax", "ymax")]
  x_double <- as.double(x2)
  names(x_double) <- names(x2)
  class(x_double) <- "bbox"
  
  x_sf <- sf::st_as_sfc(x_double)
  sf::st_crs(x_sf) <- sf::st_crs(crs_p)
  
  return(x_sf)
}

# c. Function to evaluate the best normalization method----
# run the normalization algorithm on the environmental information
best_normalization <- function(x, # data to normalize
                               n_cores=detectCores(logical=FALSE)-1 # to run the calculations in parallel 
                               #route=NULL # route to save the results
                               ){
  # Load the needed packages
  list.of.packages<-c("tidyr","parallel","data.table","dplyr","bestNormalize")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # Configure the clusters
  if(n_cores>1){
    
    cl<-makeCluster(n_cores)
    clusterExport(cl,varlist="x",envir = .GlobalEnv)
    
    # First test
    best_normalt<-bestNormalize(x,cluster = cl,
                                allow_orderNorm = TRUE,
                                standardize=FALSE) # We are goingo to center the data later   
    stopCluster(cl)
    
    return(list(method=best_normalt$chosen_transform,t.values=best_normalt[["x.t"]]))
    
    }else{
    
      best_normalt<-bestNormalize(x,allow_orderNorm = TRUE,
                                  standardize=FALSE) # We are goingo to center the data later   
      
      return(list(method=best_normalt$chosen_transform,t.values=best_normalt[["x.t"]]))
    }
  }

