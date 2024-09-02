# Range distances function
#
#' This function calculates the distance points to the perimeter and centroid of a given polygon or collection of polygons
#' and returns the information about the poitns, their distances to the different polygons and the location of the nearest
#' point in the perimiter of such polygon.
#' 
# THe function
distance_ranges<-function(range_sp, # [sf or sfc] Polygon containing the range of the species
                          points, # [sf or sfc] Object containing the points of the species
                          id_field=NULL, # [character] Colum in points with the observation identifiers, if this is not provided a new d_id field is created
                          plot=FALSE, # [logical] Should we plot the results of the spatial query? (default set to FALSE)
                          full=FALSE,  # [logical] If TRUE, the function returns the point distances along with all the complementary information
                          units_d="m" # [Character] In which units should the calculations be presented
                          # points_fig=NULL # if "soviet", things happen to the map
){
  # 0. Get the needed packages ready----
  list.of.packages<-c("sf","tidyverse")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  if(is.null(id_field)){
    points$d_id<-1:nrow(points)
    id_field<-"d_id"
  }
  
  # 1. Get the boundary and centroid of range_sp----  
  center_r <- st_centroid(range_sp)
  boundary_r <- st_boundary(range_sp)
  
  # Get the distance to the boundary of the range/study area
  dpoints <- st_distance(points,center_r)
  dboundary <- st_distance(points,boundary_r)
  
  # Adjust the units of the distance calculations
  units(dpoints) <-units_d
  units(dboundary) <-units_d
  
  # Create the basic data.frame
  names(units_d)<-"Distance calculated in:" ; print(units_d)
  dat_base<-data.frame(id=as.data.frame(points)[,colnames(points) %in% id_field],
                       Centroid_d=dpoints %>% as.numeric() ,Boundary_d=dboundary %>% as.numeric())
  
  if(full==TRUE|plot==TRUE){
    
    # Complementary information:
    # a. Are points intersecting with the polygons
    in_range<-st_intersects(points,range_sp,sparse=FALSE)
    
    # b. get the links and the points nearest points of species with the boundary of the range/study area
    boundary_links<-st_nearest_points(points,boundary_r) # Links
    boundary_points<-st_cast(boundary_links,"POINT") # Links points (point of origin and location in the boundary)
    
    db_2<-st_distance(boundary_points,boundary_r) #== dboundary
    units(db_2)<-units_d
    
    boundary_points <- boundary_points[!db_2 %in% dboundary]
    
    # Transform into a data.frame
    boundary_points <- st_coordinates(boundary_points) %>% as.data.frame()
    names(boundary_points)<-paste0("bp_",names(boundary_points))
    
    ad_dat<-boundary_points %>% mutate(in_range=in_range,.before=1)
  }
  
  # Plot the results
  if(plot==TRUE){
    par(mar=c(3,3,6,2)) #c(bottom, left, top, right)
    
    # Main plot
    plot(range_sp %>% st_geometry(),col="grey88") # range
    
    plot(center_r %>% st_geometry(),pch="\U2605",cex=4,col="skyblue",add=TRUE)
    plot(points %>% st_geometry(),
         col=ifelse(ad_dat$in_range==TRUE,"black"%>%adjustcolor(alpha.f=1),"black"%>%adjustcolor(alpha.f=0.35)),
         pch=19,cex=1.2,add=TRUE)  
    
    plot(boundary_links %>% st_geometry(),lwd=1.1,add=TRUE,lty=ifelse(ad_dat$in_range==TRUE,1,3),
         col=ifelse(ad_dat$in_range==TRUE,"black"%>%adjustcolor(alpha.f=1),"black"%>%adjustcolor(alpha.f=0.35)))
    
    
    # add the legend
    limits_l<-st_bbox(range_sp)
    
    y_m<-limits_l[c(2,4)] %>% max()
    x_m <-limits_l[c(1,3)] %>% mean()
    
    legend(#x=x_m,y=y_m,
      "top",
      xpd=TRUE,bty="n",
      legend=c("Range/Study\narea","Points in","Points out","Centroid"),
      pch=c("\U2B1B","\U25CF","\U25CF",ifelse(points_fig=="soviet","\U262D","\U2605")),
      col=c("grey88","black","black" %>% adjustcolor(alpha.f = 0.35),ifelse(points_fig=="soviet","tomato","skyblue")),
      lty=c(NA,1,3,NA),
      pt.cex = 3,inset=c(0,-0.1),
      adj=0,lwd=2,
      horiz=TRUE)
    
  }
  
  # Return the full dataset
  if(full==TRUE){
    
    boundary_links<-st_sf(id=dat_base[,"id"],geom=boundary_links)
    
    dat_base<-cbind(dat_base,ad_dat)
    return(list(units=units_d,data=dat_base,links=boundary_links))
    
  }else{
    return(list(units=units_d,data=dat_base))
  }
}


#
# Get the distance of a point to the nearest polygon and check if its overlaps ----
# It wraps the previous function around a set of indexes extracted using sf::st_nearest_feature
#
#

distance_p_pols<-function(points.d, # Points to calculated the distance
         polygons.d, # Target polygons
         crs.p=NULL, # Do you want to specify a CRS? Otherwise the crs of the polygons is used for all the spatial features
         plot.r=FALSE, # Should the results be plotted?
         id_field=NULL, # Does the points have an individual identifier? if not and 'ID' field is created based on its row number
         full=FALSE # Does links and other spatial features need to be exported (usefull for plotting)
){
  
  # 0. Load the needed packages
  list.of.packages<-c("sf","tidyverse","colorspace")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # 1. Check the CRS of the points and polygons
  if(!c("sf" %in% class(points.d) & "sf" %in% class(polygons.d))){
    print("ERROR: Please coerce both spatial objects to 'sf' objects!")
    stop()
  }
  
  # a. Are the CRS the same?
  if(is.null(crs.p)){
    if(crs(points.d)!=crs(polygons.d)){
      points.d <- points.d %>% st_transform(crs(polygons.d))
    }
  }else{
    # Coerce the spatial objects to the same CRS
    points.d <- points.d %>% st_transform(crs=crs.p)
    polygons.d <- polygons.d %>% st_transform(crs=crs.p)
  }
  
  # 2. Identified the closest polygon for each point
  index_pols<-st_nearest_feature(points.d,polygons.d) ; y <- polygons.d[index_pols,]
  
  if(is.null(id_field)){
    points.d$ID <-1:nrow(points.d)
    id_field<-"ID"
  }
  
  # Whit this index calculate the distance to the centroid and perimeter of such spatial feature
  d_pol<-lapply(1:nrow(points.d),function(i) distance_ranges(points = points.d[i,],
                                                             range_sp = y[i,],
                                                             plot=F,full=TRUE,
                                                             units_d="km",id_field=id_field)$data)
  d_pol <- do.call("rbind",d_pol)
  
  # Are they duplicated points in the dataset?
  if(nrow(points.d)<nrow(d_pol)){
    print("WARNING: Removing some duplicated distances!")
    d_pol <- d_pol[!duplicated(d_pol$id),]
  }
  
  
  # Prepare the export   
  if(full){
    # Get the linkages between the points and the perimeter of the polygons
    d_links<-lapply(1:nrow(points.d),function(i) distance_ranges(points = points.d[i,],
                                                                 range_sp = y[i,],
                                                                 plot=F,full=TRUE,
                                                                 units_d="km",id_field=id_field)$links)
    
    links<-do.call("rbind",d_links) 
    dat_end <- points.d %>% dplyr::select(id_field,"geom") ; dat_end<-cbind(dat_end,d_pol)
    
    return(list(data=dat_end,links=links))
    
  }else{
    dat_end <- points.d %>% dplyr::select(id_field,"geom") ; dat_end<-cbind(dat_end,d_pol)
    return(dat_end)
  }
}

# End of the functions.... so far