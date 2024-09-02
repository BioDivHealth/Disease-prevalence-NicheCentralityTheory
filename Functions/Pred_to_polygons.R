#
# Function to transform SDM predictions into spatial polygons based on probability thresholds
#
Pred_to_polygons <- function(x, # model prediction to process
                        pol.x, # Polygon defining the study area
                        t_value, # (numeric) threshold value
                        plot.r=FALSE,
                        export=NULL, # do we wanto to export the plot? [Character] route to the figures folder
                        name.mod=NULL # if export!=NULL, do we want to id the model?
  ){
  
  # Load the needed packages
  list.of.packages<-c("sf","terra","tidyverse","colorspace")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # Reclassify the rasters----
  threshold <- matrix(c(t_value,+Inf,1,-Inf,t_value,NA),ncol=3,byrow=TRUE)
  y <- terra::classify(x,rcl=threshold) # Delimitation of the more probable areas of presence for the species, according to the Kappa parameter
  
  # Transform the cut areas into polygons----
  y.p <- y %>% terra::as.polygons(aggregate=TRUE,round=TRUE,na.rm=TRUE,crs=crs(y)) %>% sf::st_as_sf()
  
  # Separate the different polygons----
  y.p <- y.p %>% sf::st_cast("POLYGON")
  y.p$Id <- 1:nrow(y.p)
  y.p$area <- st_area(y.p) ; units(y.p$area)<-"km^2"
  y.p$intersects <- y.p %>% st_intersects(pol.x) %>% as.numeric()
  
  y.cut <- y.p %>% st_intersection(y=pol.x)
  
  # Plot the spatial information----
  if(plot.r==TRUE){
  
    if(!is.null(export)){
      export %>% dir.create(recursive=TRUE,showWarnings = FALSE)
      png(paste(export,ifelse(is.null(name.mod),"Mod_to_polygons.png",paste0("Mod_to_polygons_",name.mod,".png")),sep="/"),
                res=600,units="cm",height=14,width=18)
    }
    
    
  lt<-layout(matrix(c(rep(c(1,2,1,2),each=2),rep(3,8)),ncol=4,nrow=4,byrow = TRUE))
  colfun_p<-colorRampPalette(c("#e77c71ff","#e7d771ff","#71e7b2ff","#71b7e7ff")%>%rev())
  
  # 1 mod predictions
  plot(x,axes=T,legend=F,bg=NA,mar=c(2,2,3,4),alpha=0.35,box="n",col=colfun_p(200))
  plot(x %>% mask(pol.x),legend=TRUE,col=colfun_p(200),add=TRUE)
  plot(pol.x %>% st_geometry(),add=TRUE,border="black",lwd=1)
  mtext(side=3,adj=0,"Model Predictions",font=2)
  
  # 2 mod selection
  plot(y,axes=T,col="#c9e771ff",legend=F,bg=NA,mar=c(2,2,3,4),alpha=0.35,box="n")
  plot(y %>% mask(pol.x),add=TRUE,col="#c9e771ff" %>% darken(0.15),legend=FALSE)
  plot(pol.x %>% st_geometry(),add=TRUE,border="black",lwd=1)
  legend("right",legend="Presence",pch=15,xpd=TRUE,
         pt.cex=1.2,col="#c9e771ff" %>% darken(0.15))
  mtext(side=3,adj=0,"Prediction Threshold Selection",font=2)
  
  # 3 polygons
  plot(y,axes=T,col=NA,legend=F,bg=NA,mar=c(2,2,3,4),alpha=0.35,box="n")
  plot(y.p %>% st_geometry(),add=TRUE,border=NA,col="grey90")
  plot(y.p %>% filter(intersects==1) %>% st_geometry(),add=TRUE,border="white",col="black")
  
  plot(y.cut %>% st_geometry(),add=TRUE,border=NA,col=colfun_p(nrow(y.cut)))
  plot(pol.x %>% st_geometry(),add=TRUE,border="grey50",lwd=1,lty=3)
  
  legend("topright",legend=c("Intersects","non-intersect","Inside range"),
         pch=15,xpd=TRUE,horiz=FALSE,
         pt.cex=1.2,col=c("black","grey90",colfun_p(1)))
  mtext(side=3,adj=0,"Threshold Polygons",font=2)
  
  if(!is.null(export)){
    dev.off()
    print(paste0("Plot exported to ",export))
     }
   }
  # Return the spatial features---- 
  return(list(trim_mod=y,pol_mod=y.p,pol_intersects=y.cut))
  
}

# End of the function

