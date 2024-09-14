########################################################
# Iterative VIF calculation (Variance Inflation Value)
########################################################
#
# Calculates the VIF for a set of variables in a dataframe and returns a table containing the 
# VIF values for the combination of variables minus the one with the highest VIF
#

VIF_vars<-function(yp, # A data.frame or matrix containing the variables from which VIF values need to be calculated
                   vars, # A list of vars we want to test
                   threshold=5, # VIF value at which to select the variables
                   return_inf = "full"# "full" or "select", select return the combination of variables that meet the threshold
                   )
                    {
  
  # 0. load the needed libraries----
  list.of.packages<-c("usdm")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # Run the calculations
   w<-yp[,colnames(yp)%in%vars]
   w[is.na(w)]<-0 
     
    VIF_dat<-list()
    Var_excluded<-list()
    yp<-data.frame(Variables=vars)
    
    # Function to calculate the VIF
    for (i in 1:ncol(w)){
      
      if (i==1){
        yp<-usdm::vif(w)
        
        # Clean up the variables (InF values are coerce to 9000 and NaN variables are exluded)
        yp[is.infinite(yp[,2]),2]<-90000
        null.var <- yp[,1][is.nan(yp[,2])]  ; yp <- yp[!is.nan(yp[,2]),]
        
        y<-yp[yp$VIF==max(yp$VIF),c(1:2)]
        
        print(y)
        
        Var_excluded[[i]]<-c(as.character(y$Variables))
        #VIF_dat[[i]]<-merge(y,yp,by="Variables",all=TRUE)
        VIF_dat[[i]]<-yp
      
      } else {
      
        if(i %in% c(ncol(w)-1,ncol(w))){
          next()
          
        }else{
          
          A<-w[,-which(names(w) %in% Var_excluded[[i-1]])]
          
          if(length(A)<2){
            Var_excluded[[i]]<-NULL
            break()
         
          }else{
          
          yp<-usdm::vif(A)
          
          # Clean up the variables (InF values are coerce to 9000 and NaN variables are exluded)
          yp[is.infinite(yp[,2]),2]<-90000
          null.var <- yp[,1][is.nan(yp[,2])]  ; yp <- yp[!is.nan(yp[,2]),]
          
          y<-yp[yp$VIF==max(yp$VIF),c(1:2)]
          
          print(y)
          
          Var_excluded[[i]]<-c(Var_excluded[[i-1]],null.var,as.character(y$Variables))
          print(y)
          
          VIF_dat[[i]]<-merge(VIF_dat[[i-1]],yp,by="Variables",all=TRUE)
          }
        }
      }
    }

    mat.vif<-VIF_dat[[length(VIF_dat)]]
    names(mat.vif)[2:ncol(mat.vif)]<-paste("vif",2:ncol(mat.vif),sep="_")
    
    if(return_inf=="full"){
      return(list(matrix=mat.vif,vars=Var_excluded))
      }
    
    if(return_inf=="select"){
        
      T.a<-apply(mat.vif[,-1],2,FUN=function(x) all(x < threshold,na.rm=TRUE))
      T.a<-names(w)[!names(w) %in% Var_excluded[[which(T.a==TRUE)[1]%>%unname()]]]
    
    return(T.a)   
    }    
 }

# End of the script

