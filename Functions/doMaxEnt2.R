#
# _____/------\    ______  /------\\\  \\
#/             \__|      \/        / \\\
#                                      //|/\/\/\/_____---
#     Function to call Maxent       \/  __  /\\\\\
#        Analysis and testing         /\  ___  /\\\\ \\\
#                                  _      \/   |/\\ \  \ \
#        ____           / |     |\  \  \\ \\\\ 
#_______/    \_________/  |____/  \\\\\  \
#
#
run_maxent<-function(train.p, # Distribution points for the species, either a spatialpoints/dataframe, or route to a dataframe with the coordinates into columns
                     test.p,
                     bk.train,
                     bk.test,
                     dp.t=TRUE, # Should we remove duplicated points?
                     n.m=1, # number of times the model is going to be repeated to average the results (default 1)
                     random_features=TRUE, # if n.m > 1, should autofeatures be choose at random or should we try all possible combiations of them. If FALSE n.m is set automatically to reflect the number of possible feature combiantions
                     predictors, # The predictors or environmental variables
                     clip=NULL,
                     results.maxent=paste(getwd(),"spp_Maxent",sep="/"), # folder to store MAXENT results
                     mod.name,
                     jackknife=NULL,        #Logical if TRUE measures the importance of each environmental variable by training 
                     plot.maxent=TRUE,
                     return=FALSE # should models be returned to the main environment?
){
  
  # 0. load the needed libraries----
  options("rgdal_show_exportToProj4_warnings"="none") # Silence packages updates warnings!
  list.of.packages<-c("dplyr","terra","sf","dismo","rJava","generics","data.table","parallel")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # 4. Add some parameters for the model----
  # 4.1.a If the model is run multiple times, change some of the arguments so the models behave differently----
  if(n.m>1){
    models <- list()
    
    if(random_features==FALSE){
      supress_features <- list()
      
      for(i in 1:4){
        supress_features[[i]]<-combinations(x=c("nolinear","noquadratic","noproduct","nothreshold"),
                                            choose = i) %>% as.data.frame()
      }
      
      supress_features[[1]]<-supress_features[[1]] %>% t() %>% as.data.frame()
      supress_features <- supress_features %>% rbindlist(fill=TRUE)
      
      supress_features <- supress_features[rep(seq_len(nrow(supress_features)), each = 15), ]
      beta.v<-rep(c(1:15),times=15)
      
      n.m<-nrow(supress_features)
      }
    
    for(w in 1:n.m){
      
      # set.seed(185)
    if(random_features==TRUE){ 
              supress_features<-c("nolinear","noquadratic","noproduct","nothreshold")
              xFeature<-sample(supress_features,size=sample(1:4,1),replace=F) # Select the number of predictive distributions to fit (this is random)
     }else{
              xFeature<-supress_features[w,] %>% as.vector() %>% unlist() %>% unname()
              xFeature <- xFeature[!is.na(xFeature)]
     }
      
      argsMX=c(
        "autofeature=false", # Prevent the creation of automatic features or relatiohips between predictors provided to the model and the response variable
        xFeature, # Fearures to be included into the model creation
        "defaultprevalence=1.00",
        paste0("betamultiplier=",ifelse(random_features==TRUE,sample(1:15,1),beta.v[w])), # random beta multiplier set between 1-15, if random_features==FALSE all possible multipliers are applied to each combination of autofeatures (beta.v is created in line 519
        "pictures=true",
        "hinge=true",
        "writeplotdata=true", 
        paste0('threads=',ifelse((detectCores(logical=FALSE)-2)>1,detectCores(logical=FALSE)-2,1)), # Numeric. The number of processor threads to use. 
        # Matching this number to the number of cores on your
        # computer speeds up some operations, especially variable jackknifing.
        "responsecurves=false",
        # "jackknife=false", #Logical if TRUE measures the importance of each environmental variable by training 
        # with each environmental variable first omitted, then used in isolation
        "askoverwrite=false"
      )
      
      if(is.null(jackknife)){
        argsMX <- c(argsMX,"jackknife=false")
      }else{
        argsMX <- c(argsMX,"jackknife=true")
      }
      
      # 4.1.b Run the model ----
      # 4.1.b.1 Extract the variables values and configure the data.frames and vectors to run MAXENT----
      # Prepare the training data: 
      tr.dat<-rbind(predictors %>% terra::extract(train.p %>% vect()),
                   predictors %>% terra::extract(bk.train %>% vect()))
      
      tr.dat<-tr.dat[colnames(tr.dat) %in% names(predictors)]
      p.tr<-c(rep(1,times=nrow(train.p)),rep(0,times=nrow(bk.train)))
      
      # Prepare the test data:  
      p.te.dat <- predictors %>% terra::extract(test.p %>% vect())  ; p.te.dat <- p.te.dat[colnames(p.te.dat) %in% names(predictors)]
      abs.te.dat <- predictors %>% terra::extract(bk.test %>% vect()) ; abs.te.dat <- abs.te.dat[colnames(abs.te.dat) %in% names(predictors)]
      
      # Run the train models:
      Train.m <- dismo::maxent(x=tr.dat,
                               p=p.tr,
                               removeDuplicates=dp.t, 
                               args=argsMX)
      
      # # 5. Test the Train model against the test data
      threshold_seq<-seq(0.01,0.99,by=0.01)
      
      Test.m <- dismo::evaluate(p = p.te.dat,
                                a = abs.te.dat,
                                model=Train.m,
                                tr=threshold_seq,
                                type="response")
      
      kappa_threshold<-threshold(Test.m)$kappa 
      
      # 5. Alternative testing
      dat.new<-rbind(terra::extract(x=predictors,y=vect(test.p)),
                     terra::extract(x=predictors,y=vect(bk.test))
      )
      
      index.values<-dat.new %>% complete.cases()
      y.new<-c(rep(1,nrow(test.p)),rep(0,nrow(bk.test)))[index.values]
      
      # 5.1 Calculate the accuracy and specificity of our model at different probability thresholds
      y.acc<-performance_model(mod.x=Train.m,
                               dat=dat.new,
                               index_vals=y.new,
                               threshold_seq = threshold_seq)  
    if(!is.null(y.acc)){
      
      # 5.2 Plot the accuracy metrics for the different tresholds
      if(plot.maxent==TRUE){
        par(mar=c(5,5,3,3)) # c(bottom, left, top, right)
        par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
        
        col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
        color.r <- col.fun(nrow(y.acc)) 
        
        plot(y=1,x=1,xlim=c(0,1.2),ylim=c(0,1.1),
             axes=F,
             ylab="Parameter performance",
             xlab="Probability threshold",
             type="n")
        
        mtext(text = "a) Model Accuracy and performance",adj=0,line=1)
        
        for(ll in 1:nrow(y.acc)){
          if(row.names(y.acc)[ll] %in% "PCC"){
            lines(y=y.acc[ll,]/100,x=colnames(y.acc) %>% as.numeric(),
                  col=color.r[ll],type="l",pch=19)  
            
          }else{
            lines(y=y.acc[ll,],x=colnames(y.acc) %>% as.numeric(),
                  col=color.r[ll],type="l",pch=19)
          }
        }
        
        abline(v=kappa_threshold,lty=3,col="tomato3",xpd=FALSE)
        
        axis(1,seq(0,1,by=0.05),labels =seq(0,1,by=0.05) %>% round(digits = 1) ,gap.axis = 2,col="white")
        axis(2,seq(0,1,by=0.05),labels =seq(0,1,by=0.05) %>% round(digits = 1) ,gap.axis = 2,las=2,col="white")
        
        legend("topleft",legend=rownames(y.acc),lty=1,pch=19,
               col=color.r,bty="n",horiz=FALSE,ncol = 5,x.intersp = 0.5,cex=0.8)
        }
    }else{
      if(plot.maxent==TRUE){
        par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
        plot(1:10,1:10,pch=NA,axes=F,ylab=NA,xlab=NA)
        legend("center",legend="No performance data :(",bty="n")
      }
    }
      
      # 5.3 Check the variable importance and variable contribution
      if(is.null(jackknife)){
        jackniffe_var_contrybution<-jack.maxent(Train.m@results,jack=FALSE,plot.j=plot.maxent)
      }else{
        jackniffe_var_contrybution<-jack.maxent(Train.m@results,jack=TRUE,plot.j=plot.maxent)
      }
      
      # 5.4 Run some predictions on the data
      gc() ; gc()
      pred.m <- terra::predict(Train.m,predictors,cores=detectCores(logical=FALSE)-2)
      
      if(!is.null(clip)){
        pred.m %>% mask(clip %>% vect())
      }
      
      # 6. Store the results in the local folder
      mx_mod<-list(Train_model=Train.m,
                   Test=y.acc,
                   Test_evaluate=Test.m,
                   Threshold_kappa=kappa_threshold,
                   Var_contrib=jackniffe_var_contrybution,
                   # prediction=pred.m, # Gives problems when the temporal folders are removed
                   MAXENT_args=argsMX)
      
      models[[w]]<-mx_mod
      
      # 7. Plot the model
      if(plot.maxent==TRUE){  
        par(mar=c(0,0,0,0)) # c(bottom, left, top, right)
        plot(pred.m,col=col.fun(550),axes=F)
        mtext("b) Model predictions",adj=0,line=1)
      }
    }
  }else{
    
    # List of model arguments
    argsMX=c(
      "autofeature=true", # Prevent the creation of automatic features or relatiohips between predictors provided to the model and the response variable
      "defaultprevalence=1.00",
      paste0("betamultiplier=",sample(1:10,1)), # random beta multiplier set between 1-10
      "pictures=true",
      "hinge=true",
      "writeplotdata=true", 
      paste0('threads=',ifelse((detectCores(logical=FALSE)-2)>1,detectCores(logical=FALSE)-2,1)), # Numeric. The number of processor threads to use. 
      # Matching this number to the number of cores on your
      # computer speeds up some operations, especially variable jackknifing.
      "responsecurves=false",
      "askoverwrite=false")
    
    if(is.null(jackknife)){
      argsMX <- c(argsMX,"jackknife=false")
    }else{
      argsMX <- c(argsMX,"jackknife=true")
    }
    
    # 4.1.b Run the model ----
    # 4.1.b.1 Extract the variables values and configure the data.frames and vectors to run MAXENT----
    # Prepare the training data: 
    tr.dat<-rbind(predictors %>% terra::extract(train.p %>% vect()),
                  predictors %>% terra::extract(bk.train %>% vect()))
    
    tr.dat<-tr.dat[colnames(tr.dat) %in% names(predictors)]
    p.tr<-c(rep(1,times=nrow(train.p)),rep(0,times=nrow(bk.train)))
    
    # Prepare the test data:  
    p.te.dat <- predictors %>% terra::extract(test.p %>% vect())  ; p.te.dat <- p.te.dat[colnames(p.te.dat) %in% names(predictors)]
    abs.te.dat <- predictors %>% terra::extract(bk.test %>% vect()) ; abs.te.dat <- abs.te.dat[colnames(abs.te.dat) %in% names(predictors)]
    
    # Run the train models:
    Train.m <- dismo::maxent(x=tr.dat,
                             p=p.tr,
                             removeDuplicates=dp.t, 
                             args=argsMX)
    
    # # 5. Test the Train model against the test data
    threshold_seq<-seq(0.01,0.99,by=0.01)
    
    Test.m <- dismo::evaluate(p = p.te.dat,
                              a = abs.te.dat,
                              model=Train.m,
                              tr=threshold_seq,
                              type="response")
    
    kappa_threshold<-threshold(Test.m)$kappa
    
    # 5. Alternative testing
    dat.new<-rbind(terra::extract(x=predictors,y=vect(test.p)),
                   terra::extract(x=predictors,y=vect(bk.test))
    )
    
    index.values<-dat.new %>% complete.cases()
    y.new<-c(rep(1,nrow(test.p)),rep(0,nrow(bk.test)))[index.values]
    
    # 5.1 Calculate the accuracy and specificity of our model at different probability thresholds
    y.acc<-performance_model(mod.x=Train.m,
                      dat=dat.new,
                      index_vals=y.new,
                      threshold_seq = threshold_seq)
    
    if(!is.null(y.acc)){
      
      # 5.2 Plot the accuracy metrics for the different tresholds
      if(plot.maxent==TRUE){
        par(mar=c(5,5,3,3)) # c(bottom, left, top, right)
        par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
        
        col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
        color.r <- col.fun(nrow(y.acc)) 
        
        plot(y=1,x=1,xlim=c(0,1.2),ylim=c(0,1.1),
             axes=F,
             ylab="Parameter performance",
             xlab="Probability threshold",
             type="n")
        
        mtext(text = "a) Model Accuracy and performance",adj=0,line=1)
        
        for(ll in 1:nrow(y.acc)){
          if(row.names(y.acc)[ll] %in% "PCC"){
            lines(y=y.acc[ll,]/100,x=colnames(y.acc) %>% as.numeric(),
                  col=color.r[ll],type="l",pch=19)  
            
          }else{
            lines(y=y.acc[ll,],x=colnames(y.acc) %>% as.numeric(),
                  col=color.r[ll],type="l",pch=19)
          }
        }
        
        abline(v=kappa_threshold,lty=3,col="tomato3",xpd=FALSE)
        
        axis(1,seq(0,1,by=0.05),labels =seq(0,1,by=0.05) %>% round(digits = 1) ,gap.axis = 2,col="white")
        axis(2,seq(0,1,by=0.05),labels =seq(0,1,by=0.05) %>% round(digits = 1) ,gap.axis = 2,las=2,col="white")
        
        legend("topleft",legend=rownames(y.acc),lty=1,pch=19,
               col=color.r,bty="n",horiz=FALSE,ncol = 5,x.intersp = 0.5,cex=0.8)
      }
    }else{
      if(plot.maxent==TRUE){
        par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
        plot(1:10,1:10,pch=NA,axes=F,ylab=NA,xlab=NA)
        legend("center",legend="No performance data :(",bty="n")
      }
    }
  
    # 5.3 Check the variable importance and variable contribution
    if(is.null(jackknife)){
      jackniffe_var_contrybution<-jack.maxent(Train.m@results,jack=FALSE,plot.j=plot.maxent)
    }else{
      jackniffe_var_contrybution<-jack.maxent(Train.m@results,jack=TRUE,plot.j=plot.maxent)
    }
    
    # 5.4 Run some predictions on the data
    pred.m<-terra::predict(Train.m,predictors)
    
    if(!is.null(clip)){
      pred.m %>% mask(clip %>% vect())
    }
    
    # 6. Store the results in the local folder
    mx_mod<-list(Train_model=Train.m,
                 Test=y.acc,
                 Test_evaluate=Test.m,
                 Threshold_kappa=kappa_threshold,
                 Var_contrib=jackniffe_var_contrybution,
                 # prediction=pred.m, # Gives problems when the temporal folder is removed
                 MAXENT_args=argsMX)
    
    models<-mx_mod
    
    # 7. Plot the model
    if(plot.maxent==TRUE){ 
      par(mar=c(0,0,0,0)) # c(bottom, left, top, right)
      plot(pred.m,col=col.fun(550),axes=F)
      mtext("b) Model predictions",adj=0)
    }
  }
  
  if(return==FALSE){
    # Creates the Temporal file folder and the folder to allocated the results----
    if(!dir.exists(results.maxent)){dir.create(results.maxent)}
    
    results.sp<-paste(results.maxent,mod.name,sep="/") ; dir.create(results.sp,recursive = TRUE,showWarnings = FALSE)
    save(models, file=paste(results.sp,paste0("Maxent_",mod.name,".Rdata"),sep="/"))
    return(paste(mod.name,"finish!"))
  }else{
    return(models)
  }
} 
#
#
# _____/------\    ______  /------\\\  \\
#/             \__|      \/        / \\\
#                                      //|/\/\/\/_____---
#     Function to Calculate model accuracy       \/  __  /\\\\\
#       and performance at different prob         /\  ___  /\\\\ \\\
#                 Thresholds          _   \/   |/\\ \  \ \
#        ____           / |     |\  \  \\ \\\\ 
#_______/    \_________/  |____/  \\\\\  \
#
#
threshold_seq<-seq(0.05,0.95,by=0.05)


performance_model<-function(mod.x, # model to test
                            dat, # Test data
                            index_vals,
                            threshold_seq # probability treshold for the classifier
                            ){

      y.acc<-data.frame(test=NA)
      
      for(l in 1:length(threshold_seq)){
        
        # Create confusion matrix
        Conf.Mat<-table(predict=ifelse(predict(mod.x,dat,type="response") >= threshold_seq[l], 1, 0),real=index_vals)
        
        if(nrow(Conf.Mat)<2){
          next()
        }
        
        Tn<-Conf.Mat[1,1] ; Fn<-Conf.Mat[1,2]
        Fp<-Conf.Mat[2,1] ; Tp<-Conf.Mat[2,2]
        
        
        A<-list(Accuracy=(Tn+Tp)/(Tn+Tp+Fp+Fn),
                TypeI=(Fp)/(Tn+Tp+Fp+Fn),
                TypeII=(Fn)/(Tn+Tp+Fp+Fn),
                TNR=(Tn)/(Tn+Fp),
                TPR=(Tp)/(Tp+Fn),
                kappa=(2*(Tp*Tn-Fn*Fp))/((Tp+Fp)*(Fp+Tn)+(Tp+Fn)*(Fn+Tn))
        )
        
        b <- A %>% unlist() %>% as.data.frame() 
        b$test<- names(A)
        
        y.acc <- merge(y.acc,b,by="test",all=TRUE) ; colnames(y.acc)[ncol(y.acc)] <- threshold_seq[l]
        
        #rm(A,b)
      }
      
      if(length(y.acc)<2){
        print("No performance information for the Model")
        return(NULL)
      }else{
        y.acc<-y.acc[!is.na(y.acc$test),] ; rownames(y.acc)<-y.acc$test
      y.acc<-y.acc[,-1] %>% as.matrix()

return(y.acc)
        }
    }

# res <- matrix(ncol=4, nrow=length(tr))
# colnames(res) <- c('tp', 'fp', 'fn', 'tn')
# xc@t <- tr
# for (i in 1:length(tr)) {
#   res[i,1] <- length(p[p>=tr[i]])  # a  true positives
#   res[i,2] <- length(a[a>=tr[i]])  # b  false positives
#   res[i,3] <- length(p[p<tr[i]])    # c  false negatives
#   res[i,4] <- length(a[a<tr[i]])    # d  true negatives
# }
# xc@confusion = res
# a = res[,1]
# b = res[,2]
# c = res[,3]
# d = res[,4]
# # after Fielding and Bell	
# xc@np <- as.integer(np)
# xc@na <- as.integer(na)
# xc@prevalence = (a + c) / N
# xc@ODP = (b + d) / N
# xc@CCR = (a + d) / N
# xc@TPR = a / (a + c)
# xc@TNR = d / (b + d)
# xc@FPR = b / (b + d)
# xc@FNR = c/(a + c)
# xc@PPP = a/(a + b)
# xc@NPP = d/(c + d)
# xc@MCR = (b + c)/N
# xc@OR = (a*d)/(c*b)
# 
# prA = (a+d)/N
# prY = (a+b)/N * (a+c)/N
# prN = (c+d)/N * (b+d)/N
# prE = prY + prN
# xc@kappa = (prA - prE) / (1-prE)
# return(xc)
# }



# _____/------\    ______  /------\\\  \\
#/             \__|      \/        / \\\
#                                      //|/\/\/\/_____---
#     Function to plot Maxent       \/  __  /\\\\\
#       Analysis and testing         /\  ___  /\\\\ \\\
#                                  _      \/   |/\\ \  \ \
#        ____           / |     |\  \  \\ \\\\ 
#_______/    \_________/  |____/  \\\\\  \
#
#
jack.maxent<-function(x,jack=TRUE,plot.j=TRUE){
  
  # Process the raw data from maxent
  y.x <- x %>% as.data.frame()
  y.x <- y.x %>% mutate(.before = 1,Term=row.names(y.x))
  
  # A) contribution and permutation performance
  # Contribution
  contribution.x <- y.x[y.x$Term %>% grepl(pattern = ".contribution$"),]
  contribution.x <- contribution.x[order(contribution.x$V1,decreasing = FALSE),]
  contribution.x$Term <- contribution.x$Term %>% gsub(pattern=".contribution$",replacement="") 
  names(contribution.x)[2] <- "Contribution"
  
  # Permutation importance
  permutation.x <- y.x[y.x$Term %>% grepl(pattern = ".permutation.importance$"),]
  permutation.x$Term <- permutation.x$Term %>% gsub(pattern = ".permutation.importance$",replacement="") 
  names(permutation.x)[2] <- "Permutation"
  
  ycontr.x <- merge(contribution.x,permutation.x,by="Term") ; ycontr.x<-ycontr.x[order(ycontr.x$Contribution),]
  rownames(ycontr.x) <- ycontr.x$Term ; ycontr.x <- ycontr.x[,-1] %>% t()
  
  # A) Variable contribution and
  # par(mfrow=c(1,2))
  if(plot.j==TRUE){
    par(mar=c(5,7,3.5,5)) # c(bottom, left, top, right)
    par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
    
    col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
    
    x<-barplot(ycontr.x,horiz=TRUE,border=NA,beside=TRUE,
               axes=F,col=col.fun(nrow(ycontr.x)),las=2,cex.names=0.8)
    
    legend("right",legend=rownames(ycontr.x),
           col=col.fun(nrow(ycontr.x)),pch=15,bty="n")
    
    axis(1,col="white") ; abline(v=5,lty=3,lwd=1,col="tomato3",xpd=F)
    axis(2,at=apply(x,2,mean),col="white",line = NA,tick=TRUE,labels=NA,lwd=0,lwd.ticks=1)
    mtext("a) Analysis of variable Contributions and Permutation Importance",side=3,adj=0.9,las=1,xpd=TRUE)
  }
  if(jack==TRUE){
    
    # Jacknife
    Training.x <- y.x[y.x$Term %>% grepl(pattern = "Training."),]
    jack01 <- Training.x[Training.x$Term %>% grepl(pattern = "only."),] ; names(jack01)[2] <- "Only"
    jack01$Term <- jack01$Term %>% gsub(pattern="Training.gain.with.only.",replacement="")
    
    jack02 <- Training.x[Training.x$Term %>% grepl(pattern = "without."),] ; names(jack02)[2] <- "Without" 
    jack02$Term <- jack02$Term %>% gsub(pattern="Training.gain.without.",replacement="")
    
    jack03 <- merge(jack01,jack02,by="Term") ; rownames(jack03) <- jack03$Term ; jack03 <- jack03[order(jack03$Only,decreasing = FALSE),]
    jack04 <- jack03[,-1] %>% t()
    
    # B) jacknife results
    if(plot.j==TRUE){
      par(mar=c(5,7,3.5,5)) # c(bottom, left, top, right)
      par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
      
      col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
      color.r <- col.fun(2)
      
      par(mar=c(5,7,3.5,5)) # c(bottom, left, top, right)
      barplot(jack04,beside=TRUE,horiz=TRUE,
              las=2,col=color.r[c(1,2)],border=NA,axes=F,
              xlab="Regularized Training Gain")
      axis(1,col="white",line=0.5)
      mtext("b) Jackknife or regularized training gain for the MAXENT variables",side=3,adj=0.9,las=1,xpd=TRUE)
      
      legend("right",legend=row.names(jack04),
             col=color.r[c(1,2)],border = F,pch=15,inset = -0.2, xpd=TRUE,cex=0.8,bty="n")
    }
  }
  # Returnt the data
  if(jack==TRUE){
    return(list(jackknife=jack04,Var_contribution=ycontr.x))
    
  }else{
    return(list(Var_contribution=ycontr.x))  
  }
}

#
#
# _____/------\    ______  /------\\\  \\
#/             \__|      \/        / \\\
#                                      //|/\/\/\/_____---
#     Function to Calculate all the possible       \/  __  /\\\\\
#       Permutation of the model auto-features         /\  ___  /\\\\ \\\
#                                         _   \/   |/\\ \  \ \
#        ____           / |     |\  \  \\ \\\\ 
#_______/    \_________/  |____/  \\\\\  \
#
#
# modified from https://stackoverflow.com/users/4564247/pierre-l

combinations <- function(x, # vector to select
                         choose # Number of elements to consider for each permutation
                         ) {
  
  d <- do.call("expand.grid", rep(list(0:1), length(x)))
  d <-d[rowSums(d) == choose,] ; d <- ifelse(d==1,TRUE,FALSE)
  
  apply(d,1,function(w) x[w]) %>% t()
}

#
# End of The function
#