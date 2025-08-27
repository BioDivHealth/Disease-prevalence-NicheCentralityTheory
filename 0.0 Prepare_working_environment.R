rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored
td <- tempdir()
#
#
#/\/\/\/\/\/\/\/\\/\/\/\/\/\/\/\/\/\/\/\/\/////\/\/\/\/\/\//\/\/\/>//\/\///\
##                                                                          \//\/|/      
#### Getting the custom functions from the different GitHub Repositories   
##                                                                          \//\/\/\
#/\/\/\/\/\/\/\/\\/\/\/\/\/\/\//\/\/\/\/\/\/////\/\/\/\/\/\///\/\\//\/\/\///\
#
# 0. load the needed libraries----
list.of.packages<-c("tidyverse","httr","tidyverse","remotes")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
# conflicts_prefer(dplyr::filter)
# conflicts_prefer(base::`%in%`)
# conflicts_prefer(base::`:`)

rm(list.of.packages,new.packages)

# 1. Connect to the AutoMaxent GitHub repository----
git_hub <- "https://api.github.com/repos/BioDivHealth/AutoMaxent/git/trees/main?recursive=1"
MaxRepo <- GET(git_hub) # Extract the repo information
MaxRepo

# 1.a Get the route to the functions----
file_path <- data.frame(unlist(lapply(content(MaxRepo)$tree, function(x) x$path)))
colnames(file_path) = c('Path')
head(file_path)

# Extract routes
file_path <- file_path %>%
  separate(Path,c('folder','filename'),'/') %>%
  filter(folder == 'Functions') %>%
  filter(str_detect(filename,'.R'))

# 1.b Configure the routes, download, and export scripts----
raw_route <- "https://raw.githubusercontent.com/BioDivHealth/AutoMaxent/refs/heads/main" #This is the raw route to the gitHub repository
MyRoute <- paste(getwd(),"Functions",sep="/")

for(i in 1:nrow(file_path)){
  write_lines(content(GET(paste(raw_route,file_path$folder[i],file_path$filename[i],sep="/"))),
              paste(MyRoute,file_path$filename[i],sep="/"))
}

# 2. Connect to the SDM_pipeline GitHub repository----
git_hub <- "https://api.github.com/repos/BioDivHealth/SDM_Pipeline/git/trees/main?recursive=1"
SDM_Repo <- GET(git_hub) # Extract the repo information
SDM_Repo

# 2.a Get the route to the functions----
sdm_path <- data.frame(unlist(lapply(content(SDM_Repo)$tree, function(x) x$path)))
colnames(sdm_path) = c('Path')
head(sdm_path)

# Extract routes
sdm_path <- sdm_path %>%
  separate(Path,c('folder','filename'),'/') %>%
  filter(folder == 'Functions') %>%
  filter(str_detect(filename,'.R'))

# 2.b Configure the routes, download, and export scripts----
sdm_route <- "https://raw.githubusercontent.com/BioDivHealth/SDM_Pipeline/refs/heads/main" #This is the raw route to the gitHub repository
             
for(i in 1:nrow(sdm_path)){
  write_lines(content(GET(paste(sdm_route,sdm_path$folder[i],sdm_path$filename[i],sep="/"))),
              paste(MyRoute,sdm_path$filename[i],sep="/"))
}

# 3. Configure the libraries in case of incompatibilities ----
# a Get the route to the SessionInfo folder from the AutoMaxEnt erpository (this also contains most libraries included into the SDM_pipeline repository)----
file_route <- data.frame(unlist(lapply(content(MaxRepo)$tree, function(x) x$path)))
MySession <- paste("./SessionInfo") ; MySession %>% dir.create(recursive = T,showWarnings = F)

# b Download and export the information
writeBin(content(GET(paste(raw_route,file_route[grep(file_route[,1],pattern="SessionInfo.rds"),1],sep="/")),"raw"),paste(MySession,"SessionInfo.rds",sep="/"))

# c Load the information and install the needed packages
Sinfo <- readRDS(MySession %>% list.files(".rds$",full.names = T))

# Older package versions can conflict with newer or older versions. Therefore, we are going to set up a new package library to host the
old.lib <- .libPaths()[1]
new.lib <- "./SessionInfo/LibComp" ; new.lib %>% dir.create(showWarnings = FALSE,recursive = TRUE)  

# d Install the packages
Packages <- Sinfo$otherPkgs
Pack.nmes <- Sinfo$otherPkgs %>% names()

# Check if packages are already installed
new.packages <- Pack.nmes[!(Pack.nmes %in% installed.packages(lib.loc = new.lib)[,"Package"])]
if(length(new.packages)>0){
  lapply(Packages[names(Packages) %in% new.packages], function(x) install_version(package = x$Package, version = x$Version, upgrade = "never", lib = new.lib))
  
}else{
  print("AutoMaxent dependencies already installed!")  
}

lapply(installed.packages()[,"Package"],require,character.only=TRUE)
rm(Packages,Pack.nmes)
