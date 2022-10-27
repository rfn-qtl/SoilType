#######################################
# SoilType R package
# Author: Roberto Fritsche-Neto
# Last update: Oct 27 2022
# email: rfn.qtl@gmail.com
#######################################

# A R package to interplay shovelomics in quantitative genetics and precision agriculture.

# Install
devtools::install_github("rfn-qtl/SoilType")

# Load
library(SoilType)

##############################################################################
# get_soil(): obtain / predict  soil data from locations around the world
##############################################################################

# Function description
?SoilType::get_soil()

#load isric dataset
soil.data <- force(soil.data)

# running the function
get_soil(env.id = "RRS", lat = 30.243208 , long = -92.353191, max.depth = 20, isric.data = soil.data)

# retrieving many samples in parallel - MET
require(foreach)
require(doParallel)
require(doMC)
# setting the number of cores that will be used
registerDoParallel(cores = detectCores())

# loading the MET
MET <- force(LSU_MET)
head(MET)
dim(MET)

# running
system.time(
  test2 <- foreach(i = 1:nrow(MET), 
                   .packages = c("caret", "stringr"), 
                   .combine = "rbind",
                   .export = c("predict", "train", "trainControl", "str_replace_all"),
                   .multicombine = TRUE, 
                   .errorhandling = "remove",
                   .verbose = TRUE    
  ) %dopar% {
    # subset the data  
    trial.MET <- droplevels.data.frame(MET[i,])
    # retrieve the data
    output <- get_soil(env.id = trial.MET$Name, 
                       lat = trial.MET$Lat, 
                       long = trial.MET$Long, 
                       max.depth = 20, 
                       isric.data = soil.data)
  }   
)

head(test2)
dim(test2)

# Soil Covariates visualization
SCov <- reshape2::dcast(test2, env ~ Trait, value.var = "Predicted.Value", mean)
head(SCov)
dim(SCov)
SCov[is.na(SCov)] <- 0
rownames(SCov) <- SCov[,1]
SCov <- scale(SCov[,2:ncol(SCov)])
SCov[is.na(SCov)] <- 0

# plots
library(heatmaply)
heatmaply(SCov, 
          fontsize_row = 6,
          fontsize_col = 6,
          file = c("SCov_heatmap.html", "SCov_heatmap.png"))

# Building a kinship (SRM - Soil Relationship Matrix)
require(EnvRtype)
SRM <- env_kernel(env.data = as.matrix(SCov), is.scaled = T, gaussian = T, sd.tol = 5)[[2]]
dim(SRM)
heatmaply(SRM, 
          fontsize_row = 6,
          fontsize_col = 6,
          file = c("SRM_heatmap.html", "SRMv_heatmap.png"))


#########################################################################################################
# get_soil_plot_level(): from a couple of samples in a trial, predict  soil characteristics at plot level
#########################################################################################################

# Description
?SoilType::get_soil_plot_level()

# loading some soil samples
data("soil.samples")

# loading the respectively shapefile
library(raster)
data("plot.polygons")
length(plot.polygons) #number of features (plots)
head(plot.polygons@data, 6)
plot.polygons[1,]@bbox

# Running the function
test3 <- get_soil_plot_level(samples = soil.samples, plot.polygons = plot.polygons)

# and take a look at the results and the average KPI per soil trait
head(test3)
tapply(test3$Rsquared, test3$Trait, mean)