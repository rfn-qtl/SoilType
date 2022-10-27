# SoilType
A R package to interplay shovelomics in quantitative genetics and precision agriculture.

# Install
devtools::install_github("rfn-qtl/SoilType")

# Load
library(SoilType)

# get_soil(): obtain / predict  soil data from locations around the world

# Function description
?SoilType::get_soil()

# an example - single location/trials
get_soil(env.id = "RRS", lat = 30.243208 , long = -92.353191, max.lower.depth = 20)

# flow control
get_soil(env.id = c("RRS", 2), lat = 30.243208 , long = -92.353191, max.lower.depth = 20)
get_soil(env.id = "RRS", lat = 130.243208 , long = -92.353191, max.lower.depth = 20)
get_soil(env.id = "RRS", lat = 30.243208 , long = -192.353191, max.lower.depth = 20)
get_soil(env.id = "RRS", lat = 30.243208 , long = -92.353191, max.lower.depth = 0)
get_soil(env.id = "RRS", lat = 30.243208 , long = -92.353191, max.lower.depth = NA)
get_soil(env.id = NaN, lat = 30.243208 , long = -92.353191, max.lower.depth = 5)

# retrieving many samples in parallel - MET
require(foreach)
require(doParallel)
require(doMC)
# setting the number of cores that will be used
registerDoParallel(cores = detectCores())

# loading the MET
MET <- read.csv((system.file("soil_data", "LSU Rice_locations.csv", package = "SoilType")))
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
                       max.lower.depth = 20)
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
          
# get_soil_plot_level(): from a couple of samples in a trial, predict  soil characteristics at plot level

# Description
?SoilType::get_soil_plot_level()

# loading some soil samples
samples <- data.table::fread(system.file("soil_data", "soil_samples.txt", package = "SoilType"))
head(samples)

# loading the respectively shapefile
library(raster)
plot.polygons <- shapefile(system.file("soil_data", "PI_R1_DAS.shp", package = "SoilType"))
length(plot.polygons) #number of features (plots)
head(plot.polygons@data, 6)
plot.polygons[1,]@bbox

# Running the function
test3 <- get_soil_plot_level(samples = samples, plot.polygons = plot.polygons)

# and take a look at the results and the average KPI per soil trait
head(test3)
tapply(test3$Rsquared, test3$Trait, mean)
