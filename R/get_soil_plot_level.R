#' Predict  soil characteristics at plot level
#'
#' This function, from a couple of samples in a trial, predict  soil characteristics at plot level via Random Forest.
#' @param samples data.frame The first columns must be sample_id, lat, and long.
#' The other columns are the soil traits to be extrapolated/predicted at plot level.
#' The minimum number of samples is 5.
#' @param plot.polygons Shapefile/polygons with the id, lat, and long plots information. The coordinate system must be the same for samples and the shapefile.
#' @return data.frame The output is a data.frame with the plot IDs, their lat and long coordinates, the number of samples used to make predictions, the RMSE and R-square for each prediction, the trait predicted, and the predicted value.
#' @examples
#' loading some soil samples
#' samples <- data.table::fread(system.file("soil_data", "soil_samples.txt", package = "SoilType"))
#' head(samples)

#' loading the respectively shapefile
#' library(raster)
#' plot.polygons <- shapefile(system.file("soil_data", "PI_R1_DAS.shp", package = "SoilType"))
#' length(plot.polygons) #number of features (plots)
#' head(plot.polygons@data, 6)
#' plot.polygons[1,]@bbox

#' Running the function
#' test3 <- get_soil_plot_level(samples = samples, plot.polygons = plot.polygons)

#' and take a look at the results and the average KPI per soil trait
#' head(test3)
#' tapply(test3$Rsquared, test3$Trait, mean)

#' @export
get_soil_plot_level <- function(samples,  plot.polygons){

  # extracting plot geo information
  plot.coord <- data.frame()
  for(i in 1:length(plot.polygons)){
  plot.coord <- rbind(plot.coord, data.frame(
    plot_id = i,
    lat = mean(plot.polygons[i,]@bbox[1,]) / 10000,
    long = mean(plot.polygons[i,]@bbox[2,]) /100000
  ))
  }

  # melting by trait
  samples.melted <- reshape2::melt(samples, id.vars = colnames(samples)[1:3])
  samples.melted$lat <- samples.melted$lat/10000
  samples.melted$long <- samples.melted$long/100000

  # check points
  unique.samples <- unique(samples.melted$sample)
  unique.lat <- unique(samples.melted$lat)
  unique.long <- unique(samples.melted$long)
  unique.values <- unique(samples.melted$value)
  if(anyNA(unique.samples) | length(unique.samples) <= 5) stop("The dataset must have more than 5 samples")
  if(anyNA(unique.lat) | is.numeric(unique.lat) == F | all(unique.lat < -90) | all(unique.lat > 90)) stop("The lat must numeric, and between -90 and 90)")
  if(anyNA(unique.long) | is.numeric(unique.long) == F | all(unique.long < -180) | all(unique.long > 180)) stop("The long must numeric, and between -180 and 180)")
  if(anyNA(unique.values) | is.numeric(unique.values) == F) stop("The values must numeric, and NA are not allowed)")


if (!requireNamespace("caret", quietly = TRUE)) {
  utils::install.packages("caret")
}
if (!requireNamespace("randomForest", quietly = TRUE)) {
  utils::install.packages("randomForest")
}

# index for traits
traits <- unique(samples.melted$variable)

out.plot.level <- data.frame()

for(i in 1:length(traits)){

  cat("------------------------------------------------ \n")
  cat("ATTENTION: This function uses your samples as training set \n")
  cat("------------------------------------------------  \n")
  cat('The sample size must be greater than 5 \n')
  cat("------------------------------------------------  \n")

  # subset the data
  subset.sample <- droplevels.data.frame(samples.melted[samples.melted$variable ==  traits[i],])
  training <- subset.sample
  trainControl <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 10)

  model_rf <- suppressWarnings(caret::train(x = training[, 2:3],
                                     metric = "Rsquared",
                                     y = training[, 5],
                                     maximize = TRUE,
                                     method = 'rf',
                                     trControl = trainControl))

  out.plot.level <- rbind(out.plot.level, data.frame(
    plot.coord,
    Trait = unique(subset.sample$variable),
    Predicted.value = predict(object = model_rf, plot.coord[,2:3]),
    Rsquared = round(as.numeric(model_rf$results[3]), 2),
    RMSE = round(as.numeric(model_rf$results[2]), 3),
    Samples = nrow(subset.sample)
    ))

  }

return(out.plot.level)

}
