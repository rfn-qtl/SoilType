#' Retrieve and predict soil samples
#'
#' This function retrieves soil samples around the world that are near to your target location. Then, predict the soil characteristics for it.
#' If the number of samples is more than 5, then the function will use Random Forest to make predictions. Otherwise, average the information.
#' @param env.id character. Identification of the site/environment (e.g. Rayne). The env_id must have length equals 1.
#' @param lat numeric. Latitude values of the site/environment (e.g. 33.65) in WGS84. The lat must numeric, length == 1, and between -90 and 90.
#' @param long numeric. Longitude values site/environment (e.g. -90.51) in WGS84. The long must numeric, length == 1, and between -180 and 180.
#' @param max.lower.depth numeric. The lowest depth soil layer to be consider, in cm (e.g. 20). The depth must numeric, length == 1, and between 5 and 160.
#' @return data.frame The output is a data.frame with the target location and its lat and long coordinates.
#' Also, the number of samples used, the RMSE, and R-square for each prediction.
#' The chemicals characteristics are:
#' CECPH7	Cation exchange capacity - buffered at pH7 cmol(c)/kg
#' ECEC	Effective cation exchange capacity cmol(c)/kg
#' NITKJD	Total nitrogen (N) g/kg
#' ORGC	Organic carbon g/kg
#' PHAQ	pH H2O unitless
#' PHPTOT	Phosphorus (P) - total mg/kg
#' TCEQ	Calcium carbonate equivalent total g/kg
#' TOTC	Total carbon (C)
#' The physical characteristics are:
#' CLAY	Clay total g/100g
#' SAND	Sand total g/100g
#' SILT	Silt total g/100g
#' WG0006	Water retention gravimetric - 6 kPa g/100g
#' WV0006	Water retention volumetric - 6 kPa cm³/100cm³
#' @examples get_soil(env.id = "RRS", lat = 30.243208 , long = -92.353191, max.lower.depth = 20)
#' @export

get_soil <- function(env.id = NULL, lat = NULL, long = NULL, max.lower.depth = 20){

# check points
if(anyNA(env.id) | length(env.id) != 1) stop("The env_id must have length equals 1")
if(anyNA(lat) | length(lat) != 1 | is.numeric(lat) == F | all(lat < -90) | all(lat > 90)) stop("The lat must numeric, length == 1, and between -90 and 90)")
if(anyNA(long) | length(long) != 1 | is.numeric(long) == F | all(long < -180) | all(long > 180)) stop("The long must numeric, length == 1, and between -180 and 180)")
if(anyNA(max.lower.depth) | length(max.lower.depth) != 1 | is.numeric(max.lower.depth) == F | all(max.lower.depth < 5) | all(max.lower.depth > 160)) stop("The depth must numeric, length == 1, and between 5 and 160)")

#load dataset
soil.data <- force(soil.data)

# first, index samples nearby the target location
sample_id <- soil.data$profiles[round(soil.data$profiles$latitude) == round(lat) & round(soil.data$profiles$longitude) == round(long),]$profile_id
sample_id1_ <- soil.data$profiles[round(soil.data$profiles$latitude) == round(lat+1) & round(soil.data$profiles$longitude) == round(long),]$profile_id
sample_id11 <- soil.data$profiles[round(soil.data$profiles$latitude) == round(lat+1) & round(soil.data$profiles$longitude) == round(long+1),]$profile_id
sample_id_1 <- soil.data$profiles[round(soil.data$profiles$latitude) == round(lat) & round(soil.data$profiles$longitude) == round(long-1),]$profile_id
sample_id__ <- soil.data$profiles[round(soil.data$profiles$latitude) == round(lat-1) & round(soil.data$profiles$longitude) == round(long-1),]$profile_id
samples <- unique(c(sample_id, sample_id1_, sample_id11, sample_id_1, sample_id__))

# if the sample size is too small, let's find points nearby
if (length(sample_id) < 10){
nearest <- function(x, your.number) {which.min(abs(x - your.number))}
lat.nearest <- soil.data$profiles$latitude[nearest(round(soil.data$profiles$latitude), lat)]
long.nearest <- soil.data$profiles$longitude[nearest(round(soil.data$profiles$longitude), long)]
sample_id <- soil.data$profiles[round(soil.data$profiles$latitude) == round(lat.nearest) & round(soil.data$profiles$longitude) == round(long.nearest),]$profile_id
sample_id1_ <- soil.data$profiles[round(soil.data$profiles$latitude) == round(lat.nearest+1) & round(soil.data$profiles$longitude) == round(long.nearest),]$profile_id
sample_id11 <- soil.data$profiles[round(soil.data$profiles$latitude) == round(lat.nearest+1) & round(soil.data$profiles$longitude) == round(long.nearest+1),]$profile_id
sample_id_1 <- soil.data$profiles[round(soil.data$profiles$latitude) == round(lat.nearest) & round(soil.data$profiles$longitude) == round(long.nearest-1),]$profile_id
sample_id__ <- soil.data$profiles[round(soil.data$profiles$latitude) == round(lat.nearest-1) & round(soil.data$profiles$longitude) == round(long.nearest-1),]$profile_id
samples <- unique(c(sample_id, sample_id1_, sample_id11, sample_id_1, sample_id__))
  }

# picking the soils data and selecting based on max.lower.depth
profile_id <- soil.data$profiles[soil.data$profiles$profile_id %in% sample_id,]
chem <- soil.data$chemicals[soil.data$chemicals$lower_depth <= max.lower.depth & soil.data$chemicals$profile_id %in% sample_id,]
chem <- merge(profile_id, chem)
phy <- soil.data$physical[soil.data$physical$lower_depth <= max.lower.depth & soil.data$physical$profile_id %in% sample_id,]
phy <- merge(profile_id, phy)

# Modeling - Regression
# First chemical traits
output.chem <- data.frame()

for(i in 9:ncol(chem)){

soil.trait <- colnames(chem)[i]
aux <- chem
aux[,i][aux[, i] == 0] <- NA
aux <- aux[!is.na(aux[, i]),]

if(nrow(aux) >= 5){
  training <- aux
  trainControl <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 10)

model_rf <- suppressWarnings(caret::train(x = training[, 4:5],
                  metric = "Rsquared",
                  y = training[, soil.trait],
                  maximize = TRUE,
                  method = 'rf',
                  trControl = trainControl))

# Predict the labels of the test set
output.chem <- rbind(output.chem, data.frame(
  env = env.id,
  Lat = lat,
  Long = long,
  Trait = soil.trait,
  Class = "chemical",
  Predicted.Value = predict(object = model_rf, data.frame(latitude = lat, longitude = long)),
  Rsquared = round(as.numeric(model_rf$results[3]), 2),
  RMSE = round(as.numeric(model_rf$results[2]), 3),
  Samples = nrow(aux)
))
  }

if(nrow(aux) >= 1 & nrow(aux) < 10){
  output.chem <- rbind(output.chem, data.frame(
    env = env.id,
    Lat = lat,
    Long = long,
    Trait = soil.trait,
    Class = "chemical",
    Predicted.Value = mean(aux[,i]),
    Rsquared = NA,
    RMSE = NA,
    Samples = nrow(aux)
  ))
}

if(nrow(aux) < 1){
  output.chem <- rbind(output.chem, data.frame(
    env = env.id,
    Lat = lat,
    Long = long,
    Trait = soil.trait,
    Class = "chemical",
    Predicted.Value = NA,
    Rsquared = NA,
    RMSE = NA,
    Samples = NA
  ))
 }
}

# Then physical traits
output.phy <- data.frame()

for(i in 9:ncol(phy)){

  soil.trait <- colnames(phy)[i]
  aux <- phy
  aux[,i][aux[, i] == 0] <- NA
  aux <- aux[!is.na(aux[, i]),]

  if(nrow(aux) >= 5){
    training <- aux
    trainControl <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 10)

    model_rf <- suppressWarnings(caret::train(x = training[, 4:5],
                                       metric = "Rsquared",
                                       y = training[, soil.trait],
                                       maximize = TRUE,
                                       method = 'rf',
                                       trControl = trainControl))

    # Predict the labels of the test set
    output.phy <- rbind(output.phy, data.frame(
      env = env.id,
      Lat = lat,
      Long = long,
      Trait = soil.trait,
      Class = "physical",
      Predicted.Value = predict(object = model_rf, data.frame(latitude = lat, longitude = long)),
      Rsquared = round(as.numeric(model_rf$results[3]), 2),
      RMSE = round(as.numeric(model_rf$results[2]), 3),
      Samples = nrow(aux)
    ))
  }

  if(nrow(aux) >= 1 & nrow(aux) < 5){
    output.phy <- rbind(output.phy, data.frame(
      env = env.id,
      Lat = lat,
      Long = long,
      Trait = soil.trait,
      Class = "physical",
      Predicted.Value = mean(aux[,i]),
      Rsquared = NA,
      RMSE = NA,
      Samples = nrow(aux)
    ))
  }

  if(nrow(aux) < 1){
    output.phy <- rbind(output.phy, data.frame(
      env = env.id,
      Lat = lat,
      Long = long,
      Trait = soil.trait,
      Class = "physical",
      Predicted.Value = NA,
      Rsquared = NA,
      RMSE = NA,
      Samples = NA
    ))
  }

}

output <- rbind(output.chem, output.phy)

output$Trait <- stringr::str_replace_all(output$Trait,
                                         c("tceq_value_avg" =  "TCEQ",
                                           "cecph7_value_avg" =  "CECPH7",
                                           "ecec_value_avg" = "ECEC",
                                           "orgc_value_avg" = "ORGC",
                                           "phaq_value_avg" = "PHAQ",
                                           "phptot_value_avg" = "PHPTOT",
                                           "totc_value_avg" = "TOTC",
                                           "nitkjd_value_avg" = "NITKJD",
                                           "clay_value_avg" = "CLAY",
                                           "sand_value_avg" = "SAND",
                                           "silt_value_avg" = "SILT",
                                           "wg0006_value_avg" = "WG0006",
                                           "wv0006_value_avg" = "WV0006"
                                         ))
return(output)

}
