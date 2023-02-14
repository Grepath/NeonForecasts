library(tidyverse)
library(neon4cast)
library(lubridate)
library(rMR)
library(glue)
library(xgboost)


# setwd("/Users/gregharrison/Documents/GitHub/NeonForecasts")

source("ignore_sigpipe.R")

forecast_date <- Sys.Date()
forecast_doy = as.numeric(format(forecast_date, '%j'))
noaa_date <- Sys.Date() - days(1)  #Need to use yesterday's NOAA forecast because today's is not available yet

# Step 0: Define a unique name which will identify your model in the leaderboard and connect it to team members info, etc
model_id <- "xgboost_temp_oxygen_chla_sequential"

# Step 1: Download latest target data and site description data
target <- readr::read_csv(paste0("https://data.ecoforecast.org/neon4cast-targets/",
                                 "aquatics/aquatics-targets.csv.gz"), guess_max = 1e6)
site_data <- readr::read_csv(paste0("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/",
                                    "main/NEON_Field_Site_Metadata_20220412.csv")) |> 
  dplyr::filter(aquatics == 1)


# Step 2: Get meterological predictions as drivers
df_past <- neon4cast::noaa_stage3()

## Helper function: for each site, average over predicted 0h horizon ensembles to get 'historic values'
noaa_mean_historical <- function(df_past, site, var) {
  df_past |>
    dplyr::filter(site_id == site,
                  variable == var) |>
    dplyr::rename(ensemble = parameter) |>
    dplyr::select(datetime, prediction, ensemble) |>
    dplyr::mutate(date = as_date(datetime)) |>
    dplyr::group_by(date) |>
    dplyr::summarize(air_temperature = mean(prediction, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::rename(datetime = date) |>
    dplyr::mutate(air_temperature = air_temperature - 273.15) |>
    dplyr::collect()
}

## Helper fn: get daily average temperature from each ensemble in future
noaa_mean_forecast <- function(site, var, reference_date) {
  endpoint = "data.ecoforecast.org"
  bucket <- glue::glue("neon4cast-drivers/noaa/gefs-v12/stage1/0/{reference_date}")
  s3 <- arrow::s3_bucket(bucket, endpoint_override = endpoint, anonymous = TRUE)
  
  # stage1 air temp is Celsius
  arrow::open_dataset(s3) |>
    dplyr::filter(site_id == site,
                  datetime >= lubridate::as_datetime(forecast_date),
                  variable == var) |>
    dplyr::select(datetime, prediction, parameter) |>
    dplyr::mutate(datetime = as_date(datetime)) |>
    dplyr::group_by(datetime, parameter) |>
    dplyr::summarize(air_temperature = mean(prediction), .groups = "drop") |>
    dplyr::select(datetime, air_temperature, parameter) |>
    dplyr::rename(ensemble = parameter) |>
    dplyr::collect()
  
}

# Step 2.5: We'll skip any site that doesn't have both temperature and oxygen
sites <- target |> na.omit() |> distinct(site_id, variable) |> 
  filter(variable %in% c("oxygen", "temperature", "chla")) |>
  count(site_id) |> filter(n==3) |> pull(site_id)

#Step 3.0: Define the forecasts model for a site
forecast_site <- function(site) {
  message(paste0("Running site: ", site))
  
  # Get site information for elevation
  #site_info <- site_data |> dplyr::filter(field_site_id == site)
  
  # historical temperatures
  noaa_past_mean <- noaa_mean_historical(df_past, site, "air_temperature")
  
  # Merge in past NOAA data into the targets file, matching by date.
  site_target <- target |>
    dplyr::select(datetime, site_id, variable, observation) |>
    dplyr::filter(variable %in% c("temperature", "oxygen", "chla"), 
                  site_id == site) |>
    tidyr::pivot_wider(names_from = "variable", values_from = "observation") |>
    dplyr::left_join(noaa_past_mean, by = c("datetime"))
  
  rm(noaa_past_mean) # save RAM 
  
  site_target[['dayofyear']] <- as.numeric( format(site_target[['datetime']], '%j'))
  site_target = na.omit(site_target)
  
  # historical <- site_target %>% filter(dayofyear == forecast_doy)
  startCheck = (forecast_doy - 7) %% 365
  
  historical <- site_target %>% 
    filter(dayofyear <= forecast_doy) %>% 
    filter(dayofyear > startCheck)
  
  forecast = data.frame()
  
  if(nrow(historical) < 1){
    noaa_future <- noaa_mean_forecast(site, "TMP", noaa_date)
    noaa_future[['dayofyear']] <- as.numeric( format(noaa_future[['datetime']], '%j'))
    future_input <- as.matrix(noaa_future[,c("air_temperature", "dayofyear")])
    
    emptyForecast = rep(NA,nrow(future_input))
    
    temperature <- 
      noaa_future |> 
      mutate(site_id = site,
             prediction = emptyForecast,
             variable = "temperature")
    
    oxygen <- 
      noaa_future |> 
      mutate(site_id = site,
             prediction = emptyForecast,
             variable = "oxygen")
    
    chla <- 
      noaa_future |> 
      mutate(site_id = site,
             prediction = emptyForecast,
             variable = "chla")
    
    forecast <- dplyr::bind_rows(temperature, oxygen, chla)
    
    
    message("Dropped Site")
  
  } else {
    sample <- sample(c(TRUE, FALSE), nrow(site_target), replace=TRUE, prob=c(0.8,0.2))
    train  <- site_target[sample, ]
    test   <- site_target[!sample, ]
    
    # Generate our labels as the current water temperature
    train.label = train$temperature
    test.label = test$temperature
    
    # Convert the input data to a matrix for xgboost
    train.data = as.matrix(train[, c("air_temperature", "dayofyear")])
    test.data = as.matrix(test[, c("air_temperature", "dayofyear")])
    
    # Generate Training Input for XGBoost
    dtrain<-xgb.DMatrix(data = train.data, label = train.label)
    # Train our model
    bst <- xgboost(data = dtrain, max.depth = 10, eta = 0.3, nthread = 2, nrounds = 15, verbose = 0)
    # Product Predictions for our Testing Dataset
    test_pred <- predict(bst, test.data)
    # Calculate Root Mean Squared Error
    rmse <- sqrt(mean((test_pred-test.label)^2))
    message('RMSE: ', rmse)
    
    
    standardDev = sd(test_pred-test.label)
    
    #  Get 30-day predicted temperature ensemble at the site
    noaa_future <- noaa_mean_forecast(site, "TMP", noaa_date)
    
    
    noaa_future[['dayofyear']] <- as.numeric( format(noaa_future[['datetime']], '%j'))
    future_input <- as.matrix(noaa_future[,c("air_temperature", "dayofyear")])
    
    temperaturePrediction = predict(bst, future_input)+rnorm(nrow(future_input))*standardDev
    
    temperature <- 
      noaa_future |> 
      mutate(site_id = site,
             prediction = temperaturePrediction,
             variable = "temperature")
    
    
    # use the same weather forecast parameters to train a model for Oxygen 
    
    # oxySample <- sample(c(TRUE, FALSE), nrow(site_target), replace=TRUE, prob=c(0.8,0.2))
    # oxyTrain  <- site_target[sample, ]
    # oxyTest   <- site_target[!sample, ]
    
    # Generate our labels as the current water temperature
    oxytrain.label = train$oxygen
    oxytest.label = test$oxygen
    
    # Convert the input data to a matrix for xgboost
    # oxytrain.data = as.matrix(train[, c("air_temperature", "dayofyear")])
    # oxytest.data = as.matrix(test[, c("air_temperature", "dayofyear")])
    # 
    # # Generate Training Input for XGBoost
    # oxydtrain<-xgb.DMatrix(data = oxytrain.data, label = oxytrain.label)
    # # Train our model
    # oxybst <- xgboost(data = oxydtrain, max.depth = 10, eta = 0.3, nthread = 2, nrounds = 15, verbose = 0)
    # # Product Predictions for our Testing Dataset
    # oxytest_pred <- predict(oxybst, oxytest.data)
    # # Calculate Root Mean Squared Error
    # rmse <- sqrt(mean((oxytest_pred-oxytest.label)^2))
    # standardDev = sd(oxytest_pred-oxytest.label)
    # message('RMSE: ', rmse)
    # 
    # noaa_future[['temperature']] <- temperaturePrediction
    # future_input <- as.matrix(noaa_future[,c("air_temperature", "dayofyear")])
    
    # Convert the input data to a matrix for xgboost
    oxytrain.data = as.matrix(train[, c("air_temperature", "dayofyear", "temperature")])
    oxytest.data = as.matrix(test[, c("air_temperature", "dayofyear", "temperature")])
    
    # Generate Training Input for XGBoost
    oxydtrain<-xgb.DMatrix(data = oxytrain.data, label = oxytrain.label)
    # Train our model
    oxybst <- xgboost(data = oxydtrain, max.depth = 10, eta = 0.3, nthread = 2, nrounds = 15, verbose = 0)
    # Product Predictions for our Testing Dataset
    oxytest_pred <- predict(oxybst, oxytest.data)
    # Calculate Root Mean Squared Error
    rmse <- sqrt(mean((oxytest_pred-oxytest.label)^2))
    message('RMSE: ', rmse)
    
    noaa_future[['temperature']] <- temperaturePrediction
    future_input <- as.matrix(noaa_future[,c("air_temperature", "dayofyear", "temperature")])
    
    forecasted_oxygen = predict(oxybst, future_input)+rnorm(nrow(future_input))*standardDev
    
    # stick bits together                  
    oxygen <- 
      noaa_future |> 
      mutate(site_id = site,
             prediction = forecasted_oxygen,
             variable = "oxygen")
    
    # Generate our labels as the current water temperature
    chlatrain.label = train$chla
    chlatest.label = test$chla
    
    # Convert the input data to a matrix for xgboost
    chlatrain.data = as.matrix(train[, c("air_temperature", "dayofyear", "temperature")])
    chlatest.data = as.matrix(test[, c("air_temperature", "dayofyear", "temperature")])
    
    # Generate Training Input for XGBoost
    chladtrain<-xgb.DMatrix(data = chlatrain.data, label = chlatrain.label)
    # Train our model
    chlabst <- xgboost(data = chladtrain, max.depth = 10, eta = 0.3, nthread = 2, nrounds = 15, verbose = 0)
    # Product Predictions for our Testing Dataset
    chlatest_pred <- predict(chlabst, chlatest.data)
    # Calculate Root Mean Squared Error
    rmse <- sqrt(mean((chlatest_pred-chlatest.label)^2))
    standardDev = sd(chlatest_pred-chlatest.label)
    message('RMSE: ', rmse)
    
    forecasted_chla = predict(chlabst, future_input)+rnorm(nrow(future_input))*standardDev
    
    # stick bits together                  
    chla <- 
      noaa_future |> 
      mutate(site_id = site,
             prediction = forecasted_chla,
             variable = "chla")
    
    forecast <- dplyr::bind_rows(temperature, oxygen, chla)
  }
  
  
  
  # Format results to EFI standard
  forecast <- forecast |>
    mutate(reference_datetime = forecast_date,
           family = "ensemble",
           model_id = model_id) |>
    rename(parameter = ensemble) |>
    select(model_id, datetime, reference_datetime,
           site_id, family, parameter, variable, prediction)
  
  

  
}

### AND HERE WE GO! We're ready to start forecasting ### 

# ## Test with a single site first!
# forecast <- forecast_site( "SUGG" )
# 
# #Visualize the ensemble predictions -- what do you think?
# forecast |> 
#   ggplot(aes(x = datetime, y = prediction, group = parameter)) +
#   geom_line(alpha=0.3) +
#   facet_wrap(~variable, scales = "free") + ggtitle("XGBoost Model")
# 
# ggplot(forecast %>% filter(datetime==noaa_date+20), aes(x=prediction)) + geom_histogram(aes(y=..density..), binwidth = 3) + geom_density(alpha=0.2, fill="#FF6666")

# Run all sites -- may be slow!
forecast <- map_dfr(sites, forecast_site)

#Forecast output file name in standards requires for Challenge.
# csv.gz means that it will be compressed
file_date <- Sys.Date() #forecast$reference_datetime[1]
forecast_file <- paste0("aquatics","-",file_date,"-",model_id,".csv.gz")

#Write csv to disk
write_csv(forecast, forecast_file)

neon4cast::forecast_output_validator(forecast_file)

# Step 4: Submit forecast!

# neon4cast::submit(forecast_file = forecast_file, metadata = NULL, ask = FALSE)

#neon4cast::check_submission(forecast_file)
