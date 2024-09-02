
library(tidyverse)
library(INLA)
library(sf)
library(spdep)
library(lubridate)


##################### data wrangle #####
filter_by_epiweek <- function(data, start_epiweek, end_epiweek, start_date) {
  data %>%
    distinct(date) %>%
    mutate(epiweek = epiweek(date)) %>%
    filter((epiweek >= start_epiweek | epiweek < end_epiweek) & date >= start_date) %>%
    pull(date) %>%
    unique()
}

read_data_model <- function(disease) {
  if (disease == "rsv_2021") {
    return(fetch_rsv_2021())
  } else if (disease == "flu_2021") {
    return(fetch_flu_2021())
  } else if (disease == "rsv_2018") {
    return(fetch_rsv_2018())
  } else if (disease == "covid_2020") {
    return(fetch_covid_2020())
  } else if (disease == "covid_2021") {
    return(fetch_covid_2021())
  } else if (disease == "rsv_2022") {
    return(fetch_rsv_2022())
  } else if (disease == "flu_2022") {
    return(fetch_flu_2022()) 
  } else if (disease == "covid_2022") {
    return(fetch_covid_2022())   
  } else {
    stop("Invalid input: please specify 'rsv_2021', 'flu_2021', 'covid_2021','rsv_2022', 'flu_2022', 'covid_2022'")
  }
}

fetch_rsv_2018 <-  function(na.rm = TRUE) { 
  rsv <-read_csv("../data/weekly-rsv-us.csv")
  
  covid_interval <- interval(ymd("2020-03-01"), ymd("2021-06-01"))
  rsv <- rsv |> 
    #filter(date >=  as.Date("2021-09-01"))
   filter(!(date %within% covid_interval))
  
  location <- read_csv("../data/location-info.csv") |>
    select(-population) |>
    rename(location_number = location,
           location = location_name)
  
  disease <- left_join(rsv, location, by = "location")
  return(disease)
}

fetch_covid_2020 <- function(na.rm = TRUE) { 
  covid <- read_csv("../data/weekly-covid-us.csv")
  
  location <- read_csv("../data/location-info.csv") |>
    select(-population) |>
    rename(location_number = location,
           location = location_name)
  
  disease <- left_join(covid, location, by = "location")
  return(disease)
}

fetch_rsv_2021 <-  function(na.rm = TRUE) { 
  rsv <-read_csv("../data/weekly-rsv-us.csv")
  
  #covid_interval <- interval(ymd("2020-03-01"), ymd("2021-06-01"))
  rsv <- rsv |> 
    filter(date >=  as.Date("2021-09-01"))
    #filter(!(date %within% covid_interval))
  
  location <- read_csv("../data/location-info.csv") |>
    select(-population) |>
    rename(location_number = location,
           location = location_name)
  
  disease <- left_join(rsv, location, by = "location")
  return(disease)
}

fetch_flu_2021 <-  function(na.rm = TRUE) { 
  flu <-read_csv("../data/weekly-flu-us.csv")
  
  #covid_interval <- interval(ymd("2020-03-01"), ymd("2021-06-01"))
  flu <- flu |> 
    filter(date >=  as.Date("2021-09-01"))
    #filter(!(date %within% covid_interval))
  
  location <- read_csv("../data/location-info.csv") |>
    select(-population) |>
    rename(location_number = location,
           location = location_name)
  
  disease <- left_join(flu, location, by = "location")
  return(disease)
}
  
fetch_covid_2021 <- function(na.rm = TRUE) { 
  covid <- read_csv("../data/weekly-covid-us.csv")
  
  covid <- covid |>
    filter(date >=  as.Date("2021-09-01"))
  
  location <- read_csv("../data/location-info.csv") |>
    select(-population) |>
    rename(location_number = location,
           location = location_name)
  
  disease <- left_join(covid, location, by = "location")
  return(disease)
}

fetch_rsv_2022 <-  function(na.rm = TRUE) {
  rsv <-read_csv("../data/weekly-rsv-us.csv")
  
  rsv <- rsv |> 
    filter(date >=  as.Date("2022-02-02"))
  
  location <- read_csv("../data/location-info.csv") |>
    select(-population) |>
    rename(location_number = location,
           location = location_name)
  
  disease <- left_join(rsv, location, by = "location")
  return(disease)
}

fetch_flu_2022 <-  function(na.rm = TRUE) { 
  flu <-read_csv("../data/weekly-flu-us.csv")
  
  flu <- flu |> 
    filter(date >=  as.Date("2022-02-02"))
  
  location <- read_csv("../data/location-info.csv") |>
    select(-population) |>
    rename(location_number = location,
           location = location_name)
  
  disease <- left_join(flu, location, by = "location")
  
  return(disease)
}

fetch_covid_2022 <-  function(na.rm = TRUE) { 
  covid <- read_csv("../data/weekly-covid-us.csv")
  
  covid <- covid |> 
    filter(date >=  as.Date("2022-02-02"))
  
  location <- read_csv("../data/location-info.csv") |>
    select(-population) |>
    rename(location_number = location,
           location = location_name)
  
  disease <- left_join(covid, location, by = "location")
  return(disease)
}

hub_formatting <- function(pred_samples) { 

  pred_samples$date <- ymd(pred_samples$date)
  
  pred_samples %>%
    group_by(location) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(
      output_type_id = quantile,
      target = "inc hosp",
      #horizon = as.integer(date - min(date)) + 1,  
      horizon=as.numeric(as.factor(date)),
      output_type = "quantile",
      value = round(value, digits = 0)
    ) %>%
    rename(target_end_date = date) %>%
    select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value, mean)
}

################# 
generate_save_weekly_forecast_rsv <- function(forecast_dates, disease_df, folder_name) {

  # Define the model
  model <- model_formula("shared", temporal="ar1", spatial="exchangeable")
  
  # Fit the model and sample predictions for each of the timepoints
  pred_quantiles <- map(forecast_dates, \(fd) {
    fit_df <- disease_df |> 
      filter(date < fd) |> 
      prep_fit_data_rsv(weeks_ahead=4, ex_lam=pop_served)
    
    fit <- fit_inla_model(fit_df, fd, model, pc_prior_u=c(1, 1))
    
    pred_summ <- sample_count_predictions(fit_df, fit, fd, nsamp=5000) |> 
      summarize_quantiles()  # Assuming summarize_quantiles is a defined function
    
    list(pred=pred_summ, fit_df=fit_df, fit=fit)
  })
  
  pred_list <- lapply(pred_quantiles, function(x) x$pred)
  
  # Add the reference_date to each dataframe in the list
  pred_list <- lapply(pred_list, function(df) {
    reference_date <- df$date[1]  # Extract the first date from the 'date' column
    df$reference_date <- ymd(reference_date) - 7
    return(df)
  })
  
  # Apply summarize_quantiles_test to each item in the list
  result_list <- lapply(pred_list, hub_formatting)
  
  # Create output directory based on folder_name
  output_dir <- paste0("../processed_data/", folder_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save each dataframe in result_list to a CSV file
  lapply(seq_along(result_list), function(i) {
    df <- result_list[[i]]
    ref_date <- as.character(df$reference_date[1])  # Extract the first reference_date
    file_path <- paste0(output_dir, "/", ref_date, "-UGA_flucast-INFLAenza.csv")
    write_csv(df, file_path)
  })
  
  return(paste("All files saved in", output_dir))
}

##### Covid and flu

generate_save_weekly_forecast <- function(forecast_dates, disease_df, folder_name) {
  
  # Define the model
  model <- model_formula("shared", temporal="ar1", spatial="exchangeable")
  
  # Fit the model and sample predictions for each of the timepoints
  pred_quantiles <- map(forecast_dates, \(fd) {
    fit_df <- disease_df |> 
      filter(date < fd) |> 
      prep_fit_data_flu_covid(weeks_ahead=4, ex_lam=population)
    
    fit <- fit_inla_model(fit_df, fd, model, pc_prior_u=c(1, 1))
    
    pred_summ <- sample_count_predictions(fit_df, fit, fd, nsamp=5000) |> 
      summarize_quantiles()  # Assuming summarize_quantiles is a defined function
    
    list(pred=pred_summ, fit_df=fit_df, fit=fit)
  })
  
  pred_list <- lapply(pred_quantiles, function(x) x$pred)
  
  # Add the reference_date to each dataframe in the list
  pred_list <- lapply(pred_list, function(df) {
    reference_date <- df$date[1]  
    df$reference_date <- ymd(reference_date) - 7
    return(df)
  })
  
  # Apply summarize_quantiles_test to each item in the list
  result_list <- lapply(pred_list, hub_formatting)
  
  # Create output directory based on folder_name
  output_dir <- paste0("../processed_data/", folder_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save each dataframe in result_list to a CSV file
  lapply(seq_along(result_list), function(i) {
    df <- result_list[[i]]
    ref_date <- as.character(df$reference_date[1])  # Extract the first reference_date
    file_path <- paste0(output_dir, "/", ref_date, "-UGA_flucast-INFLAenza.csv")
    write_csv(df, file_path)
  })
  
  return(paste("All files saved in", output_dir))
}

##### Covid_flu_besag

generate_save_weekly_forecast_besag <- function(forecast_dates, disease_df, folder_name) {

  graph <- load_us_graph(disease_df) |> sf2mat()
  
  # Define the model
  model <- model_formula("shared", temporal="ar1", spatial="besagproper")
  
  # Fit the model and sample predictions for each of the timepoints
  pred_quantiles <- map(forecast_dates, \(fd) {
    fit_df <- disease_df |> 
      filter(date < fd) |> 
      prep_fit_data_flu_covid(weeks_ahead=4, ex_lam=population)
    
    fit <- fit_inla_model(fit_df, fd, model, pc_prior_u=c(1, 1), graph=graph)
    
    pred_summ <- sample_count_predictions(fit_df, fit, fd, nsamp=5000) |> 
      summarize_quantiles()  
    
    list(pred=pred_summ, fit_df=fit_df, fit=fit)
  })
  
  pred_list <- lapply(pred_quantiles, function(x) x$pred)
  
  # Add the reference_date to each dataframe in the list
  pred_list <- lapply(pred_list, function(df) {
    reference_date <- df$date[1]  
    df$reference_date <- ymd(reference_date) - 7
    return(df)
  })
  
  # Apply summarize_quantiles_test to each item in the list
  result_list <- lapply(pred_list, hub_formatting)
  
  # Create output directory based on folder_name
  output_dir <- paste0("../processed_data/", folder_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save each dataframe in result_list to a CSV file
  lapply(seq_along(result_list), function(i) {
    df <- result_list[[i]]
    ref_date <- as.character(df$reference_date[1])  # Extract the first reference_date
    file_path <- paste0(output_dir, "/", ref_date, "-UGA_flucast-INFLAenza.csv")
    write_csv(df, file_path)
  })
  
  return(paste("All files saved in", output_dir))
}

##### RSV_besag
generate_save_weekly_forecast_besag_rsv <- function(forecast_dates, disease_df, folder_name) {
  
  graph <- load_us_graph(disease_df) |> sf2mat()
  
  # Define the model
  model <- model_formula("shared", temporal="ar1", spatial="besagproper")
  
  # Fit the model and sample predictions for each of the timepoints
  pred_quantiles <- map(forecast_dates, \(fd) {
    fit_df <- disease_df |> 
      filter(date < fd) |> 
      prep_fit_data_rsv(weeks_ahead=4, ex_lam=pop_served)
    
    fit <- fit_inla_model(fit_df, fd, model, pc_prior_u=c(1, 1), graph=graph)
    
    pred_summ <- sample_count_predictions(fit_df, fit, fd, nsamp=5000) |> 
      summarize_quantiles()  
    
    list(pred=pred_summ, fit_df=fit_df, fit=fit)
  })
  
  pred_list <- lapply(pred_quantiles, function(x) x$pred)
  
  # Add the reference_date to each dataframe in the list
  pred_list <- lapply(pred_list, function(df) {
    reference_date <- df$date[1]  
    df$reference_date <- ymd(reference_date) - 7
    return(df)
  })
  
  # Apply summarize_quantiles_test to each item in the list
  result_list <- lapply(pred_list, hub_formatting)
  
  # Create output directory based on folder_name
  output_dir <- paste0("../processed_data/", folder_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save each dataframe in result_list to a CSV file
  lapply(seq_along(result_list), function(i) {
    df <- result_list[[i]]
    ref_date <- as.character(df$reference_date[1])  # Extract the first reference_date
    file_path <- paste0(output_dir, "/", ref_date, "-UGA_flucast-INFLAenza.csv")
    write_csv(df, file_path)
  })
  
  return(paste("All files saved in", output_dir))
}

############################################### Baseline.  ######################################

##############################
##     Helper functions     ##
##############################

adjust_data <- function(data) {
  data |>
    rename(location_name = location,
           location = location_number,
           value = count) |>
    select(-abbreviation)
}

# Move to the closest saturday date
curr_else_next_date_with_ltwday <- function(date, ltwday) {
  assert_class(date, "Date")
  assert_integerish(ltwday, lower = 0L, upper = 7L)
  #
  date + (ltwday - as.POSIXlt(date)$wday) %% 7L
}

location_to_abbr <- function(location) {
  dictionary <-
    state_census %>%
    mutate(fips = sprintf("%02d", fips)) %>%
    transmute(
      location = case_match(fips, "00" ~ "US", .default = fips),
      abbr
    )
  dictionary$abbr[match(location, dictionary$location)]
}

##############################
##    Prepare input data    ##
##############################

prepare_epi_data <- function(rsv) {
  # Load required libraries
  library(dplyr)
  library(epiprocess)
  
  # Prepare target table
  target_tbl <- rsv %>%
    mutate(date = as.Date(date, format = "%Y-%m-%d"),
           location = as.character(location),
           location_name = as.character(location_name),
           value = as.numeric(value),
           weekly_rate = as.numeric(weekly_rate)) %>%
    select(-population)
  
  # Transform to epi_df from epiprocess
  target_edf <- target_tbl %>%
    transmute(
      geo_value = location_to_abbr(location),
      time_value = date,
      weekly_count = value
    ) %>%
    epiprocess::as_epi_df()
  
  return(target_edf)
}


######################
## Prepare baseline ##
######################
generate_forecasts_for_date <- function(reference_date, target_edf) {
  desired_max_time_value <- reference_date - 7L
  rng_seed <- as.integer((59460707 + as.numeric(reference_date)) %% 2e9)
  withr::with_rng_version("4.0.0", withr::with_seed(rng_seed, {
    fcst <- epipredict::cdc_baseline_forecaster(
      target_edf %>%
        filter(time_value >= as.Date("2021-12-04")) %>%
        filter(time_value <= desired_max_time_value),
      "weekly_count",
      epipredict::cdc_baseline_args_list(aheads = 1:5, nsims = 1e5)
    )
    preds <- fcst$predictions %>%
      mutate(
        forecast_date = reference_date,
        ahead = as.integer(.data$target_date - reference_date) %/% 7L
      )
    # Removed the part that binds the most recent observed data with an ahead of -1
  }))
  return(preds)
}


###################
## Format, write ##
###################
library(dplyr)
library(readr) 

format_and_save_forecasts <- function(preds, reference_date, output_dirpath) {
  # Format predictions
  preds_formatted <- preds %>%
    flusight_hub_formatter(
      target = "wk inc hosp",
      output_type = "quantile"
    ) %>%
    drop_na(output_type_id) %>%
    filter(horizon != 0) %>%
    arrange(target, horizon, location) %>%
    select(
      reference_date, horizon, target, target_end_date, location,
      output_type, output_type_id, value
    ) 
  
  # Ensure output directory exists
  if (!dir.exists(output_dirpath)) {
    dir.create(output_dirpath, recursive = TRUE)
  }
  
  # Define the path for the output file
  output_filepath <- file.path(output_dirpath, sprintf("%s-baseline.csv", reference_date))
  
  # Save formatted predictions to CSV
  write_csv(preds_formatted, output_filepath)
  
  return(output_filepath)
}


