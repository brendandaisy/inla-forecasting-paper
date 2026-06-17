library(tidyverse)
library(RSocrata)
library(lubridate)
library(INLA)
library(arrow)

source("../src/model-formulas.R")
source("../src/fit-inla-model.R")

inla.setOption(inla.mode="classic")

##################
#Functions
##################
prep_fit_data_rsv <- function(rsv, weeks_ahead=4, ex_lam=pop_served) { 
  #To Do: add Extend of population column 
  
  date_seq <- seq.Date(min(rsv$date), max(rsv$date), "1 week")
  
  date_ind <- tibble(t=seq_along(date_seq), date=date_seq) |> 
    filter(date %in% unique(rsv$date))
  
  ret <- rsv |> 
    left_join(date_ind, by=c("date")) |> 
    mutate(
      iloc=as.numeric(fct_inorder(location)), # INLA needs groups as ints starting from 1, so add numeric state code
      ex_lam={{ex_lam}}
    )
  
  if (weeks_ahead > 0) {
    # make a dataframe to hold group info for forecasting
    pred_df <- expand_grid( # makes pairs of new weeks X location
      tibble(
        date=weeks(1:weeks_ahead) + max(ret$date),  # Fixed here
        t=1:weeks_ahead + max(ret$t),
        epiweek=epiweek(date)
      ),
      distinct(ret, location, iloc)
    )
    
    # go and find most recent population/ex_lam values for each state, in case they change
    location_data <- ret |> 
      slice_max(date) |> 
      distinct(location, ex_lam)
    
    # add on most recent pop/ex_lam values
    pred_df <- pred_df |> 
      left_join(
        location_data, by=c("location"), unmatched="error", relationship="many-to-one"
      )
    
    ret <- bind_rows(ret, pred_df) # add to data for counts to be NAs
  }
  
  return(ret |> 
           mutate(t2=t) |> 
           arrange(t, iloc)
  )
}


sample_national_test3 <- function(fit_df, fit, forecast_date, nsamp = 10000) {
  # Get distinct locations and ex_lam values
  state_info <- distinct(fit_df, location, ex_lam, population) |>
    tail(12)  # Handle multiple population numbers here if necessary
  nstate <- nrow(state_info)
  
  # Prepare the return dataframe with population summed for all states
  ret_df <- fit_df |>
    filter(date > forecast_date) |>
    group_by(date, t, epiweek) |>
    summarise(population = sum(population), .groups = "drop")
  
  # Generate sampled values from the fit model
  jsamp_fvals <- exp(inla.rjmarginal(nsamp, fit$selection)$samples)
  ex_lam <- filter(fit_df, date >= forecast_date)$ex_lam
  
  # Generate sequence for indexing per state
  tslice <- map(1:nrow(ret_df), ~nstate * (.x - 1) + 1:nstate)
  
  # Now calculate the total population proportions
  total_us_population <- 345825752  # Assuming the total US population is fixed
  total_states_population <- sum(state_info$population)  # Total population of the states in the dataset
  total_exlam <- sum(state_info$ex_lam)
  
  # Calculate the refined scaling factor
  population_proportion <- total_states_population / total_us_population
  refined_scaling_factor <- 1 / population_proportion
  
  # Modify the inner function to scale using the final scaling formula
  imap_dfr(tslice, \(idx, t) {
    nat_sum_per_t <- map_dbl(1:nsamp, \(samp) {
      lambda <- jsamp_fvals[idx, samp] * ex_lam
      samp <- rpois(nstate, lambda)
      
      # Calculate the total sum for the states
      state_sums <- sum(samp)
      
      # Apply the final scaling formula directly: (value / pop_served) * population
      #scaled_sum <- (state_sums * refined_scaling_factor/ total_exlam) * total_states_population #way too big/ removing * refined does nothing 
      #scaled_sum <- (state_sums * refined_scaling_factor/ total_exlam) # way too small
      
      #return(scaled_sum)
      return(state_sums)
    })
    
    tibble_row(location = "US", count_samp = list(nat_sum_per_t))
  }) |>
    bind_cols(ret_df) |>
    select(date:epiweek, location, population, count_samp)
}

summarize_quantiles_test <- function(pred_samples, nat_samps, forecast_date, q) {
  pred_samples |> 
    filter(date > forecast_date) |> 
    bind_rows(nat_samps) |> 
    unnest(count_samp) |> 
    group_by(date, location) |> 
    summarize(qs = list(value = quantile(count_samp, probs=q)), .groups="drop") |> 
    unnest_wider(qs) |>
    pivot_longer(contains("%"), names_to="quantile") |> 
    mutate(quantile = as.numeric(gsub("[\\%,]", "", quantile))/100) |> 
    mutate(target = paste0("inc hosp"),
           #horizon=as.numeric(as.factor(date)) - 1,
           horizon=as.numeric(as.factor(date)),
           reference_date = forecast_date,
           #target_end_date = forecast_date + horizon*7,
           output_type_id = quantile,
           output_type = 'quantile',
           value = round(value)) %>%
    arrange(location, horizon, quantile) |>
    rename(target_end_date = date) %>%
    dplyr::select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value)
}

sample_count_predictions <- function(fit_df, fit, forecast_date, nsamp=1000) {
  nloc <- length(unique(fit_df$iloc))
  
  ret_df <- fit_df |>
    filter(date > forecast_date)
  
  jsamp_fvals <- exp(inla.rjmarginal(nsamp, fit$selection)$samples)
  ex_lam <- ret_df$ex_lam
  
  count_samp <- list_transpose(map(1:nsamp, \(samp) { # invert list so have sampled counts for each row
    lambda <- jsamp_fvals[,samp] * ex_lam
    rpois(length(lambda), lambda) # indep. Poisson for each spacetime series
  }))
  
  ret_df |> 
    mutate(count_samp=count_samp) |> 
    select(where(~!any(is.na(.x))))
}


generate_weekly_dates <- function(start_date, end_date) {
  # Convert input dates to Date format
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  
  # Generate sequence of dates, each 7 days apart
  seq_dates <- seq(from = start_date, to = end_date, by = "week")
  
  return(seq_dates)
}

start_date <- "2022-09-10"
#end_date <- "2022-09-30"
end_date <- "2024-04-13"
weekly_dates <- generate_weekly_dates(start_date, end_date)
#print(weekly_dates)


#######################
#start main function + data read
########################

disease_df2 <- read_csv("../data/weekly-rsv-us.csv")


library(dplyr)
library(readr)
library(lubridate)

generate_save_weekly_forecasts_rsv_folder2 <- function(forecast_dates, disease_df, folder_name) {
  # Read the population served file once
  name_key <- read_csv("../data/name_key.csv")
  pop_served <- read_csv("../data/pop_served.csv")
  
  pop_served <- pop_served |>
    select(location, location_name, pop_served, population, percentage_participation)
  
  
  disease_df <- disease_df |>
    rename(location_name = location) |>
    left_join(name_key, by = "location_name")
  
  pred_quantiles <- map(forecast_dates, \(fd) {
    fit_df <- disease_df |> 
      filter(date <= fd) |>
      prep_fit_data_rsv(weeks_ahead = 4, ex_lam = population)
    
    # Define the model
    model <- model_formula("shared", temporal = "ar1", spatial = "exchangeable")
    
    
    fit_df <- fit_df %>%
      group_by(location) %>%
      mutate(population = ifelse(is.na(population), lag(population), population)) %>%
      fill(population, .direction = "down") %>%
      ungroup()
    
    fit <- fit_inla_model(fit_df, forecast_date = fd, model, pc_prior_u = c(1, 1))
    
    nat_samps <- sample_national_test3(fit_df, fit, forecast_date = fd, nsamp = 5000)
    pred_samples <- sample_count_predictions(fit_df, fit, forecast_date = fd, nsamp = 5000)
    
    q <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
    
    cleaned_forecasts_quantiles <- summarize_quantiles_test(pred_samples, nat_samps, forecast_date = fd, q = q)
    
    unscaled_forecasts_quantiles <- cleaned_forecasts_quantiles |> 
      left_join(pop_served, by = "location") |> 
      mutate(unscaled_value = (value/pop_served) * population)
    
    # Save the unscaled version as a CSV
    unscaled_forecasts_quantiles |> 
      mutate(output_type_id = as.double(output_type_id),
             horizon = as.integer(horizon)) |> 
      relocate(location, .after = horizon) 
    
    list(pred = unscaled_forecasts_quantiles, fit_df = fit_df, fit = fit)
  })
  
  pred_list <- lapply(pred_quantiles, function(x) x$pred)
  
  # Add the reference_date to each dataframe in the list
  pred_list <- lapply(pred_list, function(df) {
    reference_date <- df$date[1]  
    #df$reference_date <- ymd(reference_date) - 7
    return(df)
  })
  
  # Create output directory based on folder_name
  output_dir <- paste0("../processed_data/", folder_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save each dataframe in result_list to a CSV file
  lapply(seq_along(pred_list), function(i) {
    df <- pred_list[[i]]
    ref_date <- as.character(df$reference_date[1])  # Extract the first reference_date
    file_path <- paste0(output_dir, "/", ref_date, "-UGA_flucast-INFLAenza.csv")
    write_csv(df, file_path)
  })
  
  return("All files processed and saved")
}


#################### Produce actual Forecast ####################
weekly_dates |> 
  map(generate_save_weekly_forecasts_rsv_folder2, 
      disease_df = disease_df2, 
      folder_name = "rsv_graphing_dec10") ##//// THIS CREATES A FOLDER


##############################################
#Reformat forecast (use these ones for figure)
###############################################
process_single_csv2 <- function(file_path) {
  # Read the CSV file
  df <- read_csv(file_path)
  
  # Apply the transformations
  df_processed <- df |> 
    mutate(unscaled_value = round(unscaled_value)) |> 
    select(reference_date, target, target_end_date, horizon, location, location_name, output_type, output_type_id, unscaled_value) |> 
    rename(value = unscaled_value,
           location_number = location,
           location = location_name)
  
  
  return(df_processed)
}

# Define the function to process a folder of CSVs
process_folder_of_csvs <- function(input_dir, output_dir) {
  # Get a list of CSV files in the input directory
  csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Apply the transformation to each CSV and save the results
  lapply(seq_along(csv_files), function(i) {
    # Process the CSV
    df <- process_single_csv2(csv_files[[i]])
    
    # Get the reference date for the file name
    ref_date <- as.character(df$reference_date[1])
    
    # Create the output file path
    file_path <- paste0(output_dir, "/", ref_date, "-UGA_flucast-INFLAenza.csv")
    
    # Save the processed CSV
    write_csv(df, file_path)
  })
}

process_folder_of_csvs("../processed_data/rsv_graphing_dec10/", "../processed_data/rsv_graphing_dec10_clean")

