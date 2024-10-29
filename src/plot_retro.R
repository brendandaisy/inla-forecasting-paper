library(tidyverse)
library(INLA)
library(lubridate)

##################### data wrangle ##################
prepare_truth_df <- function(data, date_format = "%Y-%m-%d"){ 
  data <- data |>
    mutate(location_number = as.integer(location_number)) |>
    select(location, location_number, date, count) 
  
  return(data)
}

#get files
read_and_combine_csv <- function(directory, file_pattern, strings_to_remove = NULL) {
  setwd(directory)
  
  # List files that match the pattern
  files <- list.files(pattern = file_pattern)
  
  # Optionally remove files/dates 
  if (!is.null(strings_to_remove)) {
    files <- files[!grepl(paste(strings_to_remove, collapse="|"), files)]
  }
  
  
  data_list <- lapply(files, read.csv)
  combined_data <- bind_rows(data_list)
  return(combined_data)
}


#generate dates
generate_dates <- function(start_date, end_date, weeks_increment) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  
  dates_vector <- c()
  current_date <- start_date
  
  while (current_date <= end_date) {
    # Add the formatted current date to the vector
    dates_vector <- c(dates_vector, format(current_date, "%Y-%m-%d"))
    
    # Increment the current date by the specified number of weeks
    current_date <- current_date + weeks(weeks_increment)
  }
  
  return(dates_vector)
}

#################### Plot panel 1 ###########

#for df plot
make_plot_df <- function(data, dates, truth_df) {
  
  # Ensure both location columns are of type character
  data$location <- as.character(data$location)
  truth_df$location <- as.character(truth_df$location)
  
  # Mutate to set values to NA if not in specified dates
  df_processed <- data %>%
    mutate(value = ifelse(reference_date %in% dates, value, NA)) %>%
    rename(date = target_end_date) %>%
    mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
    left_join(truth_df, by = c("date", "location")) %>%
    rename(model = reference_date, prediction = value, quantile = output_type_id) %>%
    select(model, date, horizon, location, prediction, quantile)
  
  # Further processing to take first non-NA value per group
  df_processed <- df_processed %>%
    group_by(location, date, quantile) %>%
    arrange(is.na(prediction)) %>%
    slice(1) %>%
    ungroup()
  
  return(df_processed)
}

generate_plots <- function(disease, plot_df, truth_df) {
  
  # Process the forecast data to keep only specified quantiles
  forecast_df <- plot_df %>%
    spread(quantile, prediction) 
  
  # Define color mapping based on disease
  color_mapping <- list(
    "RSV" = list(ribbon = "#fc9272", line = "#de2d26"),
    "COVID" = list(ribbon = "skyblue", line = "blue"),
    "Flu" = list(ribbon = "#addd8e", line = "#31a354")
  )
  
  colors <- color_mapping[[disease]]
  
  # Create a list to store the plots for each location
  plots_list <- list()
  
  for(curr_location in unique(forecast_df$location)) {
    
    location_forecast <- forecast_df %>%
      filter(location == curr_location, date >= "2022-09-01")
    
    location_truth <- truth_df %>%
      filter(location == curr_location, date >= "2022-09-01")
    
    # Generate a sequence of breaks for the x-axis at the start of each month
    x_breaks <- seq.Date(from = min(location_forecast$date), 
                         to = max(location_forecast$date), 
                         by = "2 months")
    
    # Suppress warnings about NAs in the plot- NAs are supposed to be in there 
    suppressWarnings({
      
      # Create the plot for the current location
      plot <- ggplot() +
        geom_ribbon(data = location_forecast, aes(x = date, ymin = `0.025`, ymax = `0.975`), fill = colors$ribbon, alpha = 0.3) +
        geom_line(data = location_forecast, aes(x = date, y = `0.5`), color = colors$line) +
        geom_point(data = location_truth, aes(x = date, y = count), color = "black") +
        theme_minimal() +
        theme_cowplot() +
        labs(x = "Date", y = "Hospital cases", title = paste(disease, "Retrospective -", curr_location)) +
        theme_minimal()
      
      # Remove background grids for the x-axis and keep y-axis grids
      plot <- plot + scale_x_date(breaks = x_breaks, date_labels = "%Y-%m") + background_grid(major = 'xy', minor = 'y')
      
      # Print the plot
      print(plot)
      
      # Store the plot in the list with the location name as the key
      plots_list[[curr_location]] <- plot
    })
  }
  
  # Return the list of plots
  return(plots_list)
}


generate_plots_save <- function(disease, plot_df, truth_df, pdf_path = NULL, width = 8, height = 6) {
  if (!is.null(pdf_path)) {
    pdf(pdf_path, width = width, height = height)
  }
  
  # Process the forecast data to keep only specified quantiles
  forecast_df <- plot_df %>%
    #filter(quantile %in% c(0.025, 0.25, 0.5, 0.75, 0.975),
    #location_name == "California") %>%
    spread(quantile, prediction) 
  
  # Define color mapping based on disease
  color_mapping <- list(
    "RSV" = list(ribbon = "#fc9272", line = "#de2d26"),
    "COVID" = list(ribbon = "skyblue", line = "blue"),
    "Flu" = list(ribbon = "#addd8e", line = "#31a354")
  )
  
  colors <- color_mapping[[disease]]
  
  # Loop through each unique location
  for(curr_location in unique(forecast_df$location)) {
    
    location_forecast <- forecast_df %>%
      filter(location == curr_location, date >= "2022-09-01")
    
    location_truth <- truth_data %>%
      filter(location == curr_location, date >= "2022-09-01")
    
    # Generate a sequence of breaks for the x-axis at the start of each month
    x_breaks <- seq.Date(from = min(location_forecast$date), 
                         to = max(location_forecast$date), 
                         by = "2 months")
    
    # Suppress warnings about NAs in the plot- NAs are supposed to be in there 
    suppressWarnings({
      
      # Create the plot for the current location
      plot <- ggplot() +
        geom_ribbon(data = location_forecast, aes(x = date, ymin = `0.025`, ymax = `0.975`), fill = colors$ribbon, alpha = 0.3) +
        geom_line(data = location_forecast, aes(x = date, y = `0.5`), color = colors$line) +
        geom_point(data = location_truth, aes(x = date, y = count), color = "black") +
        theme_minimal() +
        theme_cowplot() +
        labs(x = "Date", y = "hospital cases", title = paste(disease, "Retrospective -", curr_location)) +
        theme_minimal()
      
      # Remove background grids for the x-axis and keep y-axis grids
      plot <- plot + scale_x_date(breaks = x_breaks, date_labels = "%Y-%m") + background_grid(major = 'xy', minor = 'y')
      
      print(plot)
    })
  }
  # Close the PDF device if it was opened
  if (!is.null(pdf_path)) {
    dev.off()
  }
}

#################### Plot panel 2 ###########

# Grab retro data saved in folder- multiple folders + models  
process_data_folders <- function(folders, model_names, fullpath) {
  data_frames <- list()
  
  if (length(folders) != length(model_names)) {
    stop("The number of folders and model names must match.")
  }
  
  for (i in seq_along(folders)) {
    folder <- folders[i]
    model_name <- model_names[i]
    
    full_path <- paste0("../processed_data/", folder)
    
    files <- list.files(full_path, pattern = "\\d{4}-\\d{2}-\\d{2}-UGA_flucast-INFLAenza\\.csv", full.names = TRUE)
    
    data_list <- lapply(files, read_csv)
    
    data_combined <- bind_rows(data_list) %>%
      mutate(model = model_name)
    
    data_frames[[folder]] <- data_combined
  }
  
  return(data_frames)
}

#Grab retro data saved in folder- single folder for single model 
read_and_combine_csv <- function(directory, file_pattern, strings_to_remove = NULL) {
  setwd(directory)
  
  # List files that match the pattern
  files <- list.files(pattern = file_pattern)
  
  # Optionally remove files/dates 
  if (!is.null(strings_to_remove)) {
    files <- files[!grepl(paste(strings_to_remove, collapse="|"), files)]
  }
  
  
  data_list <- lapply(files, read.csv)
  combined_data <- bind_rows(data_list)
  return(combined_data)
}

#Format for plotting step 
format_dataframe <- function(df, truth_df) {
  if ("location_number" %in% colnames(df)) {
    # If location_number is in the dataframe, join on location_number
    formatted_df <- df |>
      mutate(location_number = as.integer(location_number)) |>
      rename(date = target_end_date) |>
      mutate(date = as.Date(date, format = "%Y-%m-%d")) |>
      left_join(truth_df, by = c("date", "location_number")) |>
      rename(prediction = value, true_value = count, quantile = output_type_id) |>
      mutate(
        quantile = as.numeric(as.character(quantile))
      ) |>
      select(model, date, location, horizon, prediction, true_value, quantile)
  } else if ("location" %in% colnames(df)) {
    # If location is in the dataframe, join on location
    formatted_df <- df |>
      rename(date = target_end_date) |>
      mutate(date = as.Date(date, format = "%Y-%m-%d")) |>
      left_join(truth_df, by = c("date", "location")) |>
      rename(prediction = value, true_value = count, quantile = output_type_id) |>
      mutate(
        quantile = as.numeric(as.character(quantile))
      ) |>
      select(model, date, location, horizon, prediction, true_value, quantile)
  } else {
    stop("Neither 'location' nor 'location_number' found in the dataframe.")
  }
  
  return(formatted_df)
}

#Allow baseline to be able to be used with format_dataframe
quantile_filter <-c(0.025, 0.25, 0.5, 0.75, 0.975)

clean_baseline <- function(bdata, disease, quantile_filter) {
  bdata <- bdata %>%
    mutate(model = paste0("baseline_", disease)) %>%
    rename(location_number = location) %>%
    filter(output_type_id %in% quantile_filter)
  
  return(bdata)
}

# For removing National which isn't needed for this plot 
remove_na_location <- function(df) {
  df <- df[!is.na(df$location), ]
  return(df)
}

#Make Relative WIS 
make_rwis <- function(df, baseline_model) {
  df <- df %>%
    group_by(location) %>%
    mutate(
      baseline_wis = interval_score[model == baseline_model],
      r_wis = ifelse(model == baseline_model, 1, interval_score / baseline_wis)
    ) %>%
    ungroup() %>%
    select(-baseline_wis) # remove temp column
  
  return(df)
}

#Grab data for making map
create_map_df <- function(map_df) {
  # Get US state boundaries including Puerto Rico
  states_sf <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
  
  # Fix state names to match with your dataframe
  states_sf <- states_sf %>%
    mutate(ID = tools::toTitleCase(ID))
  
  # Merge your data with the spatial data
  merged_sf <- states_sf %>%
    left_join(map_df, by = c("ID" = "location"))
  
  return(merged_sf)
}

#plot
generate_map <- function(data, title = "Relative WIS by State (RSV)") {
  p <- ggplot(data = data) +
    geom_sf(aes(fill = r_wis)) +
    geom_sf_text(aes(label = round(r_wis, 2)), size = 2, color = "white") +
    scale_fill_viridis_c(option = "viridis") +
    theme_minimal() +
    labs(title = title, fill = "Relative WIS") +
    theme_map()
  
  print(p)
}
