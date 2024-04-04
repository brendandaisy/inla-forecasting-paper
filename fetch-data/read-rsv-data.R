library(tidyverse)
library(RSocrata)

#with removed covid season
fetch_rsv_full <- function(na.rm = TRUE) {
  raw_csv_url <- "https://raw.githubusercontent.com/HopkinsIDD/rsv-forecast-hub/main/target-data/2024-03-08_rsvnet_hospitalization.csv" 
  #archive_csv_url <- "https://raw.githubusercontent.com/HopkinsIDD/rsv-forecast-hub/main/target-data/archive/2024-03-01_rsvnet_hospitalization.csv" 
  locations <- readr::read_csv(file = "https://raw.githubusercontent.com/HopkinsIDD/rsv-forecast-hub/main/auxiliary-data/location_census/locations.csv", col_select= 2:3)
  
  health_data <- readr::read_csv(file = raw_csv_url, na = ifelse(na.rm, "NA", ""))
  
  health_data %>%
    mutate(date = as.Date(date, format="%Y-%m-%d"),
           population = as.numeric(population),
           epiweek = lubridate::epiweek(date),
           value = ceiling(value)) %>% # Round the 'value' column
    filter(date < as.Date("2020-10-03") | date > as.Date("2021-10-02")) %>%
    filter(!is.na(value)) %>%
    filter(age_group == "0-130") %>%
    filter (target == "inc hosp") %>%
    left_join(locations, by = "location") 
}

fetch_rsv_2021 <- function(na.rm = TRUE) {
  raw_csv_url <- "https://raw.githubusercontent.com/HopkinsIDD/rsv-forecast-hub/main/target-data/2024-02-16_rsvnet_hospitalization.csv" 
  locations <- readr::read_csv(file = "https://raw.githubusercontent.com/HopkinsIDD/rsv-forecast-hub/main/auxiliary-data/location_census/locations.csv", col_select= 2:3)
  
  health_data <- readr::read_csv(file = raw_csv_url, na = ifelse(na.rm, "NA", ""))
  
  health_data %>%
    mutate(date = as.Date(date, format="%Y-%m-%d"),
           population = as.numeric(population),
           epiweek = lubridate::epiweek(date),
           value = ceiling(value)) %>% # Round the 'value' column
    filter(date >= as.Date('2021-05-01')) %>%  
    filter(!is.na(value)) %>%
    filter(age_group == "0-130") %>%
    filter (target == "inc hosp")  %>%
    left_join(locations, by = "location")  
}


#removes bad season from covid

#10/03/2020
#10/02/2021
