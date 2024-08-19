### most recent version of data fetching as of June 2024
### fetch US data for the 3 diseases and clean them to a standardized format
### save the CSV so we are working from a static product for the INLA forecasting paper
library(tidyverse)
library(RSocrata)
library(lubridate)

save_weekly_flu <- function(hd) {
    locations <- read_csv(file="https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/locations.csv", col_select= 1:4)
    
    #remove  VI and AS as they are not included for FluSight, keep only necessary vars and add epiweek and epiyear 
    flu <- hd |> 
        filter(!state %in% c("VI", "AS")) |> 
        select(state, date, epiweek, epiyear, count=previous_day_admission_influenza_confirmed) |> 
        mutate(count=parse_number(count))
    
    #summarize US Flu 
    flu_us <- flu |> 
        group_by(date, epiweek, epiyear) |> 
        summarise(count=sum(count, na.rm=TRUE)) |> # fair to ignore the few NA days in 2021
        mutate(state="US") |> 
        ungroup()
    
    #bind state and US data
    flu <- flu |> 
        bind_rows(flu_us) |> 
        left_join(locations, by=join_by("state" == "abbreviation"))
    
    #convert counts to weekly and calculates weekly rate 
    flu_wk <- flu |> 
        group_by(location_name, epiweek, epiyear, population) |> 
        summarise(
            date=max(date),
            count=sum(count, na.rm=TRUE),
            .groups="drop"
        ) |> 
        mutate(weekly_rate=(count*100000)/population)
    
    flu_wk |> 
        rename(location=location_name) |> 
        relocate(location, date, epiweek, epiyear, population, count, weekly_rate) |> 
        arrange(date, location) |> # can't hurt to just make sure everything is ordered this way
        write_csv("data/weekly-flu-us.csv")
}

save_weekly_covid <- function(hd) {
    locations <- read_csv(file="https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/locations.csv", col_select= 1:4)
    
    #remove  VI and AS as they are not included for FluSight, keep only necessary vars and add epiweek and epiyear 
    covid <- hd |> 
        filter(!state %in% c("VI", "AS")) |> 
        mutate(
            adult=parse_number(previous_day_admission_adult_covid_confirmed),
            pediatric=parse_number(previous_day_admission_pediatric_covid_confirmed),
            count=adult + pediatric
        ) |> 
        select(state, date, epiweek, epiyear, count)
    
    #summarize US Flu 
    covid_us <- covid |> 
        group_by(date, epiweek, epiyear) |> 
        summarise(count=sum(count, na.rm=TRUE)) |> # fair to ignore the few NA days in 2021
        mutate(state="US") |> 
        ungroup()
    
    #bind state and US data
    covid <- covid |> 
        bind_rows(covid_us) |> 
        left_join(locations, by=join_by("state" == "abbreviation"))
    
    #convert counts to weekly and calculates weekly rate 
    covid_wk <- covid |> 
        group_by(location_name, epiweek, epiyear, population) |> 
        summarise(
            date=max(date),
            count=sum(count, na.rm=TRUE),
            .groups="drop"
        ) |> 
        mutate(weekly_rate=(count*100000)/population)
    
    covid_wk |> 
        rename(location=location_name) |> 
        relocate(location, date, epiweek, epiyear, population, count, weekly_rate) |> 
        arrange(date, location) |>
        write_csv("data/weekly-covid-us.csv")
}

save_weekly_rsv <- function(hd) {
    # Read locations and participation data
    partic_rates <- read_csv("data/rsv-participation-rates.csv") |> # Has participation rates 
        select(percentage_participation, location_name)
    
    locations <- read_csv("https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/locations.csv", col_select=2:4)
    locations <- inner_join(locations, partic_rates, by="location_name")
    
    # Process and clean health data
    hd |> 
        filter(race == "All", sex == "All", age_category == "All") |> 
        filter(state != "RSV-NET") |> 
        select(-c('race', 'sex', 'age_category', 'cumulative_rate')) |> 
        rename(date=week_ending_date, location_name=state, weekly_rate=rate) |> 
        left_join(locations, by="location_name") |> 
        mutate(
            date=as.Date(date, format="%Y-%m-%d"),
            weekly_rate=parse_number(weekly_rate),
            pop_served=(percentage_participation/100) * population,
            count=round((pop_served*weekly_rate)/100000), 
            epiweek=epiweek(date),
            epiyear=epiyear(date)
        ) |> 
        select(-c(location, percentage_participation), location=location_name) |> 
        relocate(location, date, epiweek, season, epiyear, population, pop_served, count, weekly_rate) |> 
        arrange(date, location) |> # can't hurt to just make sure everything is ordered this way
        write_csv("data/weekly-rsv-us.csv")
}

#read data from healthdata.gov (may take a minute)
hd_flu_covid <- RSocrata::read.socrata(url="https://healthdata.gov/resource/g62h-syeh.json") |> 
    as_tibble() |> 
    # add some basic date information
    mutate(
        date=as.Date(date), 
        epiweek=epiweek(date), 
        epiyear=epiyear(date)
    ) |> 
    # filter out before mid Oct. 2020 because reporting was not mandatory yet:
    filter(date >= "2020-10-20")

hd_rsv <- RSocrata::read.socrata("https://data.cdc.gov/resource/29hc-w46k.json") |> 
    as_tibble()

save_weekly_flu(hd_flu_covid)
save_weekly_covid(hd_flu_covid)
save_weekly_rsv(hd_rsv)
    
# hd_counts |> 
#     group_by(date) |> 
#     summarize(na=sum(is.na(count))) |> 
#     ggplot(aes(date, na)) +
#     geom_point() +
#     scale_x_date(date_breaks="2 month", date_labels="%b-%y", guide=guide_axis(angle=45))
