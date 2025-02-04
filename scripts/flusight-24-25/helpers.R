library(tidyverse)
library(lubridate)
library(RSocrata)
library(usmap)

fetch_flu <- function(hd=NULL) {
    locations <- read_csv(file="https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/locations.csv", col_select=c(1, 2, 4))
    
    # TODO: if hd NULL, read from RSocrata
    
    ret <- hd |> 
        as_tibble() |> 
        filter(!jurisdiction %in% c("VI", "AS", "GU", "MP")) |> 
        mutate(
            date=as.Date(weekendingdate), 
            epiweek=epiweek(date), 
            epiyear=epiyear(date),
            abbreviation=str_replace(jurisdiction, "USA", "US"),
            count=parse_number(totalconfflunewadm),
            pct_reporting=parse_number(totalconfflunewadmperchosprep)
        ) |> 
        select(date, epiweek, epiyear, abbreviation, count, pct_reporting)
    
    ret <- left_join(ret, locations, by=c("abbreviation"))
    
    ret |> 
        mutate(
            pop_weighted=population*pct_reporting,
            weekly_rate=(count*100000)/population,
            rate_weighted=(count*100000)/pop_weighted
        ) |> 
        arrange(date, location)
}

load_us_graph <- function(disease_df) {
    us0 <- us_map(regions="states")
    
    us <- us0 |> 
        filter(abbr %in% unique(disease_df$abbreviation)) |> 
        select(state=abbr)
    
    state_order <- fct_inorder(unique(disease_df$abbreviation))
    
    # to be safe, sort order of states to ensure they match their order of appearance in data
    us |> 
        mutate(state=factor(state, levels(state_order))) |> 
        arrange(state)
}

insert_iso_loc <- function(adj_mat, idx) {
    ret <- matrix(0, nrow=nrow(adj_mat)+1, ncol=ncol(adj_mat)+1)
    ret[-idx, -idx] <- adj_mat
    ret
}

# dist_xmas <- function(dates) {
#     years <- unique(year(dates))
#     xmass <- ymd(str_c(years, "-12-25"))
#     diffs <- map_dbl(dates, \(d) {
#         min(abs(d - xmass))
#     })
#     case_when(
#         diffs < 7 ~ 3,
#         diffs < 14 ~ 2,
#         diffs < 21 ~ 1,
#         TRUE ~ 0
#     )
# }

dist_xmas <- function(dates) {
    years <- unique(year(dates))
    xmass <- ymd(str_c(years, "-12-25"))
    diffs <- map_dbl(dates, \(d) {
        wh <- which.min(abs(d - xmass))
        (d - xmass)[wh]
    })
    case_when(
        abs(diffs) >= 21 ~ "none",
        diffs > -21 & diffs <= -14 ~ "2 weeks before",
        diffs < 21 & diffs >= 14 ~ "2 weeks after",
        diffs > -14 & diffs <= -7 ~ "1 week before",
        diffs < 14 & diffs >= 7 ~ "1 week after",
        TRUE ~ "week of"
    )
}

is_holidays <- function(dist_xmas) {
    ifelse(dist_xmas == "none", 0, 1)
}

is_early_covid <- function(dates) {
    ifelse(dates <= "2022-09-01", 1, 0)
}
