library(tidyverse)
library(sf)
library(spdep)

# TODO: this isn't really general enough to go in the package
load_us_graph <- function(disease_df, f="data/us-state-boundaries.shp") {
    us0 <- read_sf(f)

    us <- us0 |>
        filter(name %in% unique(disease_df$location)) |>
        select(state=name, division, region)

    state_order <- fct_inorder(unique(disease_df$location))

    # to be safe, sort order of states to ensure they match their order of appearance in data
    us |>
        mutate(state=factor(state, levels(state_order))) |>
        arrange(state)
}

sf2mat <- function(sf_graph) {
    sf_graph |> 
        poly2nb() |> 
        nb2mat(style="W", zero.policy=TRUE)
}

# this now can work just for adding the necessary indices without adding NAs, with weeks_ahead == 0
prep_fit_data_rsv <- function(rsv, weeks_ahead=4, ex_lam=pop_served) {
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
                date=duration(1:weeks_ahead, "week") + max(ret$date),
                t=1:weeks_ahead + max(ret$t),
                epiweek=epiweek(date)
            ),
            distinct(ret, location, iloc)
        )
        
        # go and find most recent population/ex_lam values for each state, in case they change
        # note this assumes the prediction df follows most recent ex_lam
        location_data <- ret |> 
            slice_max(date) |> 
            distinct(location, ex_lam)
        
        # add on most recent pop/ex_lam values, with protection against adding/deleting unintended rows
        pred_df <- pred_df |> 
            left_join(
                location_data, by=c("location"), unmatched="error", relationship="many-to-one"
            )
        
        ret <- bind_rows(ret, pred_df) # add to data for counts to be NAs
    }
    
    ret |> 
        mutate(t2=t) |> 
        arrange(t, iloc)
}

# loc_ids which column(s) represent locations. Can include multiple one-to-one cols, if you want e.g.
# an abbreviated and full state name in the results
prep_fit_data <- function(
        disease_df, forecast_date=NULL, weeks_ahead=4, ex_lam=population,
        loc_ids=c("location")
) {
    ret <- disease_df |> 
        mutate(
            iloc=as.numeric(fct_inorder(.data[[loc_ids[1]]])),
            ex_lam={{ex_lam}}
        )
    
    if (!is.null(forecast_date)) {
        ret <- filter(ret, date <= forecast_date)
        
        pred_df <- expand_grid( # makes pairs of new weeks X location
            tibble(
                date=max(ret$date) + weeks(1:weeks_ahead),
                # t=1:weeks_ahead + max(ret$t),
                epiweek=epiweek(date)
            ),
            distinct(ret, across(all_of(loc_ids)), iloc)
        )
        
        location_data <- ret |> 
            slice_max(date) |> 
            distinct(across(all_of(loc_ids)), ex_lam)
        
        pred_df <- left_join(
            pred_df, location_data, 
            by=loc_ids, unmatched="error", relationship="many-to-one"
        )
        
        ret <- bind_rows(ret, pred_df) # add to data for counts to be NAs
    }
    
    # index dates with a sequential time index, accounting for holes in the
    # data and spacing indices properly
    date_seq <- seq.Date(min(ret$date), max(ret$date), "1 week")
    
    date_ind <- tibble(t=seq_along(date_seq), date=date_seq) |> 
        filter(date %in% unique(ret$date))
    
    ret |> 
        left_join(date_ind, by=c("date"), relationship="many-to-one") |> 
        mutate(t2=t) |> 
        arrange(t, iloc)
}
