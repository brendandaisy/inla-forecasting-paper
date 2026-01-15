library(tidyverse)
library(INLA)
library(lubridate)
library(cowplot)

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")

inla.setOption(inla.mode="classic")
theme_set(theme_cowplot())

# Read-in recent metrocast data --------------------------------------------
curr_resp_season <- 2025

quantiles_needed <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
weeks_to_drop <- 0

## Pull data in from github and format properly
flu_data0 <- read_csv('https://raw.githubusercontent.com/reichlab/flu-metrocast/refs/heads/main/target-data/latest-data.csv')
location_info <- read_csv('https://raw.githubusercontent.com/reichlab/flu-metrocast/refs/heads/main/auxiliary-data/locations.csv')
states <- unique(location_info$state_abb)

(forecast_date <- max(flu_data0$target_end_date)) # final date to include in training (Saturday before horizon -1)

flu_data <- flu_data0 |> 
    mutate(epiweek = epiweek(target_end_date),
           year = epiyear(target_end_date))  |> 
    mutate(resp_season = ifelse(epiweek>=40, year, year-1)) |> 
    filter(!resp_season %in% c(2008, 2009, 2020)) |>
    group_by(location, target, resp_season) |> 
    arrange(target_end_date) |> 
    mutate(resp_season_week = seq_along(epiweek)) |> 
    ungroup() |> 
    left_join(location_info |> 
                  select(location, state, state_abb, location_name, population), by = 'location') |> 
    rename(date=target_end_date) |> 
    mutate(observation=observation/100)

flu_data <- flu_data |> 
    left_join(location_info, keep=FALSE) |> 
    filter(state != location_name)

summarize_quantiles_na_rm <- function(pred_samples, ..., nat_samps=NULL, q=c(0.025, 0.25, 0.5, 0.75, 0.975), model=NULL) {
    pred_samples |> 
        # filter(date >= forecast_date) |> 
        bind_rows(nat_samps) |> 
        unnest(predicted) |> 
        group_by(date, location, horizon, ...) |> 
        summarize(
            mean=mean(predicted, na.rm=TRUE),
            qs=list(value=quantile(predicted, probs=q, na.rm=TRUE)), 
            .groups="drop"
        ) |> 
        unnest_wider(qs) |>
        pivot_longer(contains("%"), names_to="quantile") |> 
        mutate(quantile=parse_number(quantile)/100, model=model)
}

pred_state_metro <- function(flu_df, state_abb) {
    flu_state <- filter(flu_df, state_abb == {{state_abb}})

    fit_df <- prep_fit_data(flu_state, forecast_date, ex_lam=1)
    
    if (length(unique(flu_state$location)) == 1) {
        model <- model_formula(
            response="observation", covars=c(), seasonal="shared", temporal="ar1", spatial="none"
        )
    } else {
        model <- model_formula(
            response="observation", seasonal="shared", temporal="ar1", spatial="exchangeable"
        )
    }
    
    cens <- min(fit_df$observation[fit_df$observation > 0], na.rm=TRUE)/2
    fit <- fit_inla_model(
        fit_df, model, family="beta", response=observation,
        control.family=list(beta.censor.value=cens)
    )
    
    pred_samp <- forecast_samples(fit_df, fit, nsamp=5000, family="beta", response=observation)
    summarize_quantiles_na_rm(pred_samp, q=quantiles_needed) |> 
        mutate(value=100*value)
}

# main line to fit and make predictions for each state
pred_summ_all <- map_dfr(states, \(s) pred_state_metro(flu_data, s))

reference_date <- forecast_date + weeks(1)

# plot all the county level predictions
pred_summ_all |> 
    pivot_wider(names_from=quantile) |> 
    ggplot(aes(date)) +
    geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), fill="skyblue", alpha=0.5, col=NA) +
    geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`), fill="blue1", alpha=0.5, col=NA) +
    geom_line(aes(y=`0.5`), col="blue4") +
    geom_point(aes(date, y=100*observation), filter(flu_data, date >= forecast_date - weeks(7)), size=0.8) +
    facet_wrap(~location, scales="free_y") +
    scale_x_date(date_breaks="2 weeks", date_labels="%b %d", guide=guide_axis(angle=45)) +
    labs(x=NULL, y="metro percent ed/ili visits") +
    theme_half_open()

ggsave(
    paste0("figs/metrocast-25-26/INFLAenza-prop-ed-metro", reference_date, ".pdf"), 
    width=12, height=8.5
)

# steps to produce output file for pull request
metro_output <- pred_summ_all |> 
    mutate(
        target=ifelse(location == "nyc", "ILI ED visits pct", "Flu ED visits pct"),
        reference_date=reference_date,
        horizon=horizon-1,
        target_end_date=reference_date + horizon*7,
        output_type_id=quantile,
        output_type='quantile'
    ) |> 
    # arrange(abbreviation, horizon, quantile) |>
    select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value)

# grab state level predictions from the flusight output file----------------------

# Rajath: you'll need to change this path to wherever you keep your local fork
# of the FluSight repository:
path_flusight_repo <- "output/flusight/FluSight-forecast-hub/"

# Rajath: you'll need to change this path to wherever you keep your local fork
# of the metrocast repository:
path_metrocast_repo <- "output/metrocast/flu-metrocast"

flusight_output <- read_csv(paste0(path_flusight_repo, "model-output/UGA_flucast-INFLAenza/", reference_date, "-UGA_flucast-INFLAenza.csv"))
location_info_flusight <- read_csv(file="scripts/flusight-25-26/locations.csv", col_select=c(1, 2, 4))

# remove NY since state predictions aren't needed
abb_states_needed <- setdiff(location_info$state_abb, "NY")

flusight_output_sub <- flusight_output |> 
    filter(target == "wk inc flu prop ed visits") |> 
    left_join(location_info_flusight) |> 
    filter(
        abbreviation %in% abb_states_needed,
        output_type_id %in% quantiles_needed
    ) |> 
    select(-location) |> 
    left_join(
        filter(location_info, state == location_name), 
        by=join_by(abbreviation == state_abb)
    ) |> 
    mutate(
        target="Flu ED visits pct", 
        value=100*value
    ) |> 
    select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value)

# Combine forecasts for metro areas and states------------------------------------
comb_output <- bind_rows(metro_output, flusight_output_sub)

write_csv(comb_output, paste0(path_metrocast_repo, "model-output/NAU-INFLAenza/", reference_date, "-NAU-INFLAenza.csv"))
