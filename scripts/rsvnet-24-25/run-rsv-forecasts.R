# --------------------------------------------------------------------------------
# run-rsv-forecasts.R-------------------------------------------------------------
# notes for running forecasts: ---------------------------------------------------
# 1) `reference_date` is always the Saturday following the submission-------------
# 2) `last_date` is the date to begin forecasting, this is usually the last date--
# of available data, but you will need to subtract weeks if dropping additional---
# weeks of data-------------------------------------------------------------------

library(tidyverse)
library(INLA)
library(lubridate)
library(cowplot)

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")

inla.setOption(inla.mode="classic")

plot_seasonal <- function(fit, forecast_date) {
    summ <- fit$summary.random$epiweek |> 
        as_tibble()
    
    ggplot(summ, aes(ID)) +
        geom_vline(xintercept=epiweek(forecast_date), linetype="dashed", col="gray70") +
        geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), alpha=0.5, col=NA) +
        geom_ribbon(aes(ymin=`0.25quant`, ymax=`0.75quant`), alpha=0.5, col=NA) +
        geom_line(aes(y=mean)) +
        labs(x="epiweek", y="seasonal effect") +
        theme_half_open() +
        theme(legend.position="none")
}

# download from GitHub each week
rsv0 <- read_csv("scripts/rsvnet-24-25/2025-01-24_rsvnet_hospitalization.csv")

locations <- read_csv("scripts/rsvnet-24-25/rsvnet-locations.csv") |> 
    select(-population) # population already provided so use theirs to be safe

# Process and clean data
covid_interval <- interval(ymd("2020-03-01"), ymd("2021-06-01"))

rsv <- rsv0 |> 
    filter(age_group == "0-130", target == "inc hosp") |> 
    select(-age_group, -target) |> 
    left_join(locations, by="location") |> 
    mutate(epiweek=epiweek(date)) |> 
    filter(date >= "2018-10-06", !(date %within% covid_interval)) |> 
    mutate(count=round(value), .keep="unused") |> 
    relocate(location, abbreviation, location_name, everything()) |> 
    arrange(date, location_name) # can't hurt to just make sure everything is ordered this way

(last_date <- max(rsv$date)) # last date of available training data
reference_date <- ymd("2025-02-01") # Saturday following submission
(weeks_ahead <- as.numeric(reference_date - last_date) / 7 - 1 + 4)

fit_df <- prep_fit_data(
    rsv, last_date, weeks_ahead=weeks_ahead, 
    ex_lam=population, loc_ids=c("location", "abbreviation", "location_name")
) |> 
    mutate(
        is_state=ifelse(location == "US", 0, 1), 
        is_nat=ifelse(location == "US", 1, 0), 
        t3=t
    )

# add upper bound of `low` for truncated Poisson
fit_df <- fit_df |> 
    group_by(location) |> 
    mutate(low=1.5*max(count, na.rm=TRUE), high=Inf) |> 
    ungroup()

model <- 'Y ~ 1 + location_name + 
f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, scale.model=TRUE) + 
f(t, is_state, model="ar1", hyper=hyper_wk) + 
f(t3, is_nat, model="ar1", hyper=hyper_wk) + 
f(t2, is_state, model="ar1", hyper=hyper_wk, group=iloc, control.group=list(model="exchangeable"))'

# notice I made the prediction uncertainty param smaller than the 0.2 in the paper,
# since they use state population as an offset (doesn't seem to have a real difference tho)
fit <- fit_inla_model(fit_df, model, pc_prior_u=c(1, 0.1), forecast_date=last_date)

pred_samp <- forecast_samples(fit_df, fit, nsamp=5000)

q_wis <- c(0.01, 0.025, seq(0.05, 0.95, by=0.05), 0.975, 0.99)

# adjust the predictions for the censured Poisson, since forecast_samples
# just uses Poisson
pred_summ <- pred_samp |> 
    rowwise() |>
    mutate(predicted=list(ifelse(predicted > low, low, predicted))) |>
    ungroup() |>
    # note that these first two cols just say to keep these two var names for convenience
    summarize_quantiles(abbreviation, location_name, q=q_wis)
    

pred_summ |> 
    pivot_wider(names_from=quantile) |> 
    ggplot(aes(date)) +
    geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), fill="skyblue", alpha=0.5, col=NA) +
    geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`), fill="blue1", alpha=0.5, col=NA) +
    geom_line(aes(y=`0.5`), col="blue4") +
    geom_vline(xintercept=ymd(reference_date), col="tomato") +
    geom_point(aes(date, y=count), filter(rsv, date >= ymd(last_date) - weeks(7)), size=0.8) +
    facet_wrap(~location_name, scales="free_y") +
    scale_x_date(date_breaks="2 weeks", date_labels="%b %d", guide=guide_axis(angle=45)) +
    labs(x=NULL, y="hospitalizations") +
    theme_half_open()

ggsave(
    paste0("scripts/rsvnet-24-25/prediction-figs/rsv-predictions-", reference_date, "-drop.pdf"), 
    width=12, height=6.5
)

# Save the forecasts in format required by the hub:
pred_summ |> 
    filter(horizon > weeks_ahead - 4) |> 
    mutate(horizon=horizon - (weeks_ahead - 4))
# Maya to complete this section

