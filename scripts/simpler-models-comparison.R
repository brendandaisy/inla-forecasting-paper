# --------------------------------------------------------------------------------
# Table S1: comparing the proposed model to other variations removing or----------
# adding seasonal or spatial components-------------------------------------------
# --------------------------------------------------------------------------------
# usage: each code block is set up for a single disease, so needs to be edited----
# slightly to run for other diseases----------------------------------------------
# --------------------------------------------------------------------------------

library(tidyverse)
library(INLA)
library(lubridate)
library(furrr)

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")
source("src/scoring.R")

inla.setOption(inla.mode="classic")

flu <- read_csv("data/weekly-flu-us.csv") |> 
    filter(date >= "2021-09-01", location != "US") # exclude peak COVID-19 period

covid <- read_csv("data/weekly-covid-us.csv") |> 
    filter(location != "US")

covid_interval <- interval(ymd("2020-03-01"), ymd("2021-06-01"))
rsv <- read_csv("data/weekly-rsv-us.csv") |> 
    filter(!(date %within% covid_interval))

graph <- load_us_graph(flu) |> sf2mat()

retro_forecast_dates_flu_rsv <- c(
    seq.Date(ymd("2022-09-10"), ymd("2023-04-20"), by="1 weeks"),
    seq.Date(ymd("2023-09-09"), ymd("2024-04-20"), by="1 weeks")
)

# include summer months for COVID
retro_forecast_dates_covid <- seq.Date(ymd("2022-09-10"), ymd("2024-04-20"), by="1 weeks")

q_wis <- c(0.01, 0.025, seq(0.05, 0.95, by=0.05), 0.975, 0.99)

#  work to send to separate R process---------------------------------------------
future::plan(multicore, workers=4)

# full model
model <- model_formula(seasonal="shared", temporal="ar1", spatial="exchangeable")
fr <- future_map(retro_forecast_dates_flu_rsv, \(fd) {
    fit_df <- prep_fit_data(rsv, fd, weeks_ahead=4, ex_lam=pop_served)
    fit <- fit_inla_model(fit_df, model, graph=graph, pc_prior_u=c(1, 0.2))

    forecast_samples(fit_df, fit, nsamp=5000) |>
        summarize_quantiles(q=q_wis, model="INFLAenza")
}, .options=furrr_options(seed=TRUE))

saveRDS(fr, "results/simpler-mod-comp/rsv-inflaenza.rds")

# model without seasonal component
mod_no_seasonal <- model_formula(seasonal="none", temporal="ar1", spatial="exchangeable")
fr_no_seasonal <- future_map(retro_forecast_dates_flu_rsv, \(fd) {
    fit_df <- prep_fit_data(rsv, fd, weeks_ahead=4, ex_lam=pop_served)
    fit <- fit_inla_model(fit_df, mod_no_seasonal, graph=graph, pc_prior_u=c(1, 0.2))

    forecast_samples(fit_df, fit, nsamp=5000) |>
        summarize_quantiles(q=q_wis, model="no seasonal")
}, .options=furrr_options(seed=TRUE))

saveRDS(fr_no_seasonal, "results/simpler-mod-comp/rsv-no-seasonal.rds")

# model with no spatial sharing (fit separate to each state)
mod_no_spatial <- model_formula(covars=c(), seasonal="shared", temporal="ar1", spatial="none")
fr_no_spatial <- future_map(retro_forecast_dates_flu_rsv, \(fd) {
    fit_df <- prep_fit_data(rsv, fd, weeks_ahead=4, ex_lam=pop_served)

    forecast_baseline(fit_df, model=mod_no_spatial, pc_prior_u=c(1, 0.2)) |>
        summarize_quantiles(q=q_wis, model="no spatial")
}, .options=furrr_options(seed=TRUE))

saveRDS(fr_no_spatial, "results/simpler-mod-comp/rsv-no-spatial.rds")

# model with independent seasonal effects for each state
mod_iid_seasonal <- model_formula(seasonal="iid", temporal="ar1", spatial="besagproper")
fr_iid_seasonal <- future_map(retro_forecast_dates_covid, \(fd) {
    fit_df <- prep_fit_data(covid, fd, weeks_ahead=4, ex_lam=population)
    fit <- fit_inla_model(fit_df, mod_iid_seasonal, graph=graph, pc_prior_u=c(1, 1))
    
    forecast_samples(fit_df, fit, nsamp=5000) |>
        summarize_quantiles(q=q_wis, model="iid seasonal")
}, .options=furrr_options(seed=TRUE))

saveRDS(fr_iid_seasonal, "results/simpler-mod-comp/covid-iid-seasonal.rds")

# model with geographic seasonality based on HHS regions and outlier states
hhs_lookup <- load_us_graph(flu) |>
    st_drop_geometry() |> 
    mutate(region=ifelse(region == 9, 3, region)) |> 
    mutate(region=factor(region, labels=c("Northeast", "Midwest", "South", "West/Pacific"))) |> 
    select(-division) |> 
    nest(locs=c(state))

mod_geo_seasonal <- "count ~ 1 + location + 
f(epiweek, model='rw2', scale.model=TRUE, cyclic=TRUE, hyper=hyper_epwk, group=season_group, control.group=list(model='iid')) + 
f(t, model='ar1', hyper=hyper_wk) + f(t2, model='ar1', hyper=hyper_wk, group=iloc, control.group=list(model='exchangeable'))"
fr_geo_seasonal <- future_map(retro_forecast_dates_flu_rsv, \(fd) {
    fit_df <- prep_fit_data(rsv, fd, weeks_ahead=4, ex_lam=pop_served) |> 
        mutate(season_group=case_when(
            location == "Puerto Rico" ~ 1,
            location == "Hawaii" ~ 2,
            location == "Alaska" ~ 3,
            location %in% hhs_lookup$locs[[1]]$state ~ 4,
            location %in% hhs_lookup$locs[[2]]$state ~ 5,
            location %in% hhs_lookup$locs[[3]]$state ~ 6,
            location %in% hhs_lookup$locs[[4]]$state ~ 7
        ))
    
    fit <- fit_inla_model(fit_df, mod_geo_seasonal, pc_prior_u=c(1, 0.2))
    
    forecast_samples(fit_df, fit, nsamp=5000) |>
        summarize_quantiles(q=q_wis, model="geo seasonal")
}, .options=furrr_options(seed=TRUE))

saveRDS(fr_geo_seasonal, "results/simpler-mod-comp/rsv-geo-seasonal.rds")

# processing results-------------------------------------------------------------
library(scoringutils)

# note that for flu/RSV and the no spatial model, some forecast tended to diverge, so 
# take two independent runs and choose the one that fit successfully
pred_flu <- bind_rows(
    readRDS("results/simpler-mod-comp/flu-inflaenza.rds"),
    readRDS("results/simpler-mod-comp/flu-no-seasonal.rds"),
    readRDS("results/simpler-mod-comp/flu-no-spatial-1.rds"), 
    readRDS("results/simpler-mod-comp/flu-iid-seasonal.rds"),
    readRDS("results/simpler-mod-comp/flu-geo-seasonal.rds")
)

pred_flu_2 <- bind_rows(
    readRDS("results/simpler-mod-comp/flu-inflaenza.rds"),
    readRDS("results/simpler-mod-comp/flu-no-seasonal.rds"),
    readRDS("results/simpler-mod-comp/flu-no-spatial-2.rds"),
    readRDS("results/simpler-mod-comp/flu-iid-seasonal.rds"),
    readRDS("results/simpler-mod-comp/flu-geo-seasonal.rds")
)

scores_flu <- score_pred_quantiles(flu, pred_flu) |> mutate(disease="flu")
scores_flu_2 <- score_pred_quantiles(flu, pred_flu_2) |> mutate(disease="flu")

bad_scores <- scores_flu |> 
    slice_max(wis, n=54)

fixed_scores <- bad_scores |> 
    left_join(scores_flu_2, by=c("date", "location", "horizon", "model")) |> 
    mutate(wis=pmin(wis.x, wis.y)) |> 
    select(date:model, wis, pic_50=pic_50.x, pic_95=pic_95.x) |> 
    filter(wis < 1000) # there were still 9 scores that were way higher than the fixed ones, so just remove from the average

scores_flu <- scores_flu |> 
    anti_join(bad_scores) |> 
    bind_rows(fixed_scores)

# for the table
summarize_scores(scores_flu, by=c("model")) |> 
    mutate(pct_better=(wis-first(wis))/wis)

pred_rsv <- bind_rows(
    readRDS("results/simpler-mod-comp/rsv-inflaenza.rds"),
    readRDS("results/simpler-mod-comp/rsv-no-seasonal.rds"),
    readRDS("results/simpler-mod-comp/rsv-no-spatial.rds"),
    readRDS("results/simpler-mod-comp/rsv-iid-seasonal.rds"),
    readRDS("results/simpler-mod-comp/rsv-geo-seasonal.rds")
)

scores_rsv <- score_pred_quantiles(rsv, pred_rsv) |> mutate(disease="rsv")

# remove one "no seasonal" forecast that tends to diverge for some reason
bad_scores <- scores_rsv |> 
    filter(wis > 1000)

scores_rsv <- anti_join(scores_rsv, bad_scores)

summarize_scores(scores_rsv, by=c("model")) |> 
    mutate(pct_better=(wis-first(wis))/wis)

pred_covid <- bind_rows(
    readRDS("results/simpler-mod-comp/covid-inflaenza.rds"),
    readRDS("results/simpler-mod-comp/covid-no-seasonal.rds"),
    readRDS("results/simpler-mod-comp/covid-no-spatial.rds"),
    readRDS("results/simpler-mod-comp/covid-iid-seasonal.rds"),
    readRDS("results/simpler-mod-comp/covid-geo-seasonal.rds")
)

scores_covid <- score_pred_quantiles(covid, pred_covid) |> mutate(disease="covid")

summarize_scores(scores_covid, by=c("model")) |> 
    mutate(pct_better=(wis-first(wis))/wis)
