library(tidyverse)
library(INLA)
library(lubridate)
library(furrr)

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")

inla.setOption(inla.mode="classic")

flu <- read_csv("data/weekly-flu-us.csv") |> 
    filter(date >= "2021-09-01") # exclude peak COVID-19 period

graph <- load_us_graph(flu) |> sf2mat()

retro_forecast_dates <- c(
    seq.Date(ymd("2022-09-03"), ymd("2023-04-29"), by="3 weeks"),
    seq.Date(ymd("2023-09-02"), ymd("2024-04-27"), by="3 weeks")
)

q_wis <- c(0.01, 0.025, seq(0.05, 0.95, by=0.05), 0.975, 0.99)

### work to send to separate R process
future::plan(multicore, workers=4)

# TODO: you could have a save fit option that saves each fit_df and fit within map into a folder
# with name matching the main rds

model <- model_formula(seasonal="shared", temporal="ar1", spatial="besagproper")
fr <- future_map(retro_forecast_dates, \(fd) {
    fit_df <- prep_fit_data(flu, fd, weeks_ahead=4, ex_lam=population)
    fit <- fit_inla_model(fit_df, model, graph=graph)

    forecast_samples(fit_df, fit, nsamp=5000) |>
        summarize_quantiles(q=q_wis, model="INFLAenza")
}, .options=furrr_options(seed=TRUE))

saveRDS(fr, "retro-files/inflaenza-11-7.rds")

mod_no_seasonal <- model_formula(seasonal="none", temporal="ar1", spatial="besagproper")
fr_no_seasonal <- future_map(retro_forecast_dates, \(fd) {
    fit_df <- prep_fit_data(flu, fd, weeks_ahead=4, ex_lam=population)
    fit <- fit_inla_model(fit_df, mod_no_seasonal, graph=graph)

    forecast_samples(fit_df, fit, nsamp=5000) |>
        summarize_quantiles(q=q_wis, model="no seasonal")
}, .options=furrr_options(seed=TRUE))

saveRDS(fr_no_seasonal, "retro-files/no-seasonal-11-7.rds")

mod_no_spatial <- model_formula(covars=c(), seasonal="shared", temporal="ar1", spatial="none")
fr_no_spatial <- future_map(retro_forecast_dates, \(fd) {
    fit_df <- prep_fit_data(flu, fd, weeks_ahead=4, ex_lam=population)

    forecast_baseline(fit_df, model=mod_no_spatial) |>
        summarize_quantiles(q=q_wis, model="no spatial")
}, .options=furrr_options(seed=TRUE))

saveRDS(fr_no_spatial, "retro-files/no-spatial-11-8.rds")

fr_base <- future_map(retro_forecast_dates, \(fd) {
    fit_df <- prep_fit_data(flu, fd, weeks_ahead=4, ex_lam=population)

    forecast_baseline(fit_df) |>
        summarize_quantiles(q=q_wis, model="RW1")
}, .options=furrr_options(seed=TRUE))

saveRDS(fr_base, "retro-files/rw1-11-7.rds")

#  processing results-------------------------------------------------------------
library(scoringutils)

source("src/scoring.R")

plot_retro <- function(fr, disease_df, loc_sub=NULL) {
    date_min <- min(map_vec(fr, ~min(.$date))) - weeks(4)
    date_max <- max(map_vec(fr, ~max(.$date))) + weeks(4)
    
    if (is.null(loc_sub)) {
        loc_sub <- unique(fr[[1]]$location)
    }
    else {
        loc_sub <- intersect(unique(fr[[1]]$location), loc_sub)
    }
    
    df_sub <- filter(
        disease_df,
        date >= date_min,
        date <= date_max,
        location %in% loc_sub
    )
    
    gg <- ggplot(df_sub, aes(date)) +
        facet_wrap(~location, scales="free_y")
    # coord_cartesian_panels(panel_limits=panel_limits) +
    # scale_x_date(breaks=recent_dates[seq(1, length(recent_dates), 2)], date_labels="%d %b", guide=guide_axis(angle=45)) +
    
    for (pr in fr) {
        pred_sub <- filter(pr, location %in% loc_sub) |>
            pivot_wider(names_from=quantile)
        
        gg <- gg +
            geom_ribbon(
                aes(ymin=`0.025`, ymax=`0.975`),
                pred_sub,
                fill="skyblue", alpha=0.5, col=NA
            ) +
            geom_ribbon(
                aes(ymin=`0.25`, ymax=`0.75`),
                pred_sub,
                fill="blue1", alpha=0.5, col=NA
            ) +
            geom_line(aes(y=mean), pred_sub, col="blue4")
    }
    gg + geom_point(aes(y=count), size=0.79, alpha=0.8)
}

plot_retro(readRDS("retro-files/no-spatial-11-7.rds"), flu, c("California", "Oregon", "Washington", "Nevada"))

# for the no spatial model, take two runs and replace where the first diverged
nsp1 <- readRDS("retro-files/no-spatial-11-7.rds") |> 
    bind_rows() |> 
    filter(location != "California")
    
nsp2 <- readRDS("retro-files/no-spatial-11-8.rds") |> 
    bind_rows() |> 
    filter(location == "California")

nsp_res <- bind_rows(nsp1, nsp2)

all_pred <- bind_rows(
    readRDS("retro-files/inflaenza-11-7.rds"),
    readRDS("retro-files/no-seasonal-11-7.rds"),
    nsp_res,
    readRDS("retro-files/rw1-11-7.rds")
)

scores <- score_pred_quantiles(flu, all_pred) |> 
    mutate()

summarize_scores(scores) |>
    pivot_longer(overprediction:dispersion) |>
    ggplot(aes(value, fct_reorder(model, wis), fill=name)) +
    # ggplot(aes(interaction(model, location), value, fill=name)) +
    geom_col()

#' forecast_retro <- function(
#'         fit_df, model, forecast_dates, weeks_ahead=4, loc_sub=NULL, 
#'         q=c(0.025, 0.25, 0.5, 0.75, 0.975), nsamp=5000,...) {
#'     
#'     if (is.null(loc_sub))
#'         loc_sub <- unique(fit_df$location)
#'     
#'     map(forecast_dates, \(fd) {
#'         fd_exact <- min(filter(fit_df, date >= fd)$date)
#'         fd_end <- fd_exact + weeks(weeks_ahead-1)
#'         
#'         # trim the fit_df to contain the training data and just the exact dates/locations to forecast
#'         pred_df <- fit_df |> 
#'             filter(date < fd | (date <= fd_end & location %in% loc_sub))
#'         
#'         # replace the dates to forecast with NA
#'         pred_idx <- which(pred_df$date >= fd_exact)
#'         pred_df$count[pred_idx] <- NA
#'         
#'         fit <- fit_inla_model(
#'             pred_df, fd, model, 
#'             pred_idx=pred_idx, joint_forecast=FALSE, ...
#'         )
#'         
#'         pred_summ <- sample_count_predictions(pred_df, fit, fd_exact, nsamp=nsamp) |> 
#'             summarize_quantiles(q=q) |> 
#'             pivot_wider(names_from=quantile)
#'         
#'         list(pred=pred_summ, fit_df=pred_df, fit=fit)
#'     })
#' }
#' 
#' update_pred_quantiles <- function(fr, q, nsamp=5000) {
#'     map(fr, \(el) {
#'         fd <- min(el$pred$date)
#'         
#'         pred_summ <- sample_count_predictions(el$fit_df, el$fit, fd, nsamp=nsamp) |> 
#'             summarize_quantiles(q=q) |> 
#'             pivot_wider(names_from=quantile)
#'         
#'         list(pred=pred_summ, fit_df=el$fit_df, fit=el$fit)
#'     })
#' }
#' 

#' 
#' 
#' #' coverage_timeseries
#' #' @description
#' #' Produce a continuous sequence of forecast quantiles within an interval for a specific forecast
#' #' horizon. Useful for coverage plots.
#' #'
#' #' @param fit_df a dataframe of data as produced by `prep_fit_data_<disease>`
#' #' @param model an INLA model formula
#' #' @param start start date to begin forecasting. This does not have to be an exact date
#' #' within `fit_df` - truncation of `fit_df` begins for dates >= `start`
#' #' @param end end date to stop forecasting. This does not have to be an exact date
#' #' within `fit_df`
#' #' @param horizon number of weeks ahead to do each forecast
#' #' @param loc_sub optionally, limit forecasts to just a vector of locations. Note that model
#' #' fitting/training still uses all locations in this case
#' #' @param nsamp number of samples to draw with `sample_count_predictions`
#' #' @param ... additional arguments to passed to `fit_inla_model`
#' #'
#' #' @return a tibble containing forecast quantiles for each date between `start` and `end`
#' coverage_timeseries <- function(fit_df, model, start, end, horizon=2, loc_sub=NULL, nsamp=2000,...) {
#'     forecast_dates <- unique(filter(fit_df, date >= start, date <= end)$date)
#'     if (is.null(loc_sub))
#'         loc_sub <- unique(fit_df$location)
#'     
#'     map_dfr(forecast_dates, \(fd) {
#'         obs_end_date <- fd - weeks(horizon)
#'         
#'         # trim the fit_df to contain the training data and just the exact date to forecast
#'         pred_df <- fit_df |> 
#'             filter(date <= obs_end_date | (date == fd & location %in% loc_sub))
#'         
#'         # replace the date to forecast with NA
#'         pred_idx <- which(pred_df$date == fd)
#'         pred_df$count[pred_idx] <- NA
#'         
#'         fit <- fit_inla_model(
#'             pred_df, fd, model, 
#'             pred_idx=pred_idx, joint_forecast=FALSE, ...
#'         )
#'         
#'         # TODO: change to a list including fit object
#'         sample_count_predictions(pred_df, fit, fd, nsamp=nsamp) |> 
#'             summarize_quantiles() |> 
#'             pivot_wider(names_from=quantile) |> 
#'             mutate(last_date=obs_end_date, .after=date)
#'     })
#' }
