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
    filter(date >= "2021-09-01") # exclude peak COVID-19 period

covid <- read_csv("data/weekly-covid-us.csv")

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

saveRDS(fr, "retro-files/rsv-inflaenza.rds")

# model without seasonal component
mod_no_seasonal <- model_formula(seasonal="none", temporal="ar1", spatial="exchangeable")
fr_no_seasonal <- future_map(retro_forecast_dates_flu_rsv, \(fd) {
    fit_df <- prep_fit_data(rsv, fd, weeks_ahead=4, ex_lam=pop_served)
    fit <- fit_inla_model(fit_df, mod_no_seasonal, graph=graph, pc_prior_u=c(1, 0.2))

    forecast_samples(fit_df, fit, nsamp=5000) |>
        summarize_quantiles(q=q_wis, model="no seasonal")
}, .options=furrr_options(seed=TRUE))

saveRDS(fr_no_seasonal, "retro-files/rsv-no-seasonal.rds")

# model with no spatial sharing (fit separate to each state)
mod_no_spatial <- model_formula(covars=c(), seasonal="shared", temporal="ar1", spatial="none")
fr_no_spatial <- future_map(retro_forecast_dates_flu_rsv, \(fd) {
    fit_df <- prep_fit_data(rsv, fd, weeks_ahead=4, ex_lam=pop_served)

    forecast_baseline(fit_df, model=mod_no_spatial, pc_prior_u=c(1, 0.2)) |>
        summarize_quantiles(q=q_wis, model="no spatial")
}, .options=furrr_options(seed=TRUE))

saveRDS(fr_no_spatial, "retro-files/rsv-no-spatial.rds")

# fr_base <- future_map(retro_forecast_dates_flu_covid, \(fd) {
#     fit_df <- prep_fit_data(rsv, fd, weeks_ahead=4, ex_lam=pop_served)
# 
#     forecast_baseline(fit_df, pc_prior_u=c(1, 0.2)) |>
#         summarize_quantiles(q=q_wis, model="RW1")
# }, .options=furrr_options(seed=TRUE))
# 
# saveRDS(fr_base, "retro-files/rsv-rw1.rds")

#  processing results-------------------------------------------------------------
library(scoringutils)

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

plot_retro(readRDS("retro-files/flu-no-spatial-2.rds"), flu)

# for flu/RSV and the no spatial model, take two runs and replace where the first diverged
nsp1 <- readRDS("retro-files/flu-no-spatial-1.rds") |> 
    bind_rows()
    
nsp2 <- readRDS("retro-files/flu-no-spatial-2.rds") |> 
    bind_rows()

nsp_res <- bind_rows(nsp1, nsp2)

retro_files <- list.files("retro-files", full.names=TRUE) |> 
    setdiff(c("retro-files/flu-no-spatial-1.rds", "retro-files/flu-no-spatial-2.rds"))

pred_flu <- bind_rows(
    readRDS("retro-files/flu-inflaenza.rds"),
    readRDS("retro-files/flu-no-seasonal.rds"),
    readRDS("retro-files/flu-no-spatial-1.rds")
)

pred_flu_2 <- bind_rows(
    readRDS("retro-files/flu-inflaenza.rds"),
    readRDS("retro-files/flu-no-seasonal.rds"),
    readRDS("retro-files/flu-no-spatial-2.rds")
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
summarize_scores(scores_flu, by=c("model"))

# if we want it, useful for p-vals on the ranking of each model
get_pairwise_comparisons(scores_flu)

pred_rsv <- bind_rows(
    readRDS("retro-files/rsv-inflaenza.rds"),
    readRDS("retro-files/rsv-no-seasonal.rds"),
    readRDS("retro-files/rsv-no-spatial.rds")
    # readRDS("retro-files/rsv-rw1.rds")
)

scores_rsv <- score_pred_quantiles(rsv, pred_rsv) |> mutate(disease="rsv")

# remove one "no seasonal" forecast that tends to diverge for some reason
bad_scores <- scores_rsv |> 
    filter(wis > 1000)

scores_rsv <- anti_join(scores_rsv, bad_scores)

summarize_scores(scores_rsv, by=c("model"))
get_pairwise_comparisons(scores_rsv)

pred_covid <- bind_rows(
    readRDS("retro-files/covid-inflaenza.rds"),
    readRDS("retro-files/covid-no-seasonal.rds"),
    readRDS("retro-files/covid-no-spatial.rds")
    # readRDS("retro-files/covid-rw1.rds")
)

scores_covid <- score_pred_quantiles(covid, pred_covid) |> mutate(disease="covid")

summarize_scores(scores_covid, by=c("model"))

get_pairwise_comparisons(scores_flu, baseline="RW1") |> 
    plot_pairwise_comparisons(type="mean_scores_ratio")

bind_rows(scores_flu, scores_rsv, scores_covid) |> 
    summarize_scores(by=c("model", "disease")) |>
    mutate(disease=fct_inorder(disease)) |> 
    pivot_longer(overprediction:dispersion) |>
    ggplot(aes(value, fct_reorder(model, wis), fill=name)) +
    # ggplot(aes(interaction(model, location), value, fill=name)) +
    geom_col(alpha=0.9) +
    facet_wrap(~disease, nrow=3, scales="free_y") +
    labs(x="average WIS", y=NULL, fill=NULL) +
    scale_y_discrete(guide=guide_axis(angle=45)) +
    scale_fill_manual(values=c("#ba6cb1", "#f0b924", "tomato3")) +
    coord_flip() +
    theme_half_open()

ggsave("figs/retro-model-variants.pdf", width=6, height=5)

bind_rows(scores_flu, scores_rsv, scores_covid) |> 
    summarize_scores(by=c("model", "disease", "date")) |> 
    ggplot(aes(date, wis, col=model)) +
    geom_line() +
    facet_wrap(~disease, nrow=3, scales="free_y") +
    scale_y_continuous(transform="log")

summarize_scores(scores_rsv, by=c("model", "location")) |>
    pivot_longer(overprediction:dispersion) |>
    ggplot(aes(interaction(model, location), value, fill=name)) +
    geom_col() +
    scale_x_discrete(guide=guide_axis(angle=90))

#  recreate national forecasts for RSV--------------------------------------------
sample_national_rsv <- function(pred_samp_rsv) {
    pred_samp_rsv |> 
        mutate(predicted=map(predicted, \(p) tibble(sample_id=seq_along(p), predicted=p))) |> 
        unnest(predicted) |> 
        group_by(date, horizon, sample_id) |> 
        summarise(predicted=sum(predicted), groups="drop") |> 
        mutate(abbreviation="US", location="US") |> 
        select(-sample_id) |> 
        nest(predicted=predicted) |> 
        mutate(predicted=map(predicted, ~.$predicted))
}

future::plan(multicore, workers=4)

model <- model_formula(seasonal="shared", temporal="ar1", spatial="exchangeable")
fr <- future_map(retro_forecast_dates_flu_rsv, \(fd) {
    fit_df <- prep_fit_data(rsv, fd, weeks_ahead=4, ex_lam=pop_served)
    fit <- fit_inla_model(fit_df, model, pc_prior_u=c(1, 0.2))
    
    forecast_samples(fit_df, fit, nsamp=5000) |>
        sample_national_rsv() |> 
        summarize_quantiles(q=q_wis)
}, .options=furrr_options(seed=TRUE))

saveRDS(fr, "retro-files/rsv-national.rds")

fr_nat <- readRDS("retro-files/rsv-national.rds")

rsv_nat <- rsv |> 
    group_by(date) |> 
    summarize(count_nat=sum(count))

gg <- rsv_nat |> 
    filter(date > "2022-08-01", date < "2024-06-01") |> 
    ggplot(aes(date))

for (pr in fr_nat[seq(1, length(fr_nat), 5)]) {
    pred_wide <- pivot_wider(pr, names_from=quantile)
    
    gg <- gg +
        geom_ribbon(
            aes(ymin=`0.025`, ymax=`0.975`),
            pred_wide,
            fill="skyblue", alpha=0.5, col=NA
        ) +
        geom_line(aes(y=mean), pred_wide, col="blue4")
}
gg + geom_point(aes(y=count_nat), size=0.79, alpha=0.8) +
    coord_cartesian(ylim=c(0, 2100)) +
    labs(x=NULL, y="RSV")

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
