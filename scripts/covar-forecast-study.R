# preliminary work on showcasing how different covariance models for the----------
# short-term temporal component have a big effect on forecasts--------------------


library(tidyverse)
library(INLA)
library(sf)
library(spdep)
library(lubridate)
library(cowplot)

source("fetch-data/read-flu-data.R")
source("model-prep-and-fit.R")
source("model-summarize-and-plot.R")

fit_forecast_rep <- function(fit_df, pred_date, model) {
    hyper_wk <- list(prec=list(prior="pc.prec", param=c(1, 0.01)))
    
    pred_idx <- which(fit_df$date >= pred_date)
    mod <- as.formula(model)
    
    inla(
        mod, family="poisson", data=fit_df,
        E=fit_df$ex_lam,
        quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975),
        selection=list(Predictor=pred_idx),
        control.compute=list(return.marginals.predictor=TRUE),
        control.predictor=list(link=1)
    )
}

forecast_rep_plots <- function(training_weeks, target_date, model) {
    pred_date <- target_date - duration(min(training_weeks, 8), "weeks") # date to start making count predictions
    
    fit_df <- flu0 |> 
        select(-location, location=location_name, count=value) |> 
        filter(
            date < target_date, 
            date >= (target_date - duration(training_weeks, "weeks")), 
            location == "Georgia"
        ) |> 
        prep_fit_data(weeks_ahead=4)
    
    fit <- fit_forecast_rep(fit_df, pred_date, model)
    
    list(
        seasonal=plot_seasonal(fit, target_date),
        weekly=plot_latent_pred(fit_df, fit, target_date),
        plot_pred=plot_count_pred(fit_df, fit, target_date, pred_date)
    )
}

plot_seasonal <- function(fit, target_date) {
    fit$summary.random$epiweek |>
        as_tibble() |> 
        ggplot(aes(ID, mean)) +
        geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col=NA, alpha=0.3) +
        geom_line(col="gray30") +
        geom_vline(xintercept=epiweek(target_date-duration(1, "week")), col="plum2", linetype="31", linewidth=1.2, alpha=0.7) +
        # facet_wrap(~state) +
        labs(x="epiweek", title="seasonal effect", y=NULL) +
        theme_half_open()
}

plot_latent_pred <- function(fit_df, fit, target_date) {
    fit_df |> 
        mutate(
            mean=fit$summary.random$t$mean,
            ymin=fit$summary.random$t$`0.025quant`,
            ymax=fit$summary.random$t$`0.975quant`
        ) |> 
        # filter(date >= "2023-12-01") |> 
        ggplot(aes(date)) +
        geom_ribbon(aes(ymin=ymin, ymax=ymax), col=NA, alpha=0.3) +
        geom_line(aes(y=mean), col="gray30") +
        geom_vline(xintercept=target_date-duration(1, "week"), col="plum2", linetype="31", linewidth=1.2, alpha=0.7) +
        labs(y=NULL, x="week", title="weekly effect") +
        scale_x_date(guide=guide_axis(angle=30)) +
        theme_half_open() +
        theme(legend.position="none")
}

plot_count_pred <- function(fit_df, fit, target_date, pred_date) {
    pred_samples <- sample_count_predictions(fit_df, fit, pred_date, nsamp=2000)
    
    pred_summ <- summarize_quantiles(pred_samples, NULL, q=c(0.025, 0.5, 0.975)) |> 
        pivot_wider(names_from=quantile, values_from=value)
    
    fit_df |> 
        filter(date >= pred_date) |>
        ggplot(aes(date)) +
        geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), pred_summ, col=NA, alpha=0.3) +
        geom_line(aes(y=mean), pred_summ, col="gray30") +
        geom_point(aes(y=count), size=1.1, col="black") +
        geom_vline(xintercept=target_date-duration(1, "week"), col="plum2", linetype="31", linewidth=1.2, alpha=0.7) +
        # scale_color_manual(values=c("#0c2c84", "#1d91c0", "#7fcdbb")) +
        # scale_fill_manual(values=c("#0c2c84", "#1d91c0", "#7fcdbb")) +
        scale_x_date(date_labels="%Y-%m", guide=guide_axis(angle=30)) +
        labs(x="Week", y="Hospitalizations") +
        theme_half_open() +
        theme(legend.position="none")
}

# Data prep and model fitting-----------------------------------------------------
# Read in current data
flu0 <- fetch_flu()
max(flu0$date)
training_weeks <- 150
target_date <- ymd("2024-01-13") # first "future" date to start forecasting
# TODO: using a seasonal effect for IDIG, but consider redoing without
form <- 'count ~ f(epiweek, model="rw2", cyclic=TRUE, scale.model=TRUE) + 
    f(t, model="rw1", cyclic=FALSE, scale.model=TRUE)'

fplots <- forecast_rep_plots(training_weeks, target_date, form)


fplots$seasonal
