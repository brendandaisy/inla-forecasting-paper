library(tidyverse)
library(INLA)
library(lubridate)

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")

inla.setOption(inla.mode="classic")

covid <- read_csv("data/weekly-covid-us.csv")

# ok so lets try to work with all data going to mid 2020
ggplot(covid, aes(date, count)) +
    geom_point() +
    facet_wrap(~location, scales="free")

# TODO check that the short-term effect looks similar for AR1 in the exchangeable vs. besag versions,
# maybe with a subset of 6 or so states

# Inspecting forecasts with the base model----------------------------------------
fd_22_23 <- c("2022-02-01", "2022-11-01", "2023-02-15", "2023-06-01") # some older dates 1st
fd_23_24 <- c("2023-10-01", "2023-11-15", "2023-12-30", "2024-02-01")

# new way to specify the current "base" model
model <- model_formula("shared", temporal="ar1", spatial="besagproper")
graph <- load_us_graph(covid) |> sf2mat()

# fit the model and sample predictions for each of the timepoints
pred_quantiles <- map(fd_23_24, \(fd) {
    fit_df <- covid |> 
        filter(date <= fd) |> 
        prep_fit_data_flu_covid(weeks_ahead=4, ex_lam=population)
    
    fit <- fit_inla_model(fit_df, fd, model, pc_prior_u=c(1, 1), graph=graph)
    
    pred_summ <- sample_count_predictions(fit_df, fit, fd, nsamp=5000) |> 
        summarize_quantiles() |> 
        pivot_wider(names_from=quantile)
    
    list(pred=pred_summ, fit_df=fit_df, fit=fit)
})

# Plot the retrospective forecasts compared to truth-----------------------------
sub_states <- c(
    "Alaska", "California", "Florida", "Georgia", "Louisiana", "Iowa", 
    "Maine", "New York", "Vermont", "Virginia"
)

# TODO add US to your considerations

covid_sub <- filter(
    covid, 
    date >= (ymd(min(fd_23_24))-weeks(8)), 
    date < (ymd(max(fd_23_24))+weeks(8)),
    location %in% sub_states
)

gg <- ggplot(covid_sub, aes(date)) +
    geom_point(aes(y=count), size=0.79, alpha=0.8) +
    facet_wrap(~location, scales="free_y") +
    # coord_cartesian_panels(panel_limits=panel_limits) +
    # scale_x_date(breaks=recent_dates[seq(1, length(recent_dates), 2)], date_labels="%d %b", guide=guide_axis(angle=45)) +
    theme(title=element_text(size=12))

for (el in pred_quantiles) {
    pred_sub <- filter(el$pred, location %in% sub_states)
    
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
gg

# Other plots for model outputs---------------------------------------------------
last_forecast <- last(pred_quantiles)

seasonal <- last_forecast$fit$summary.random$epiweek

weekly_main <- last(pred_quantiles)$fit$summary.random$t |> 
    as_tibble() |> 
    mutate(date=unique(last(pred_quantiles)$fit_df$date)) |> 
    # filter(date %in% recent_dates) |> 
    select(date, mean, contains("quant"))

ggplot(weekly_main, aes(date, mean)) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col=NA, fill="gray80", alpha=0.6) +
    geom_ribbon(aes(ymin=`0.25quant`, ymax=`0.75quant`), col=NA, fill="gray60", alpha=0.6) +
    # annotate(
    #     "segment",
    #     x=min(weekly_summ$date), y=0, xend=max(weekly_summ$date), yend=0,
    #     linetype="13", linewidth=1.05
    # ) +
    geom_line(col="gray40", alpha=0.8) +
    # geom_vline(xintercept=forecast_date, linetype="13", linewidth=1.05) +
    # scale_x_date(breaks=recent_dates, date_labels="%d %b", guide=guide_axis(angle=45)) +
    labs(x=NULL, y=NULL, title="**short-term main effect**")

weekly_inter <- last(pred_quantiles)$fit_df |> 
    arrange(location, date) |> 
    bind_cols(as_tibble(last(pred_quantiles)$fit$summary.random$t2)) |> 
    filter(date >= min(recent_dates)) |> 
    select(location, date, mean, contains("quant"))

ggplot(weekly_inter, aes(date, mean)) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col=NA, fill="gray80", alpha=0.6) +
    geom_ribbon(aes(ymin=`0.25quant`, ymax=`0.75quant`), col=NA, fill="gray60", alpha=0.6) +
    annotate(
        "segment",
        x=min(weekly_summ$date), y=0, xend=max(weekly_summ$date), yend=0,
        linetype="13", linewidth=1.05
    ) +
    geom_line(col="gray40", alpha=0.8) +
    # annotate("text", x=epiweek(forecast_date)+1, y=-0.6, label="latest date", col=hc1) +
    geom_vline(xintercept=forecast_date, linetype="13", linewidth=1.05) +
    facet_wrap(~location, nrow=3) +
    scale_x_date(breaks=recent_dates, date_labels="%d %b", guide=guide_axis(angle=45)) +
    labs(x=NULL, y=NULL, title="**short-term interaction effect**")

# Scratchwork trying to get besag + RW2 to work-----------------------------------

# Things youve tried, keeping in mind the slight difference from before where 
# Alaska etc. are disconnect, and there's a temporal main effect:
# RW1 + besagproper doesnt work with EITHER location included or not
# AR1 + besagproper works
# AR1 + besag (scale.model=TRUE for besag) works, though was maybe slower


graph <- load_us_graph(covid) |> graph2mat()
model <- model_formula("shared", temporal="ar1", spatial="besagproper")

fit4 <- fit_inla_model(fit_df, forecast_date, model, graph=graph)
summary(fit4)

plot(fit4)

hyper_epwk <- list(prec=list(prior="pc.prec", param=c(1, 0.01)))
hyper_wk <- list(prec=list(prior="pc.prec", param=c(1, 0.01)))
ftemp <- count ~ 1 + location + f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, scale.model=TRUE) + f(t, model="rw1", hyper=hyper_wk, scale.model=TRUE)

rtemp <- inla.knmodels(
    ftemp,
    control.st=list(time=t, space=iloc, spacetime=1:nrow(fit_df), graph=graph, type=4),
    data=fit_df,
    family="poisson",
    E=fit_df$ex_lam,
    control.compute=list(dic=FALSE, mlik=FALSE, return.marginals.predictor=TRUE)
)
