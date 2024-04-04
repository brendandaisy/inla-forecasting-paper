library(tidyverse)
library(INLA)
library(sf)
library(spdep)
library(lubridate)
library(cowplot)

source("fetch-data/read-rsv-data.R")
source("model-prep-and-fit.R")
source("model-summarize-and-plot.R")

# Read in current data
rsv0 <- fetch_rsv_full()
location_info <- distinct(rsv0, location, location_name) # get the location coding for saving final results
forecast_date <- ymd("2022-11-05") # first week where forecasting will take place
prediction_date <- forecast_date - duration(6, "weeks") # week to start prediction of counts, for comparing fit to recent training data

rsv <- rsv0 |> 
    select(-c(location, age_group, target), location=location_name, count=value) |> 
    filter(date < forecast_date)

fit_df <- prep_fit_data(rsv, weeks_ahead=4)
model <- rsv_model_new()
fit <- fit_inla_model(fit_df, prediction_date, model, joint=TRUE, pc_prior_u=c(0.5, 1))

summary(fit)
plot_seasonal(fit_df, fit, forecast_date)
plot_seasonal_states(fit_df, fit, forecast_date)

fit_df |> 
    arrange(snum, t) |> 
    mutate(
        mean=fit$summary.random$t$mean,
        ymin=fit$summary.random$t$`0.025quant`,
        ymax=fit$summary.random$t$`0.975quant`
    ) |> 
    filter(date >= "2022-06-01") |> 
    ggplot(aes(date)) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), col=NA, alpha=0.3) +
    geom_line(aes(y=mean), col="gray30") +
    facet_wrap(~location) +
    geom_vline(xintercept=forecast_date-duration(1, "week"), col="plum2", linetype="31", linewidth=1.2, alpha=0.7) +
    labs(y=NULL, x="week", title="weekly effect") +
    scale_x_date(guide=guide_axis(angle=30)) +
    theme_half_open() +
    theme(legend.position="none")

# fitted values
fit_df |> 
    mutate(
        mean=fit$summary.fitted.values$mean,
        ymin=fit$summary.fitted.values$`0.025quant`,
        ymax=fit$summary.fitted.values$`0.975quant`
    ) |> 
    filter(date >= "2022-06-01") |> 
    ggplot(aes(date)) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), col=NA, alpha=0.3) +
    geom_line(aes(y=mean), col="gray30") +
    facet_wrap(~location) +
    geom_vline(xintercept=forecast_date-duration(1, "week"), col="plum2", linetype="31", linewidth=1.2, alpha=0.7) +
    labs(y=NULL, x="week", title="fitted values") +
    scale_x_date(guide=guide_axis(angle=30)) +
    theme_half_open() +
    theme(legend.position="none")

pred_samples <- sample_count_predictions(fit_df, fit, prediction_date, nsamp=2000)

pred_summ <- summarize_quantiles(pred_samples, NULL, q=c(0.025, 0.25, 0.5, 0.75, 0.975)) |> 
    pivot_wider(names_from=quantile, values_from=value)

fit_df |> 
    filter(date >= "2021-04-01") |>
    ggplot(aes(date)) +
    geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), pred_summ, col=NA, alpha=0.3, fill="lightblue") +
    geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`), pred_summ, col=NA, alpha=0.9, fill="lightblue") +
    geom_line(aes(y=mean), pred_summ, col="blue") +
    geom_point(aes(y=count), size=1.1, col="black") +
    facet_wrap(~location, scales="free_y") +
    scale_x_date(date_labels="%Y-%m", guide=guide_axis(angle=30)) +
    # scale_y_continuous(trans="log1p") +
    labs(x="Week", y="RSV Cases") +
    theme_half_open() +
    theme(legend.position="none")

ggsave("figs/rsv-model-3-15-reg-priors.pdf", width=8, height=6)

# Compare forecast variance more thoroughly--------------------------------------
qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)

mod0 <- rsv_model_exchangeable()
fit0 <- fit_inla_model(fit_df, prediction_date, mod0, q=qs, joint=FALSE)
fit0b <- fit_inla_model(fit_df, prediction_date, mod0, q=qs, joint=FALSE, pc_prior_u=c(0.5, 0.05))
fit0c <- fit_inla_model(fit_df, prediction_date, mod0, q=qs, joint=FALSE, pc_prior_u=c(0.05, 0.05))

mod1 <- 'count ~ 1 + location +
    f(epiweek, model="rw2", hyper=hyper_epwk, scale.model=TRUE, cyclic=TRUE) +
    f(t, model="ar1", hyper=hyper_wk, group=snum, control.group=list(model="iid"))'
fit1 <- fit_inla_model(fit_df, prediction_date, mod1, q=qs, joint=FALSE)
fit1b <- fit_inla_model(fit_df, prediction_date, mod1, q=qs, joint=FALSE, pc_prior_u=c(0.5, 0.05))
fit1c <- fit_inla_model(fit_df, prediction_date, mod1, q=qs, joint=FALSE, pc_prior_u=c(0.05, 0.05))

mod2 <- 'count ~ 1 + location +
    f(epiweek, model="rw2", hyper=hyper_epwk, scale.model=TRUE, cyclic=TRUE) +
    f(t, model="ar", hyper=hyper_wk, group=snum, order=4, control.group=list(model="exchangeable"))'
fit2 <- fit_inla_model(fit_df, prediction_date, mod2, q=qs, joint=FALSE)
fit2b <- fit_inla_model(fit_df, prediction_date, mod2, q=qs, joint=FALSE, pc_prior_u=c(0.5, 0.05))
fit2c <- fit_inla_model(fit_df, prediction_date, mod2, q=qs, joint=FALSE, pc_prior_u=c(0.05, 0.05))
# mod1 <- rsv_model_exchangeable(epiweek=FALSE)
# fit1 <- fit_inla_model(fit_df, prediction_date, mod1, q=qs, joint=FALSE)
# fit1b <- fit_inla_model(fit_df, prediction_date, mod1, q=qs, joint=FALSE, pc_prior_u=c(NA, 0.05))


mod3 <- 'count ~ 1 + location +
    f(epiweek, model="rw2", hyper=hyper_epwk, group=snum, control.group=list(model="iid"), scale.model=TRUE) +
    f(t, model="ar1", hyper=hyper_wk, group=snum, control.group=list(model="exchangeable"))'
fit3 <- fit_inla_model(fit_df, prediction_date, mod3, q=qs, joint=FALSE)
fit3b <- fit_inla_model(fit_df, prediction_date, mod3, q=qs, joint=FALSE, pc_prior_u=c(0.5, 0.05))
fit3c <- fit_inla_model(fit_df, prediction_date, mod3, q=qs, joint=FALSE, pc_prior_u=c(0.05, 0.05))

mod4 <- 'count ~ 1 + location +
    f(epiweek, model="rw2", hyper=hyper_epwk, group=snum, control.group=list(model="iid"), scale.model=TRUE) +
    f(t, model="ar", hyper=hyper_wk, order=4, group=snum, control.group=list(model="exchangeable"))'
fit4 <- fit_inla_model(fit_df, prediction_date, mod4, q=qs, joint=FALSE)
fit4b <- fit_inla_model(fit_df, prediction_date, mod4, q=qs, joint=FALSE, pc_prior_u=c(0.5, 0.05))
fit4c <- fit_inla_model(fit_df, prediction_date, mod4, q=qs, joint=FALSE, pc_prior_u=c(0.05, 0.05))

fits <- list(fit0, fit0b, fit0c, fit1, fit1b, fit1c, fit2, fit2b, fit2c, fit3, fit3b, fit3c, fit4, fit4b, fit4c)
names(fits) <- c("fit0", "fit0b", "fit0c", "fit1", "fit1b", "fit1c", "fit2", "fit2b", "fit2c", "fit3", "fit3b", "fit3c", "fit4", "fit4b", "fit4c")

pred_idx <- which(is.na(fit_df$count))
pred_uppers <- imap_dfr(fits, \(fit_obj, name) {
    fit_df |> 
        slice(pred_idx) |> 
        select(date, location) |> 
        mutate(
            type=str_extract(name, "\\d"),
            model=name,
            lin_pred1=fit_obj$summary.linear.predictor$`0.75quant`[pred_idx],
            lin_pred2=fit_obj$summary.linear.predictor$`0.975quant`[pred_idx],
            fit_val1=fit_obj$summary.fitted.values$`0.75quant`[pred_idx],
            fit_val2=fit_obj$summary.fitted.values$`0.975quant`[pred_idx]
        )
})

# pred_uppers |> 
#     ggplot(aes(model, fit_val2, fill=date)) +
#     geom_col() +
#     facet_wrap(~location, scales="free")

# remove model 1 since it was strictly worse than model 0
pred_uppers |> 
    filter(type != "1") |> 
    ggplot(aes(model, fit_val2, fill=as.factor(date))) +
    geom_col() +
    facet_wrap(~location, scales="free") +
    scale_x_discrete(guide=guide_axis(angle=45))

ggsave("figs/rsv-fit-val-upper-3-15.pdf", width=10.5, height=6.5)
