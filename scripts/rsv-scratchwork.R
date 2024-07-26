library(tidyverse)
library(INLA)
library(lubridate)
library(cowplot)

source("coord_cartesian_panels.R")

plot_inla <- function(fit, plot_hyper=FALSE) {
    plot(fit, plot.hyperparameters=plot_hyper, plot.random.effects=F, plot.fixed.effects=F, plot.lincomb=F)
    par(mfrow=c(1, 1))
}   

inla.setOption(inla.mode="classic")

# TODO 7/15: ask Spencer: is there evidence that seasonality for flu etc. follows regional
# patterns in the US? Would be strong justification to try sep. seasonality just for 
# each region

rsv <- read_csv("data/weekly-rsv-us.csv")

start_date <- ymd("2021-06-01")
forecast_date <- ymd("2023-11-20") # when to start forecasts for this example

forecast_dates <- c("2023-09-20", "2023-11-15", "2023-12-30", "2024-02-01")
model <- model_formula("shared", temporal="ar1", spatial="exchangeable")

covid_interval <- interval(ymd("2020-03-01"), ymd("2021-06-01"))

pred_quantiles <- map(forecast_dates, \(fd) {
    fit_df <- rsv |> 
        filter(!(date %within% covid_interval), date <= fd) |> 
        prep_fit_data_rsv(weeks_ahead=4, ex_lam=pop_served)
    
    fit <- fit_inla_model(fit_df, fd, model, pc_prior_u=c(1, 0.2))
    
    pred_summ <- sample_count_predictions(fit_df, fit, fd, nsamp=5000) |> 
        summarize_quantiles() |> 
        pivot_wider(names_from=quantile)
    
    list(pred=pred_summ, fit_df=fit_df, fit=fit)
})

# panel_limits <- pred_quantiles |> 
#     group_by(location) |> 
#     summarise(ymin=0, ymax=0.75*max(`0.975`))

last_date <- ymd(max(forecast_dates))
rsv_recent <- filter(rsv, date >= (last_date-weeks(20)), date < (last_date+weeks(8)))

recent_dates <- unique(rsv_recent$date)

gg <- ggplot(rsv_recent, aes(date)) +
    geom_point(aes(y=count), size=0.83) +
    facet_wrap(~location, scales="free_y") +
    # coord_cartesian_panels(panel_limits=panel_limits) +
    scale_x_date(breaks=recent_dates[seq(1, length(recent_dates), 2)], date_labels="%d %b", guide=guide_axis(angle=45)) +
    theme(title=element_text(size=12))

for (el in pred_quantiles) {
    gg <- gg +
        geom_ribbon(
            aes(ymin=`0.025`, ymax=`0.975`), 
            el$pred,
            fill="skyblue", alpha=0.5, col=NA
        ) +
        geom_ribbon(
            aes(ymin=`0.25`, ymax=`0.75`), 
            el$pred,
            fill="blue1", alpha=0.5, col=NA
        ) +
        geom_line(aes(y=mean), el$pred, col="blue4")
}
gg

plot_seasonal(fit_df, fit)

# TODO make plot function

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
    # theme_mod_output()

###
fit_df_ca <- rsv |> 
    filter(date >= start_date, date < forecast_date) |> 
    prep_fit_data_rsv(ex_lam=pop_served) |> 
    filter(location == "California")

ggplot(fit_df_ca, aes(date, count)) +
    geom_point()

ggplot(tibble(sd=sqrt(1/inla.pc.rprec(10000, 0.01, 0.01))), aes(sd)) +
    geom_density()
    # coord_cartesian(xlim=c(NA, 100000))

ggplot(tibble(sd=sqrt(1/rgamma(10000, 100, 5e-5))), aes(sd)) +
    geom_density()

hyper <- list(prec=list(prior="pc.prec", param=c(1, 0.01)))
hyper_fix <- list(prec=list(initial=0.1, fixed=TRUE))

mod_rw <- count ~ 1 + f(epiweek, model="rw2", cyclic=TRUE) + f(t, model="rw2", hyper=hyper)
mod_ar <- count ~ 1 + f(t, model="ar", order=3, hyper=hyper)

mod_sea <- count ~ 1 + f(t, model="seasonal", season.length=52, scale.model=TRUE, hyper=hyper) + f(t2, model="rw1", scale.model=TRUE)

# TODO 6/27: add weights to, for example, turn off the COVID years for the seasonal effect but use
# previous years. Perhaps: f(epiweek, non_covid_year, model="rw2") + f(epiweek2, covid_year, ...)
# Although ultimately it's probably better to just remove the covid years and index with t

pred_idx <- which(fit_df_ca$date >= forecast_date)

fit <- inla(
    mod_sea, family="poisson", data=mutate(fit_df_ca, t2=t),
    E=fit_df_ca$ex_lam,
    # selection=list(Predictor=pred_idx),
    control.compute=list(dic=TRUE, mlik=FALSE, return.marginals.predictor=TRUE),
    control.predictor=list(link=1) # compute quantiles for NA dates (i.e. do the forecasting)
)

summary(fit)
plot_inla(fit)

# TODO: this has nothing to do with posterior prediction!!!
# the posterior covariance matrix for AR _absolutely_ has correlation all through the timeseries
inla_rar1 <- function(x0, kappa, rho, t=1:5, n=5000) {
    tau <- kappa / (1 - rho^2)
    ret <- map(1:n, \(i) {
        as.vector(arima.sim(list(order=c(1, 0, 0), ar=rho), n=4, sd=sqrt(1/tau)))
        # accumulate(t, \(u, i) rnorm(1, rho*u, 1/sqrt(tau)), .init=x0)
    })
    do.call(rbind, ret) |> as_tibble()
    # tau <- kappa / (1 - rho^2)
    # accumulate(t, \(u, i) rnorm(1, rho*u, 1/sqrt(tau)), .init=x0)
}

u <- inla_rar1(-13, 2, 0.89, n=10000)
summarize(u, across(everything(), ~median(.)))
inla_rar1(2e-5, 0.339, 0.976)

# fit$summary.random$t |> 
#     as_tibble() |> 
#     mutate(t=fit_df$t, location=fit_df$location) |> 
#     ggplot(aes(t)) +
#     geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), fill="gray70", col=NA, alpha=0.7) +
#     geom_line(aes(y=mean)) +
#     coord_cartesian(xlim=c(100, NA)) +
#     facet_wrap(~location, scales="free_y") +
#     theme_cowplot()
