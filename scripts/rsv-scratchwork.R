library(tidyverse)
library(INLA)
library(lubridate)
library(cowplot)

# TODO 7/15: ask Spencer: is there evidence that seasonality for flu etc. follows regional
# patterns in the US? Would be strong justification to try sep. seasonality just for 
# each region

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")

inla.setOption(inla.mode="classic")

rsv <- read_csv("data/weekly-rsv-us.csv")

forecast_dates <- c("2023-09-20", "2023-11-15", "2023-12-30", "2024-02-01")
covid_interval <- interval(ymd("2020-03-01"), ymd("2021-06-01")) # period during covid to remove

# new way to specify the current "base" model
model <- model_formula("shared", temporal="ar1", spatial="exchangeable")

# fit the model and sample predictions for each of the timepoints
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

#  Plot the retrospective forecasts compared to truth-----------------------------
# panel_limits <- pred_quantiles |> 
#     group_by(location) |> 
#     summarise(ymin=0, ymax=0.75*max(`0.975`))

last_date <- ymd(max(forecast_dates))
rsv_sub <- filter(rsv, date >= (last_date-weeks(20)), date < (last_date+weeks(8)))

recent_dates <- unique(rsv_sub$date)

gg <- ggplot(rsv_sub, aes(date)) +
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

# Other plots for model outputs---------------------------------------------------
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

#  Old scratchwork for a single state---------------------------------------------
plot_inla <- function(fit, plot_hyper=FALSE) {
    plot(fit, plot.hyperparameters=plot_hyper, plot.random.effects=F, plot.fixed.effects=F, plot.lincomb=F)
    par(mfrow=c(1, 1))
}

fit_df_ca <- rsv |> 
    filter(date >= "2021-06-01", date < "2023-12-01") |> 
    prep_fit_data_rsv(weeks_ahead=4, ex_lam=pop_served) |> 
    filter(location == "California")

ggplot(fit_df_ca, aes(date, count)) +
    geom_point()

ggplot(tibble(sd=sqrt(1/inla.pc.rprec(10000, 0.01, 0.01))), aes(sd)) +
    geom_density()
    # coord_cartesian(xlim=c(NA, 100000))

ggplot(tibble(sd=sqrt(1/rgamma(10000, 100, 5e-5))), aes(sd)) +
    geom_density()

hyper <- list(prec=list(prior="pc.prec", param=c(0.05, 0.01)))
hyper_fix <- list(prec=list(initial=0.1, fixed=TRUE))

mod_rw <- count ~ 1 + f(t, model="rw2")
mod_ar <- count ~ 1 + f(t, model="ar1", hyper=hyper)

mod_sea <- count ~ 1 + f(t, model="seasonal", season.length=52, scale.model=TRUE, hyper=hyper) + f(t2, model="rw1", scale.model=TRUE)

pred_idx <- which(fit_df_ca$date >= "2023-12-01")

fit <- inla(
    mod_ar, family="poisson", data=fit_df_ca,
    E=fit_df_ca$ex_lam,
    # selection=list(Predictor=pred_idx),
    control.compute=list(dic=TRUE, mlik=FALSE, return.marginals.predictor=TRUE),
    control.predictor=list(link=1) # compute quantiles for NA dates (i.e. do the forecasting)
)

summary(fit)
plot_inla(fit)

# TODO 8/21: the univariate case still clearly shows high prediction uncertainty with the 
# default prior, after data change and the problem mostly went away for the full model.
# Maybe we have the tools now to illustrate the issue as a Supplementary Note?

post_sd_upt05 <- inla.tmarginal(function(x) sqrt(1/x), fit$marginals.hyperpar$`Precision for t`)

post_sd <- bind_rows(
    mutate(as_tibble(post_sd_upt05), u="u=0.05"),
    mutate(as_tibble(post_sd_upt2), u="u=0.2"),
    mutate(as_tibble(post_sd_u1), u="u=1"),
    mutate(as_tibble(post_sd_u10), u="u=10"),
)

(pal <- c("gray70", "black", "goldenrod1", "tomato"))

p1 <- ggplot(post_sd) +
    geom_line(aes(x, y, col=u), linewidth=1.03) +
    stat_function(fun=~inla.pc.dprec(1/.x^2, 0.05, 0.1)*abs(-2*.x^(-3)), col=pal[1], linetype="31", linewidth=1.03) +
    stat_function(fun=~inla.pc.dprec(1/.x^2, 0.2, 0.1)*abs(-2*.x^(-3)), col=pal[2], linetype="31", linewidth=1.03) +
    stat_function(fun=~inla.pc.dprec(1/.x^2, 1, 0.1)*abs(-2*.x^(-3)), col=pal[3], linetype="31", linewidth=1.03) +
    stat_function(fun=~inla.pc.dprec(1/.x^2, 10, 0.1)*abs(-2*.x^(-3)), col=pal[4], linetype="31", linewidth=1.03) +
    scale_color_manual(values=pal) +
    xlim(0, 3) +
    coord_cartesian(ylim=c(0, 5)) +
    labs(x="AR standard deviation", y=NULL, col=NULL) +
    theme_half_open()

covid <- read_csv("data/weekly-covid-us.csv")

ca_log_rates <- covid |> 
    filter(location == "California") |> 
    mutate(disease="COVID-19", log_rate=log(count/population + 1/(2*population)))

ca_log_rates <- rsv |> 
    filter(location == "California") |> 
    mutate(disease="RSV", log_rate=log(count/pop_served + 1/(2*pop_served))) |> 
    bind_rows(ca_log_rates)

log_summ <- ca_log_rates |> 
    group_by(disease) |> 
    summarize(mean=mean(log_rate), med=median(log_rate), sd=sd(log_rate))

library(ggtext)

p2 <- ggplot(ca_log_rates) +
    stat_function(fun=log, xlim=c(1e-8, 0.05), n=500, linewidth=1.03) +
    geom_jitter(aes(x=exp(log_rate), y=log_rate, col=disease), height=0, width=0.005, shape=1, alpha=0.35) +
    geom_errorbar(aes(x=0.014, ymin=mean-sd, ymax=mean+sd, col=disease), log_summ, width=0.004) +
    geom_point(aes(x=0.014, y=mean, col=disease), log_summ) +
    geom_text(
        aes(x=0.018, y=mean, label=lab, col=disease), 
        mutate(log_summ, lab=c("sigma == 0.7", "sigma == 1.8")),
        size=6, hjust="left", show.legend=FALSE, parse=TRUE
    ) +
    scale_color_manual(values=c("gray70", "goldenrod1")) +
    coord_cartesian(xlim=c(-0.02, 0.05), ylim=c(-16, -2)) +
    labs(x="rate per capita", y="log rate", col=NULL) +
    theme_half_open()
    # theme(axis.text.x=element_text(size=7))

plot_grid(p1, p2, nrow=1, rel_widths=c(1, 1), align="h", axis="b", labels="AUTO")

ggsave("figs/prior-precision-rsv.pdf", width=9.7, height=4.3)
