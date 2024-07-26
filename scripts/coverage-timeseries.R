library(tidyverse)
library(INLA)
library(lubridate)
library(cowplot)

source("src/prep-fit-dfs.R")

inla.setOption(inla.mode="classic")

# TODO 6/4: one interesting thing we could do in some settings is to use a RW2d rather than Kronecker
# product. This could be appropriate if there is a linear ordering on groups, like with age groups for
# Paraguay. Then we have a 2d regular time x age grid and there's probably more flexibility than the 
# Type IV interaction model

# TODO 6/4: why not just try fitting INLA's seasonal model and see what happens? Could clear some
# things up for ya

# fit_rsv_model <- function(fit_df, pred_idx) {
#     hyper <- list(prec=list(prior="pc.prec", param=c(0.5, 0.01)))
#     hyper_seas <- list(prec=list(prior="pc.prec", param=c(1, 0.01)))
#     
#     mod <- count ~ 1 + location + 
#         # TODO 7/5: epiweek needs to have been cyclic!
#         f(
#             epiweek, model="rw2", scale.model=TRUE, cyclic=TRUE, 
#             hyper=hyper_seas, group=iloc, control.group=list(model="iid")
#         ) +
#         # f(t, model="rw1", scale.model=TRUE, hyper=hyper) +
#         f(t2, model="rw2", scale.model=TRUE, hyper=hyper, group=iloc)
#     
#     fit <- inla(
#         mod, family="poisson", data=mutate(fit_df, t2=t),
#         E=fit_df$ex_lam,
#         selection=list(Predictor=pred_idx),
#         control.compute=list(dic=TRUE, mlik=FALSE, return.marginals.predictor=TRUE),
#         control.predictor=list(link=1) # compute quantiles for NA dates (i.e. do the forecasting)
#     )
# }

coverage_timeseries <- function(fit_df, start, end, horizon=2) {
    forecast_dates <- unique(filter(fit_df, date >= start, date <= end)$date)
    map_dfr(forecast_dates, \(fd) {
        pred_df <- filter(fit_df, date <= fd)
        obs_end_date <- fd - weeks(horizon-1)
        pred_idx <- which(pred_df$date == fd)
        pred_df$count[which(pred_df$date >= obs_end_date)] <- NA
        fit <- fit_rsv_model(pred_df, pred_idx)
        
        sample_count_predictions(pred_df, fit, fd, nsamp=2000) |> 
            summarize_quantiles() |> 
            pivot_wider(names_from=quantile)
    })
}

rsv <- read_csv("data/weekly-rsv-us.csv")
start_date <- ymd("2021-06-01")
forecast_date <- ymd("2023-11-15") # when to start forecasts for this example

fit_df <- rsv |> 
    filter(date >= start_date) |> 
    # filter(location %in% c("California", "New York", "Georgia", "Utah")) |> 
    prep_fit_data_rsv()

ct <- coverage_timeseries(fit_df, "2023-09-15", "2024-02-01", horizon=4)

rsv_ct_truth <- rsv |> 
    filter(date >= ymd("2023-09-15")-weeks(4), date < ymd("2024-02-01")+weeks(4)) |> 
    filter(location %in% c("California", "New York", "Georgia", "Utah"))

ct |> 
    filter(location %in% c("California", "New York", "Georgia", "Utah")) |> 
    ggplot(aes(date)) +
    geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), fill="skyblue", alpha=0.5, col=NA) +
    geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`), fill="blue1", alpha=0.5, col=NA) +
    geom_line(aes(y=mean), col="blue4") +
    geom_point(aes(y=count), rsv_ct_sub, size=0.83) +
    facet_wrap(~location, scales="free_y")

###
summary(fit)

pred <- sample_count_predictions(fit_df_2state, fit, forecast_date-weeks(4), nsamp=5000)
pred_quantiles <- summarize_quantiles(pred) |> 
    pivot_wider(names_from=quantile)

rsv |> 
    filter(date >= forecast_date-weeks(16), date < forecast_date+weeks(6)) |> 
    filter(location %in% c("California", "New York")) |> 
    ggplot(aes(date)) +
    geom_ribbon(
        aes(ymin=`0.025`, ymax=`0.975`), 
        pred_quantiles,
        fill="skyblue", alpha=0.5, col=NA
    ) +
    geom_ribbon(
        aes(ymin=`0.25`, ymax=`0.75`), 
        pred_quantiles,
        fill="blue1", alpha=0.5, col=NA
    ) +
    geom_line(aes(y=mean), pred_quantiles, col="blue4") +
    geom_point(aes(y=count), size=0.83) +
    facet_wrap(~location, scales="free_y")
