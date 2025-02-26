# superseded by retrospective-new.R

library(tidyverse)
library(INLA)
library(lubridate)
library(cowplot)

walk(list.files("src", full.names=TRUE), source)
inla.setOption(inla.mode="classic")

coverage_timeseries <- function(fit_df, model, start, end, horizon=2) {
    forecast_dates <- unique(filter(fit_df, date >= start, date <= end)$date)
    
    map_dfr(forecast_dates, \(fd) {
        pred_df <- filter(fit_df, date <= fd)
        obs_end_date <- fd - weeks(horizon-1)
        pred_idx <- which(pred_df$date == fd)
        pred_df$count[which(pred_df$date >= obs_end_date)] <- NA
        
        fit <- fit_inla_model(
            pred_df, fd, model, 
            pred_idx=pred_idx, joint_forecast=FALSE, pc_prior_u=c(1, 0.2)
        )
        
        sample_count_predictions(pred_df, fit, fd, nsamp=2000) |> 
            summarize_quantiles() |> 
            pivot_wider(names_from=quantile)
    })
}

rsv <- read_csv("data/weekly-rsv-us.csv")
covid_interval <- interval(ymd("2020-03-01"), ymd("2021-06-01")) # period during covid to remove

mod_seas <- model_formula("shared", temporal="ar1", spatial="none")
mod_short <- model_formula("none", temporal="ar1", spatial="exchangeable")
mod_both <- model_formula("shared", temporal="ar1", spatial="exchangeable")

fit_df <- rsv |> 
    filter(!(date %within% covid_interval)) |> 
    prep_fit_data_rsv(weeks_ahead=0, ex_lam=pop_served)

ct_seas <- coverage_timeseries(fit_df, mod_seas, "2023-10-15", "2024-04-01", horizon=3)
ct_short <- coverage_timeseries(fit_df, mod_short, "2023-10-15", "2024-04-01", horizon=3)
ct_both <- coverage_timeseries(fit_df, mod_both, "2023-10-15", "2024-04-01", horizon=3)

ct <- bind_rows(
    mutate(ct_seas, model="no interaction term"),
    mutate(ct_short, model="no seasonal term"),
    mutate(ct_both, model="full model")
)

cov_states <- c("California", "New York", "Georgia", "Utah", "Maryland", "Tennessee")

rsv_ct_truth <- rsv |> 
    filter(date >= ymd("2023-10-15")-weeks(4), date < ymd("2024-04-01")+weeks(4)) |> 
    filter(location %in% cov_states)

ct |> 
    filter(location %in% cov_states) |> 
    ggplot(aes(date, fill=model)) +
    # geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), alpha=0.1, col=NA) +
    geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`), alpha=0.5, col=NA) +
    # geom_line(aes(y=mean), col="blue4") +
    geom_point(aes(date, count), rsv_ct_truth, size=0.83, inherit.aes=FALSE) +
    facet_wrap(~location, scales="free_y") +
    theme_minimal()

ggsave("figs/rsv-coverage-3wk.pdf", width=6.5, height=3.8)

# Do it for covid-----------------------------------------------------------------
covid <- read_csv("data/weekly-covid-us.csv")
graph <- load_us_graph(covid) |> sf2mat()

mod_seas <- model_formula("shared", temporal="ar1", spatial="none")
mod_short <- model_formula("none", temporal="ar1", spatial="besagproper")
mod_both <- model_formula("shared", temporal="ar1", spatial="besagproper")

fit_df <- rsv |> 
    filter(!(date %within% covid_interval)) |> 
    prep_fit_data_rsv(weeks_ahead=0, ex_lam=pop_served)

ct_seas <- coverage_timeseries(fit_df, mod_seas, "2023-10-15", "2024-04-01", horizon=3)
ct_short <- coverage_timeseries(fit_df, mod_short, "2023-10-15", "2024-04-01", horizon=3)
ct_both <- coverage_timeseries(fit_df, mod_both, "2023-10-15", "2024-04-01", horizon=3)

ct <- bind_rows(
    mutate(ct_seas, model="no interaction term"),
    mutate(ct_short, model="no seasonal term"),
    mutate(ct_both, model="full model")
)

cov_states <- c("California", "New York", "Georgia", "Utah", "Maryland", "Tennessee")

rsv_ct_truth <- rsv |> 
    filter(date >= ymd("2023-10-15")-weeks(4), date < ymd("2024-04-01")+weeks(4)) |> 
    filter(location %in% cov_states)

ct |> 
    filter(location %in% cov_states) |> 
    ggplot(aes(date, fill=model)) +
    # geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), alpha=0.1, col=NA) +
    geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`), alpha=0.5, col=NA) +
    # geom_line(aes(y=mean), col="blue4") +
    geom_point(aes(date, count), rsv_ct_truth, size=0.83, inherit.aes=FALSE) +
    facet_wrap(~location, scales="free_y") +
    theme_minimal()

ggsave("figs/rsv-coverage-3wk.pdf", width=6.5, height=3.8)
