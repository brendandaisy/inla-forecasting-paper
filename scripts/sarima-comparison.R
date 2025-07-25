# SARIMA parameters are optimized using the default auto.arima settings, which----
###uses a stepwise search procedure. 

library(tidyverse)
library(forecast)
library(lubridate)

source("src/scoring.R")

forecast_score_sarima <- function(disease_df, fd, loc) {
    ddf_sub <- filter(disease_df, location == loc, date <= fd)
    
    sarima_fit <- auto.arima(ddf_sub$count)
    
    sarima_pred <- forecast(sarima_fit, h=4, level=1-2*q_wis[q_wis < 0.5])
    
    sarima_forecast <- bind_rows(
        format_sarima_pred(sarima_pred$lower, fd, loc, lower=TRUE),
        format_sarima_pred(sarima_pred$upper, fd, loc, lower=FALSE),
        tibble(
            date=fd+weeks(1:4), location=loc, horizon=1:4, model="sarima",
            value=sarima_pred$mean, conf_level=0, quantile=0.5
        )
    ) |> 
        arrange(horizon, quantile) |> 
        mutate(quantile=round(quantile, 3))
    
    score_pred_quantiles(disease_df, sarima_forecast)
}

format_sarima_pred <- function(pred, fd, loc, lower=TRUE) {
    pred |> 
        as_tibble() |> 
        mutate(date=fd+weeks(1:4), location=loc, horizon=1:4, model="sarima") |> 
        pivot_longer(contains("%")) |> 
        mutate(
            conf_level=parse_number(name)/100, 
            quantile=if(lower) (1-conf_level)/2 else 1-(1-conf_level)/2,
            .keep="unused"
        )
}

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

scores_rsv <- expand_grid(fd=retro_forecast_dates_flu_rsv, loc=unique(rsv$location)) |> 
    rowwise() |>
    group_map(~forecast_score_sarima(rsv, .$fd, .$loc)) |> 
    bind_rows()

(summ_rsv <- summarize_scores(scores_rsv))
summarize_scores(scores_rsv, by=c("location"))
(summ_rsv$wis-9.0)/summ_rsv$wis

# "rows containing NA values" warning is just from the final forecast dates
# where not all horizons can be made, but nbd for scoring
scores_flu <- expand_grid(fd=retro_forecast_dates_flu_rsv, loc=unique(flu$location)) |> 
    rowwise() |>
    group_map(~forecast_score_sarima(flu, .$fd, .$loc)) |> 
    bind_rows()

(summ_flu <- summarize_scores(scores_flu))
summarize_scores(scores_flu, by=c("location"))
(summ_flu$wis-40.7)/summ_flu$wis

scores_covid <- expand_grid(fd=retro_forecast_dates_covid, loc=unique(covid$location)) |> 
    rowwise() |>
    group_map(~forecast_score_sarima(covid, .$fd, .$loc)) |> 
    bind_rows()

(summ_covid <- summarize_scores(scores_covid))
summarize_scores(scores_covid, by=c("location"))
(summ_covid$wis-48.7)/summ_covid$wis

### manually checking some forecasts
forecast_date <- ymd("2023-01-15")
ddf <- filter(flu, location == "Georgia", date <= forecast_date)

sarima_fit <- auto.arima(ddf$count)
summary(sarima_fit)

sarima_pred <- forecast(sarima_fit, h=4, level=1-2*q_wis[q_wis < 0.5])

sarima_pred$lower # corresponds to (1-X)/2% quantile
sarima_pred$upper # corresponds to 100-(1-X)/2% quantile

sarima_forecast <- bind_rows(
    format_sarima_pred(sarima_pred$lower, forecast_date, "Georgia", lower=TRUE),
    format_sarima_pred(sarima_pred$upper, forecast_date, "Georgia", lower=FALSE),
    tibble(
        date=forecast_date+weeks(1:4), location="Georgia", horizon=1:4, model="sarima",
        value=sarima_pred$mean, conf_level=0, quantile=0.5
    )
) |> 
    arrange(horizon, quantile) |> 
    mutate(quantile=round(quantile, 3))

sarima_forecast |> 
    select(-conf_level) |> 
    pivot_wider(names_from=quantile, values_from=value) |> 
    ggplot(aes(x=date)) +
    geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), fill="skyblue", alpha=0.5, col=NA) +
    geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`), fill="blue1", alpha=0.5, col=NA) +
    geom_line(aes(y=`0.5`), col="blue4") +
    geom_point(aes(y=count), ddf)

