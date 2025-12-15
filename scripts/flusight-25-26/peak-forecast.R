library(tidyverse)
library(INLA)
library(lubridate)
library(cowplot)

inla.setOption(inla.mode="classic")

hd <- RSocrata::read.socrata(url="https://data.cdc.gov/resource/mpgq-jmmr.json")
flu0 <- fetch_flu(hd)
flu <- filter(flu0, date >= "2021-09-01")

forecast_date <- ymd("2024-11-23") # final date to include in training

flu_seas <- flu |> 
    mutate(
        season=case_when(
            date < "2021-09-01" ~ "2020-21",
            date >= "2021-09-01" & date < "2022-09-01" ~ "2021-22",
            date >= "2022-09-01" & date < "2023-09-01" ~ "2022-23",
            date >= "2023-09-01" & date < "2024-09-01" ~ "2023-24",
            date >= "2024-09-01" ~ "2024-25"
        ),
        season_week=(epiweek-35) %% 52 + 1
    )

flu_peak <- flu_seas |> 
    group_by(season, location) |> 
    filter(date < "2024-09-01") |> 
    summarise(
        peak_intensity=max(rate_weighted), 
        peak_time=epiweek[which.max(rate_weighted)],
        peak_season_time=season_week[which.max(rate_weighted)]
    )

epi_labs <- seq(min(flu_peak$peak_time), max(flu_peak$peak_time), 5)

ggplot(flu_peak, aes(peak_season_time, (peak_intensity^.5 - 1)/.5)) +
    geom_point()

                       