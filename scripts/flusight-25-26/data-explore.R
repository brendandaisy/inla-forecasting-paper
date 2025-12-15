library(tidyverse)
library(lubridate)
library(usmap)

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")
source("scripts/flusight-24-25/helpers.R")

inla.setOption(inla.mode="classic")

hd <- read_csv("~/Downloads/Weekly_Hospital_Respiratory_Data__HRD__Metrics_by_Jurisdiction__National_Healthcare_Safety_Network__NHSN___Preliminary__20241120.csv")
flu0 <- fetch_flu()
flu <- filter(flu0, date >= "2021-09-01")

flu |> 
    filter(date >= "2023-12-01") |> 
    ggplot(aes(date, pct_reporting)) +
    geom_line() +
    facet_wrap(~location, scales="free_y")

flu |> 
    pivot_longer(c(weekly_rate, rate_weighted)) |> 
    ggplot(aes(date, value, col=name)) +
    geom_line() +
    facet_wrap(~location, scales="free")

ggplot(flu, aes(pct_reporting)) +
    geom_histogram() +
    facet_wrap(~epiyear)