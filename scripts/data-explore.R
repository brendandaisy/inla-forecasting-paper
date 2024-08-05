library(tidyverse)
library(lubridate)
library(INLA)

rsv <- read_csv("data/weekly-rsv-us.csv")
flu <- read_csv("data/weekly-flu-us.csv")

start_flu <- "2021-09-01" # confirm what we did for FluSight
flu <- filter(flu, date >= start_flu)

ggplot(flu, aes(date, count)) +
    geom_point() +
    facet_wrap(~location, scales="free")

# TODO: figure: departure from average national data, after controlling for seasonality
# basically Figure 1 from Osthus and Moran except starting with residuals from the seasonal effect

mod_season <- count ~ 1 + location + f(epiweek, model="rw2", scale.model=TRUE)
    
    fit <- inla(
        mod, family="poisson", data=fit_df, # poisson regression link
        E=fit_df$ex_lam,
