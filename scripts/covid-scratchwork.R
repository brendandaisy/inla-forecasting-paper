library(tidyverse)
library(INLA)
library(lubridate)

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")

inla.setOption(inla.mode="classic")

covid <- read_csv("data/weekly-covid-us.csv")
start_date <- ymd("2021-06-01")
forecast_date <- ymd("2023-11-15") # when to start forecasts for this example

fit_df_ca <- covid |> 
    # filter(date >= start_date, date < forecast_date) |> 
    prep_fit_data_flu_covid(ex_lam=population) |> 
    filter(location == "California")

ggplot(fit_df_ca, aes(date, count)) +
    geom_point()

# TODO check that the short-term effect looks similar for AR1 in the exchangeable vs. besag versions,
# maybe with a subset of 6 or so states

###
fit_df <- covid |> 
    # filter(date >= start_date, date < forecast_date) |> 
    prep_fit_data_flu_covid(ex_lam=population)

forecast_date <- max(covid$date)

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
