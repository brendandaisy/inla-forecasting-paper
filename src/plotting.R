library(tidyverse)
library(INLA)
library(lubridate)

plot_inla <- function(fit, plot_hyper=FALSE) {
    plot(fit, plot.hyperparameters=plot_hyper, plot.random.effects=F, plot.fixed.effects=F, plot.lincomb=F)
    par(mfrow=c(1, 1))
} 