library(tidyverse)
library(INLA)
library(lubridate)
library(ggthemes)

plot_inla <- function(fit, plot_hyper=FALSE) {
    plot(fit, plot.hyperparameters=plot_hyper, plot.random.effects=F, plot.fixed.effects=F, plot.lincomb=F)
    par(mfrow=c(1, 1))
}

theme_mod_output <- function() {
    theme_minimal() +
        theme(
            panel.grid.minor=element_blank()
            # plot.title=element_markdown(size=13),
            # axis.title.y=element_markdown(),
            # strip.text=element_text(size=12)
        )
}

plot_seasonal <- function(fit, forecast_date, seasonal=c("shared", "iid")) {
    seasonal_summ <- fit$summary.random$epiweek |> 
        as_tibble() |> 
        select(epiweek=ID, mean, contains("quant"))
    
    ret <- ggplot(seasonal_summ, aes(epiweek, mean)) +
        geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col=NA, fill="gray80", alpha=0.6) +
        geom_ribbon(aes(ymin=`0.25quant`, ymax=`0.75quant`), col=NA, fill="gray60", alpha=0.6) +
        annotate(
            "segment", 
            x=0, y=0, xend=52, yend=0, 
            col="#F4A247", linetype="13", linewidth=1.05
        ) +
        geom_line(col="gray40", alpha=0.8) +
        geom_vline(xintercept=epiweek(forecast_date), col="#046CE9", linetype="13", linewidth=1.05) +
        labs(x="epiweek", y="seasonal effect") +
        theme_mod_output()
    
    ret
}
