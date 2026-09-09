# --------------------------------------------------------------------------------
# Supplemental figure showing adjustment of the precision of the temporal---------
# effects for RSV-----------------------------------------------------------------
# --------------------------------------------------------------------------------
library(tidyverse)
library(INLA)
library(lubridate)
library(cowplot)

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")

inla.setOption(inla.mode="classic")

rsv <- read_csv("data/weekly-rsv-us.csv")

fit_df_ca <- rsv |> 
    filter(date >= "2021-06-01", date < "2023-12-01") |> 
    prep_fit_data_rsv(weeks_ahead=4, ex_lam=pop_served) |> 
    filter(location == "California")

post_sd_ar <- function(u) {
    hyper <- list(prec=list(prior="pc.prec", param=c(u, 0.01)))
    
    mod_ar <- count ~ 1 + f(t, model="ar1", hyper=hyper)
    
    pred_idx <- which(fit_df_ca$date >= "2023-12-01")
    
    fit <- inla(
        mod_ar, family="poisson", data=fit_df_ca,
        E=fit_df_ca$ex_lam,
        # selection=list(Predictor=pred_idx),
        control.compute=list(dic=TRUE, mlik=FALSE, return.marginals.predictor=TRUE),
        control.predictor=list(link=1) # compute quantiles for NA dates (i.e. do the forecasting)
    )
    
    inla.tmarginal(function(x) sqrt(1/x), fit$marginals.hyperpar$`Precision for t`)   
}

post_sd <- bind_rows(
    mutate(as_tibble(post_sd_ar(0.05)), u="u=0.05"),
    mutate(as_tibble(post_sd_ar(0.2)), u="u=0.2"),
    mutate(as_tibble(post_sd_ar(1)), u="u=1"),
    mutate(as_tibble(post_sd_ar(10)), u="u=10"),
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
    labs(x="Rate per capita", y="Log rate", col=NULL) +
    theme_half_open()
# theme(axis.text.x=element_text(size=7))

plot_grid(p1, p2, nrow=1, rel_widths=c(1, 1), align="h", axis="b", labels="auto")

ggsave("figs/prior-precision-rsv.png", width=9.7, height=4.3, bg="white")
