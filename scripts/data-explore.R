library(tidyverse)
library(lubridate)
library(INLA)
library(corrr)
library(sf)
library(igraph)
library(cowplot)

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")

inla.setOption(inla.mode="classic")

us_dist_mat <- function(us_sf) {
    ig <- poly2nb(us_sf) |> 
        c() |> 
        map(\(nb) if (any(nb == 0)) integer(0) else nb) |>  # replace empty neighbors with an integer(0)
        igraph::graph_from_adj_list()
    
    ret <- distances(ig)
    rownames(ret) <- unique(us_sf$state)
    colnames(ret) <- unique(us_sf$state)
    
    ret |> 
        as_tibble(rownames="state_from") |> 
        pivot_longer(-state_from, names_to="state_to", values_to="dist")
}

rsv <- read_csv("data/weekly-rsv-us.csv") |> 
    filter(!(date %within% interval(ymd("2020-03-01"), ymd("2021-06-01")))) |> 
    mutate(season_week=(epiweek-40) %% 52 + 1) # epiweek 40 is start of RSV season, going off provided season column

flu <- read_csv("data/weekly-flu-us.csv") |> 
    filter(
        location != "US", # won't need nat'l for these figures
        date > "2021-09-01" # No real point in looking at data before around here
    ) |> 
    mutate(
        season=case_when(
            date < "2021-09-01" ~ "2020-21",
            date >= "2021-09-01" & date < "2022-09-01" ~ "2021-22",
            date >= "2022-09-01" & date < "2023-09-01" ~ "2022-23",
            date >= "2023-09-01" ~ "2023-24"
        ),
        season_week=(epiweek-35) %% 52 + 1
    )

covid <- read_csv("data/weekly-covid-us.csv") |> 
    filter(location != "US") |> 
    mutate(
        season=case_when( # just define as same as flu season
            date < "2021-09-01" ~ "2020-21",
            date >= "2021-09-01" & date < "2022-09-01" ~ "2021-22",
            date >= "2022-09-01" & date < "2023-09-01" ~ "2022-23",
            date >= "2023-09-01" ~ "2023-24"
        ),
        season_week=(epiweek-35) %% 52 + 1
    )

# distinct(flu, season, epiweek, season_week) |> filter(season_week %in% c(1, 52))

ggplot(covid, aes(date, count)) +
    geom_point() +
    facet_wrap(~location, scales="free")

# For Flu: want to illustrate that the US, and to a larger extent many individual states, there
# was a much sharper peak in hospitalizations in 22-23 following lull from COVID,
# contrary to other findings about flu severity
# Importantly for this paper: is there structure to these differences? The southeast coastal states
# are consistently the ones which break this pattern of the two peaks

decompose_timeseries <- function(data, us_dist) {
    nat_seas <- data |> 
        group_by(epiweek, season_week) |> 
        summarise(mean=mean(weekly_rate), .groups="drop")
    
    resid_seas <- data |> 
        left_join(nat_seas, by=c("epiweek", "season_week"), relationship="many-to-one") |> 
        mutate(resid_seas=weekly_rate-mean)
    
    # resid_seas_diff <- resid_seas |> 
    #     group_by(location, season) |> 
    #     summarise(seas_mae=mean(abs(resid_seas)), .groups="drop") |> 
    #     slice_max(seas_mae, n=4)
    
    resid_seas_summ <- resid_seas |> 
        group_by(date, season) |> 
        summarise(week_mean=mean(resid_seas))
    
    resid_seas_corr <- resid_seas |> 
        pivot_wider(id_cols=date, names_from=location, values_from=resid_seas) |> 
        correlate() |> 
        stretch(na.rm=TRUE, remove.dups=TRUE) |> 
        rename(state_from=x, state_to=y) |> 
        left_join(us_dist, by=c("state_from", "state_to"))
    
    return(list(
        nat_seas=nat_seas,
        resid_seas=resid_seas,
        resid_seas_summ=resid_seas_summ,
        resid_seas_corr=resid_seas_corr
    ))
}

plot_disease_summary <- function(dts, data, disease=c("RSV", "flu", "COVID-19"), col="tomato") {
    data_sub <- filter(data, season_week <= 35)
    
    p1 <- ggplot(data_sub, aes(season_week, weekly_rate, group=interaction(season, location))) +
        geom_line(alpha=0.25) +
        # geom_line(aes(col=interaction(season, location)), data=inner_join(flu_sub, frs_diff, by=c("location", "season"))) +
        geom_line(aes(season_week, mean), filter(dts$nat_seas, season_week <= 35), col=col, linewidth=1.02, inherit.aes=FALSE) +
        labs(x="season week", y=glue("weekly rate {disease} / 100k")) +
        scale_x_continuous(breaks=seq(0, 35 , 5)) +
        theme_half_open() +
        background_grid(major="xy")
    
    
    if (disease == "RSV") {
        # grouping for disconnected timeseries (e.g. early RSV seasons)
        aes1 <- aes(date, resid_seas, group=interaction(location, season))
        aes2 <- aes(date, week_mean, group=season)
    }
    else {
        aes1 <- aes(date, resid_seas, group=location)
        aes2 <- aes(date, week_mean)
    }
    
    p2 <- ggplot(dts$resid_seas, aes1) + 
        geom_line(alpha=0.25) +
        geom_line(aes2, dts$resid_seas_summ, col=col, linewidth=1.02, inherit.aes=FALSE) +
        # geom_line(aes(col=location), filter(flu_resid_seas, location %in% hl_states), linewidth=1.01)
        scale_x_date(date_breaks="6 months", date_labels="%b '%y", guide=guide_axis(angle=45)) +
        labs(x=NULL, y="residual rate / 100k") +
        theme_half_open() +
        background_grid(major="xy")
    
    corr_summ <- dts$resid_seas_corr |> 
        filter(!is.na(r), is.finite(dist)) |> 
        group_by(dist) |> 
        summarise(med=median(r), l=med - 1.58*IQR(r)/sqrt(n()), u=med + 1.58*IQR(r)/sqrt(n()))
    
    # TODO 8/20: using loess is a bit disingenuous but looks nicer
    p3 <- dts$resid_seas_corr |> 
        filter(!is.na(r), is.finite(dist)) |> 
        ggplot(aes(factor(dist), r)) +
        geom_point(alpha=0.25, shape=1) +
        # geom_boxplot(col="tomato", fill=NA, outlier.shape=NA, alpha=0.5) +
        geom_line(aes(dist, med), corr_summ, col=col, linewidth=1.02) +
        geom_ribbon(aes(dist, ymin=l, ymax=u), corr_summ, linetype="dotted", col=col, fill=NA, inherit.aes=FALSE) +
        # stat_boxplot(aes(y=after_stat(xlower)), geom="line", linetype="dotted", col=col, linewidth=1.02) +
        # stat_boxplot(aes(x=dist, y=after_stat(notchupper)), geom="line", linetype="dotted", col=col, linewidth=1.02) +
        # geom_smooth(aes(dist, r), se=FALSE, col="tomato", method="loess") +
        labs(x="neighborhood distance", y="residual correlation") +
        theme_half_open()
    
    plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1, 1.4, 1), align="h", axis="b")
}

box_ribbon <- function(x) {
    n <- length(x)
    tibble(
        y=median(x),
        ymin=quantile(x, 0.1),
        ymax=quantile(x, 0.9)
        # ymin=y - 1.58*IQR(x)/sqrt(n),
        # ymax=y + 1.58*IQR(x)/sqrt(n),
    )
}

# hl_states <- c("Indiana", "Arkansas", "South Carolina", "North Carolina")

us <- load_us_graph(flu)
us_dist <- us_dist_mat(us)

dts_flu <- decompose_timeseries(flu, us_dist)
p1 <- plot_disease_summary(dts_flu, flu, "flu")

dts_rsv <- decompose_timeseries(rsv, us_dist)
p2 <- plot_disease_summary(dts_rsv, rsv, "RSV")

dts_covid <- decompose_timeseries(covid, us_dist)
p3 <- plot_disease_summary(dts_covid, covid, "COVID-19")

plot_grid(p1, p2, p3, nrow=3)
ggsave("figs/figure-1-data-summ.pdf", width=11, height=7.5)

# frs_corr |> 
#     filter(r >= 0.85) |> 
#     arrange(dist) |> 
#     View()

# Just do covid separately:
nat_seas_c19 <- covid |> 
    group_by(epiweek) |> 
    summarise(mean=mean(weekly_rate), .groups="drop")

resid_seas_c19 <- covid |> 
    left_join(nat_seas_c19, by="epiweek", relationship="many-to-one") |> 
    mutate(resid_seas=weekly_rate-mean)

resid_seas_summ_c19 <- resid_seas_c19 |> 
    group_by(date, epiyear) |> 
    summarise(week_mean=mean(resid_seas))

resid_seas_corr_c19 <- resid_seas_c19 |> 
    pivot_wider(id_cols=date, names_from=location, values_from=resid_seas) |> 
    correlate() |> 
    stretch(na.rm=TRUE, remove.dups=TRUE) |> 
    rename(state_from=x, state_to=y) |> 
    left_join(us_dist, by=c("state_from", "state_to"))

covid |> 
    mutate(epiweek2=(epiweek - 35) %% 53) |> 
    ggplot(aes(epiweek2, weekly_rate, group=interaction(location, epiyear))) +
    geom_line(alpha=0.25) +
    # geom_line(aes(col=interaction(season, location)), data=inner_join(flu_sub, frs_diff, by=c("location", "season"))) +
    geom_line(aes(epiweek2, mean), mutate(nat_seas_c19, epiweek2=(epiweek - 35) %% 53), col="tomato", inherit.aes=FALSE) +
    labs(x="epiweek", y="weekly rate covid / 100k") +
    scale_x_continuous(breaks=seq(0, 50, 5), labels=as.character((seq(0, 50, 5) + 35) %% 53)) +
    theme_half_open() +
    background_grid(major="xy")

# grouping for disconnected timeseries (e.g. early RSV seasons)
p2 <- ggplot(dts$resid_seas, aes(date, resid_seas, group=interaction(location, season))) + 
    geom_line(alpha=0.25) +
    geom_line(aes(date, week_mean, group=season), dts$resid_seas_summ, col="tomato", inherit.aes=FALSE) +
    # geom_line(aes(col=location), filter(flu_resid_seas, location %in% hl_states), linewidth=1.01)
    scale_x_date(date_breaks="3 months", date_labels="%b '%y", guide=guide_axis(angle=45)) +
    labs(x=NULL, y="residual rate / 100k") +
    theme_half_open() +
    background_grid(major="xy")

# TODO 8/20: using loess is a bit disingenuous but looks nicer
p3 <- dts$resid_seas_corr |> 
    filter(!is.na(r), is.finite(dist)) |> 
    ggplot(aes(factor(dist), r)) +
    geom_point(alpha=0.25, shape=1) +
    # geom_boxplot(col="tomato", fill=NA, outlier.shape=NA, alpha=0.5) +
    # stat_summary(aes(x=dist), fun="median", col="tomato", geom="line") +
    # stat_boxplot(aes(x=dist, y=after_stat(lower)), geom="line", linetype="dotted", col="tomato") +
    geom_smooth(aes(dist, r), se=FALSE, col="tomato", method="loess") +
    labs(x="neighborhood distance", y="residual correlation") +
    theme_half_open()

plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1, 1.4, 1), align="h", axis="b")

###
fit_df <- flu |> 
    prep_fit_data_flu_covid(weeks_ahead=0)

hyper_epwk <- list(prec=list(prior="pc.prec", param=c(1, 0.01)))
mod_season <- count ~ 1 + location + 
    f(epiweek, model="rw2", hyper=hyper_epwk, scale.model=TRUE, cyclic=TRUE)

# if (joint_forecast)
#     pred_idx <- which(fit_df$date >= forecast_date)

fit <- inla(
    mod_season, family="poisson", data=fit_df,
    E=fit_df$ex_lam,
    # selection=if (is.null(pred_idx)) NULL else list(Predictor=pred_idx),
    control.compute=list(mlik=FALSE, return.marginals.predictor=TRUE),
    control.predictor=list(link=1) # produce marginal fitted values with default (log) link function
)

plot_seasonal(fit_df, fit)


###

# fit the final flu model to the full training set
model <- model_formula("shared", temporal="ar1", spatial="besagproper")
graph <- load_us_graph(covid) |> graph2mat()

# fit the model and sample predictions for each of the timepoints
fit_df <- flu |> 
    prep_fit_data_flu_covid(weeks_ahead=0, ex_lam=population)
    
fit <- fit_inla_model(fit_df, NULL, model, pc_prior_u=c(1, 1), graph=graph, joint_forecast=FALSE)

summary(fit)
seasonal_mean <- tibble(epiweek=1:52, mean=fit$summary.random$epiweek$mean)

flu |> 
    filter(location != "US") |> 
    left_join(seasonal_mean, by="epiweek") |> 
    filter(season_week <= 35) |> 
    mutate(trans_mean=100000*exp(-12.065 + mean)) |> 
    ggplot(aes(season_week, weekly_rate, group=interaction(season, location))) +
    geom_line(alpha=0.1) +
    geom_line(aes(season_week, trans_mean), col="red", inherit.aes=FALSE) +
    scale_y_continuous(sec.axis=sec_axis(trans=~log(100000) - 12.065 + .x))

peaks_flu <- flu |> 
    group_by(location, epiyear) |> 
    slice_max(count, n=1) |> 
    slice_min(season_week, n=1) |> # break ties according to earlier in the season
    ungroup()

# TODO would need to think about how to add the peak seasonal effect to y axis
# e.g. take exp(global intercept + max seasonal) and trans. to weekly rate
ggplot(peaks_flu, aes(season_week, weekly_rate, col=as.factor(epiyear))) +
    geom_point()
