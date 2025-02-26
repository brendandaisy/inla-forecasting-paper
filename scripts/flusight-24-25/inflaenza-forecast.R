library(tidyverse)
library(INLA)
library(lubridate)
library(cowplot)

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")
source("scripts/flusight-24-25/helpers.R")

inla.setOption(inla.mode="classic")

prep_fit_data_flusight <- function(
        disease_df, forecast_date=NULL, weeks_ahead=4, ex_lam=population, pct_reporting_weeks=0
) {
    ret <- disease_df |> 
        filter(abbreviation != "US") |> # make sure US is not in training data
        # group_by(date) |> 
        # mutate(t=cur_group_id(), .after=date) |> # label each date sequentially for time index
        # ungroup() |> 
        mutate(
            iloc=as.numeric(fct_inorder(abbreviation)),
            ex_lam={{ex_lam}}
        )
    
    if (!is.null(forecast_date)) {
        ret <- filter(ret, date <= forecast_date)
        
        pred_df <- expand_grid( # makes pairs of new weeks X location
            tibble(
                date=max(ret$date) + weeks(1:4),
                # t=1:weeks_ahead + max(ret$t),
                epiweek=epiweek(date)
            ),
            distinct(ret, abbreviation, iloc)
        )
        
        location_data <- ret |> 
            filter(date >= forecast_date - weeks(pct_reporting_weeks)) |>
            group_by(abbreviation, location) |> 
            summarize(ex_lam=mean(ex_lam))
        
        pred_df <- left_join(
            pred_df, location_data, 
            by=c("abbreviation"), unmatched="error", relationship="many-to-one"
        )
        
        ret <- bind_rows(ret, pred_df) # add to data for counts to be NAs
    }
    
    # index dates with a sequential time index, accounting for holes in the
    # data and spacing indices properly
    date_seq <- seq.Date(min(ret$date), max(ret$date), "1 week")
    
    date_ind <- tibble(t=seq_along(date_seq), date=date_seq) |> 
        filter(date %in% unique(ret$date))
    
    ret |> 
        left_join(date_ind, by=c("date"), relationship="many-to-one") |> 
        mutate(t2=t) |> 
        arrange(t, iloc)
}

plot_seasonal <- function(fit, forecast_date) {
    summ <- fit$summary.random$epiweek |> 
        as_tibble() |> 
        mutate(season_group=rep(c("Contiguous", "AK", "HI", "PR"), each=52))
    
    ggplot(summ, aes(ID, col=season_group, fill=season_group)) +
        geom_vline(xintercept=epiweek(forecast_date), linetype="dashed", col="gray70") +
        geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), alpha=0.5, col=NA) +
        geom_ribbon(aes(ymin=`0.25quant`, ymax=`0.75quant`), alpha=0.5, col=NA) +
        geom_line(aes(y=mean)) +
        facet_wrap(~season_group, nrow=2) +
        labs(x="epiweek", y="seasonal effect") +
        theme_half_open() +
        theme(legend.position="none")
}

plot_holiday <- function(fit, labels=fct_inorder(levels(fit_df$dist_xmas))) {
    summ <- fit$summary.random$ixmas |> 
        as_tibble() |> 
        mutate(`distance xmas`=labels)
    
    ggplot(summ, aes(`distance xmas`)) +
        geom_linerange(aes(ymin=`0.025quant`, ymax=`0.975quant`)) +
        geom_point(aes(y=mean), size=1.2) +
        scale_x_discrete(guide=guide_axis(angle=25)) +
        labs(y=NULL) +
        theme_half_open()
}

hd <- RSocrata::read.socrata(url="https://data.cdc.gov/resource/mpgq-jmmr.json")
flu0 <- fetch_flu(hd)
flu <- filter(flu0, date >= "2021-09-01")

(forecast_date <- max(flu$date)) # final date to include in training (Saturday before submission day)

# check the percentage reporting in recent weeks:
flu |> 
    filter(date >= forecast_date - weeks(4)) |> 
    ggplot(aes(date, pct_reporting)) +
    geom_col() +
    facet_wrap(~abbreviation, scales="free_y")

fit_df <- prep_fit_data_flusight(
    flu, forecast_date, 
    weeks_ahead=4, ex_lam=pop_weighted, pct_reporting_weeks=0
) |> 
    mutate(
        season_group=case_when(
            abbreviation == "PR" ~ 4,
            abbreviation == "HI" ~ 3,
            abbreviation == "AK" ~ 2,
            TRUE ~ 1
        ),
        dist_xmas=fct_inorder(dist_xmas(date)),
        ixmas=as.numeric(dist_xmas),
        holidays=is_holidays(dist_xmas),
        early_covid=is_early_covid(date)
    ) |> 
    filter(!(is.na(count) & date <= forecast_date)) # TODO: need to remove historic NAs or `forecast_samples` breaks

graph <- load_us_graph(flu) |> 
    sf2mat() |> 
    insert_iso_loc(which(unique(flu$abbreviation) == "PR"))

model <- 'count ~ 1 + abbreviation + early_covid +
f(ixmas, model="rw1", hyper=hyper_epwk, scale.model=TRUE) +
f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, scale.model=TRUE, group=season_group, control.group=list(model="iid")) + 
f(t, model="ar1", hyper=hyper_wk) + f(iloc, model="besagproper", hyper=hyper_wk, graph=graph,
group=t2, control.group=list(model="ar1"))'

fit <- fit_inla_model(fit_df, model, graph=graph)

pred_samp <- forecast_samples(fit_df, fit, nsamp=5000)
pred_samp_us <- aggregate_forecast(pred_samp, tags=tibble_row(abbreviation="US", location="US"))

# pred_samp_us <- pred_samp |> 
#     mutate(predicted=map(predicted, \(p) tibble(sample_id=seq_along(p), predicted=p))) |> 
#     unnest(predicted) |> 
#     group_by(date, horizon, sample_id) |> 
#     summarise(predicted=sum(predicted), .groups="drop") |> 
#     mutate(abbreviation="US", location="US") |> 
#     select(-sample_id) |> 
#     nest(predicted=predicted) |> 
#     mutate(predicted=map(predicted, ~.$predicted))

q_wis <- c(0.01, 0.025, seq(0.05, 0.95, by=0.05), 0.975, 0.99)

(reference_date <- forecast_date + weeks(1)) # official definition of the "current" forecast date

pred_summ <- pred_samp |> 
    bind_rows(pred_samp_us) |> 
    summarize_quantiles(abbreviation, q=q_wis) |> 
    mutate(target = paste0("wk inc flu hosp"),
           reference_date=reference_date,
           horizon=horizon-1,
           target_end_date=reference_date + horizon*7,
           output_type_id=quantile,
           output_type='quantile'
    ) |> 
    # arrange(abbreviation, horizon, quantile) |>
    select(reference_date, target, horizon, target_end_date, abbreviation, location, output_type, output_type_id, value)

write_csv(
    select(pred_summ, -abbreviation), 
    paste0("flusight-predictions/FluSight-forecast-hub/model-output/UGA_flucast-INFLAenza/", reference_date, "-UGA_flucast-INFLAenza.csv")
)

# prev_counts <- flu |> 
#     filter(
#         epiyear == epiyear(forecast_date)-1,
#         epiweek >= epiweek(forecast_date)-7,
#         epiweek <= epiweek(forecast_date)+4
#     ) |> 
#     select(epiweek, location, count) |> 
#     mutate(matching_date=rep(forecast_date + weeks(-7:4), each=53))

prev_counts_epiweek <- fit_df |> 
    filter(date >= (forecast_date - weeks(7))) |> 
    mutate(epiyear=epiyear(date)) |> 
    distinct(epiyear, epiweek)

prev_counts <- flu |> 
    inner_join(mutate(prev_counts_epiweek, epiyear=epiyear-1)) |> 
    mutate(matching_date=rep(forecast_date + weeks(-7:4), each=53))

pred_summ |> 
    pivot_wider(names_from=output_type_id) |> 
    ggplot(aes(target_end_date)) +
    geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), fill="skyblue", alpha=0.5, col=NA) +
    geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`), fill="blue1", alpha=0.5, col=NA) +
    geom_line(aes(y=`0.5`), col="blue4") +
    geom_point(aes(matching_date, count), prev_counts, size=0.8, col="tomato") +
    geom_point(aes(date, y=count), filter(flu, date >= forecast_date - weeks(7)), size=0.8) +
    facet_wrap(~abbreviation, scales="free_y") +
    scale_x_date(date_breaks="2 weeks", date_labels="%b %d", guide=guide_axis(angle=45)) +
    labs(x=NULL, y="hospitalizations") +
    theme_half_open()

ggsave(
    paste0("flusight-predictions/INFLAenza-predictions-", reference_date, ".pdf"), 
    width=12, height=8.5
)

###
plot_seasonal(fit, forecast_date)
plot_holiday(fit)
