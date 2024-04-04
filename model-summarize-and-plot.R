## Work with Bren's forecasts
library(tidyverse)

# assumes that `forecast_date` matches the forecast date of `fit` (i.e. in `fit_inla_model` call)
sample_count_predictions <- function(fit_df, fit, forecast_date, nsamp=1000) {
    # samp_counts <- map2_dfr(fit$marginals.fitted.values, fit_df$ex_lam, \(m, ex) {
    #     msamp <- pmax(0, inla.rmarginal(nsamp, m)) # sampling sometimes produces very small neg. numbers
    #     ct_samp <- rpois(nsamp, msamp * ex)
    #     tibble_row(count_samp=list(ct_samp))
    # })
    # 
    # return(bind_cols(fit_df, samp_counts))
    nstate <- length(unique(fit_df$snum))
    
    ret_df <- fit_df |>
        filter(date >= forecast_date)
    
    jsamp_fvals <- exp(inla.rjmarginal(nsamp, fit$selection)$samples)
    ex_lam <- ret_df$ex_lam
    
    count_samp <- list_transpose(map(1:nsamp, \(samp) { # invert list so have sampled counts for each row
        lambda <- jsamp_fvals[,samp] * ex_lam
        rpois(length(lambda), lambda) # ind. Poisson for each spacetime series
    }))
    return(mutate(ret_df, count_samp=count_samp))
}

sample_national <- function(fit_df, fit, forecast_date, nsamp=1000) {
    state_info <- distinct(fit_df, location, ex_lam)
    nstate <- nrow(state_info)
    
    ret_df <- fit_df |> 
        filter(date >= forecast_date) |> 
        group_by(date, t, epiweek) |> 
        summarise(population=sum(population), .groups="drop")
    
    jsamp_fvals <- exp(inla.rjmarginal(nsamp, fit$selection)$samples)
    ex_lam <- filter(fit_df, date >= forecast_date)$ex_lam
    
    tslice <- map(1:nrow(ret_df), ~nstate*(.x-1) + 1:nstate) # produce sequence [1:nstate, nstate+1:2nstate, ...]
    
    map_dfr(tslice, \(slice) {
        nat_sum_per_t <- map_dbl(1:nsamp, \(samp) {
            lambda <- jsamp_fvals[slice, samp] * ex_lam
            samp <- rpois(nstate, lambda)
            sum(samp)
        })
        # qs <- quantile(nat_sum_per_t, q)
        # names(qs) <- str_c("q", names(qs))
        tibble_row(location="US", count_samp=list(nat_sum_per_t))
    }) |> 
        bind_cols(ret_df) |> 
        select(date:epiweek, location, population, count_samp)
}

# A simpler quantile summary not in the FluSight format
summarize_quantiles <- function(pred_samples, nat_samps, q) {
    pred_samples |> 
        # filter(date >= forecast_date) |> 
        bind_rows(nat_samps) |> 
        unnest(count_samp) |> 
        group_by(date, location) |> 
        summarize(
            mean=mean(count_samp),
            qs = list(value = quantile(count_samp, probs=q)), 
            .groups="drop"
        ) |> 
        unnest_wider(qs) |>
        pivot_longer(contains("%"), names_to="quantile") |> 
        mutate(quantile=parse_number(quantile)/100)
}

summarize_quantiles_cdc <- function(pred_samples, nat_samps, forecast_date, q) {
    pred_samples |> 
        filter(date >= forecast_date) |> 
        bind_rows(nat_samps) |> 
        unnest(count_samp) |> 
        group_by(date, location) |> 
        summarize(qs = list(value = quantile(count_samp, probs=q)), .groups="drop") |> 
        unnest_wider(qs) |>
        pivot_longer(contains("%"), names_to="quantile") |> 
        mutate(quantile = as.numeric(gsub("[\\%,]", "", quantile))/100) |> 
        mutate(target = paste0("wk inc flu hosp"),
               horizon=as.numeric(as.factor(date)) - 1,
               reference_date = forecast_date,
               target_end_date = forecast_date + horizon*7,
               output_type_id = quantile,
               output_type = 'quantile',
               value = round(value)) %>%
        arrange(location, horizon, quantile) |>
        dplyr::select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value)
}

summarize_rate_change_cdc <- function(pred_samples, nat_samps, forecast_date) {
    pred_samples <- bind_rows(pred_samples, nat_samps)
    
    last_values <- pred_samples |> 
        # get from the maximum date of observed data
        filter(!is.na(count)) |> 
        filter(date == max(date)) |> 
        select(location, count_latest=count)
    
    # get last values for US
    last_values <- bind_rows(last_values, tibble(location="US", count_latest=sum(last_values$count_latest)))
    
    pred_samples |> 
        filter(date >= forecast_date) |> 
        bind_rows(nat_samps) |> 
        unnest(count_samp) |> 
        left_join(last_values, by=c("location")) |> 
        mutate(horizon=as.numeric(as.factor(date)) - 1,
               change = count_samp - count_latest,
               pop_change = change/population * 100000) |> 
        mutate(classification = case_when(horizon == -1 & (abs(change) < 10 | abs(pop_change) < 1) ~ 'stable',
                                          horizon == -1 & (pop_change >= 1 & pop_change < 2) ~ 'increase',
                                          horizon == -1 & (pop_change >= 2) ~ 'large_increase',
                                          horizon == -1 & (pop_change <= -1 & pop_change > -2) ~ 'decrease',
                                          horizon == -1 & (pop_change <= -2) ~ 'large_decrease',
                                          horizon == 0 & (abs(change) < 10 | abs(pop_change) < 1) ~ 'stable',
                                          horizon == 0 & (pop_change >= 1 & pop_change < 3) ~ 'increase',
                                          horizon == 0 & (pop_change >= 3) ~ 'large_increase',
                                          horizon == 0 & (pop_change <= -1 & pop_change > -3) ~ 'decrease',
                                          horizon == 0 & (pop_change <= -3) ~ 'large_decrease',
                                          horizon == 1 & (abs(change) < 10 | abs(pop_change) < 2) ~ 'stable',
                                          horizon == 1 & (pop_change >= 2 & pop_change < 4) ~ 'increase',
                                          horizon == 1 & (pop_change >= 4) ~ 'large_increase',
                                          horizon == 1 & (pop_change <= -2 & pop_change > -4) ~ 'decrease',
                                          horizon == 1 & (pop_change <= -4) ~ 'large_decrease',
                                          horizon == 2 & (abs(change) < 10 | abs(pop_change) < 2.5) ~ 'stable',
                                          horizon == 2 & (pop_change >= 2.5 & pop_change < 5) ~ 'increase',
                                          horizon == 2 & (pop_change >= 5) ~ 'large_increase',
                                          horizon == 2 & (pop_change <= -2.5 & pop_change > -5) ~ 'decrease',
                                          horizon == 2 & (pop_change <= -5) ~ 'large_decrease',
                                          horizon == 3 & (abs(change) < 10 | abs(pop_change) < 2.5) ~ 'stable',
                                          horizon == 3 & (pop_change >= 2.5 & pop_change < 5) ~ 'increase',
                                          horizon == 3 & (pop_change >= 5) ~ 'large_increase',
                                          horizon == 3 & (pop_change <= -2.5 & pop_change > -5) ~ 'decrease',
                                          horizon == 3 & (pop_change <= -5) ~ 'large_decrease',
                                          T ~ NA
        )) |> 
        count(location, horizon, classification) |> 
        group_by(location, horizon) |>
        mutate(value = n/sum(n)) |>
        ungroup() |>
        complete(location, horizon, classification, fill=list(value=0)) |>
        mutate(target = 'wk flu hosp rate change',
               reference_date = forecast_date,
               target_end_date = forecast_date + horizon*7,
               output_type_id = classification,
               output_type = 'pmf',
               value = value) %>%
        dplyr::select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value)
}

# Make plots --------------------------------------------------------------
plot_state_forecast <- function(location_name, curr_season_data, forecast_df) {
    curr_df <- curr_season_data |> 
        filter(location == location_name)
    
    forecast_df <- filter(forecast_df, location == location_name)
    # ts_df <- filter(sampled_ts, location == location_name)
    
    ggplot(forecast_df, aes(target_end_date, `0.5`)) +
        geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), alpha = .2) +
        geom_ribbon(aes(ymin = `0.25`, ymax = `0.75`), alpha = .2) +
        geom_line() +
        # geom_line(
        #     aes(date, count_samp, group=samp_id, col=as.factor(samp_id)), ts_df, 
        #     alpha=0.5, show.legend=FALSE
        # ) +
        geom_point(data = curr_df, aes(date, count)) +
        labs(title = location_name, x = NULL, y ='Admits') +
        background_grid(major = 'xy', minor = 'y') +
        coord_cartesian(ylim = c(0, max(c(curr_df$count, forecast_df$`0.75`))))
}

plot_state_pmf <- function(location_name, curr_season_data, ratechange_df) {
    curr_df <- curr_season_data |> 
        filter(location == location_name)
    
    ratechange_df <- filter(ratechange_df, location == location_name)
    
    pred_date_start <- min(ratechange_df$target_end_date)
    pred_date_end <- max(ratechange_df$target_end_date)
    
    p1 <- ggplot(curr_df, aes(date, count)) +
        geom_point() +
        guides(x=guide_axis(angle=40)) +
        labs(title = location_name, x = NULL, y ='Admits') +
        background_grid(major = 'xy', minor = 'y') +
        theme(plot.margin=margin(7, 0, 7, 7))
    
    p2 <- ggplot(ratechange_df, aes(target_end_date, value, fill=output_type_id)) +
        geom_col() +
        labs(x=NULL, y=NULL) +
        scale_y_continuous(expand=expansion()) +
        scale_x_date(breaks=c(pred_date_start, pred_date_end), guide=guide_axis(angle=40), date_labels="%b %d") +
        background_grid(major="none") +
        theme(
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            # axis.line.y=element_line(color="gray85"),
            axis.line.y=element_blank(),
            plot.margin=margin(7, 7, 7, 0),
            legend.position="none"
        )
    
    plot_grid(p1, plot_grid(NULL, p2, rel_heights=c(0.1, 1), nrow=2), nrow=1, rel_widths=c(1, 0.3), axis="b")
}

pmf_legend <- function(ratechange_df) {
    ratechange_df <- filter(ratechange_df, location == "Texas")
    
    p <- ggplot(ratechange_df, aes(target_end_date, value, fill=output_type_id)) +
        geom_col() +
        labs(fill=NULL) +
        theme(legend.position="right", legend.box.margin=margin())
    
    get_legend(p)
}

# assumes the most recent week of data is in first row, for current epiweek
plot_seasonal <- function(df, fit, forecast_date) {
    ep_wk_fx <- fit$summary.random$epiweek |>
        as_tibble()
    # mutate(state=rep(unique(dat_pred$state), each=53)) # assumes states are in correct order in dat_pred
    
    ggplot(ep_wk_fx, aes(ID, mean)) +
        geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col="gray70", alpha=0.6) +
        geom_line(col="tomato3") +
        geom_vline(xintercept=epiweek(forecast_date), col="steelblue", linetype="dashed") +
        # facet_wrap(~state) +
        labs(x="Week of the year", y="Seasonal effect (log scale)") +
        theme_bw() +
        theme(legend.position="none")
}

plot_seasonal_states <- function(df, fit, forecast_date) {
    states <- unique(df$location)
    nepi <- length(unique(df$epiweek))
    
    ep_wk_fx <- fit$summary.random$epiweek |>
        as_tibble() |> 
        mutate(state=rep(states, each=nepi)) # assumes states are in correct order in dat_pred
    
    ggplot(ep_wk_fx, aes(ID, mean)) +
        geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col="gray70", alpha=0.6) +
        geom_line(col="tomato3") +
        geom_vline(xintercept=epiweek(forecast_date), col="steelblue", linetype="dashed") +
        facet_wrap(~state) +
        labs(x="Week of the year", y="Seasonal effect (log scale)") +
        theme_bw() +
        theme(legend.position="none")
}
