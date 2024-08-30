library(tidyverse)
library(lubridate)

# TODO: the current way of specifying ret_df with the forecast data is pretty 
# unsafe, since sometimes forecast_date is used as the first date to start
# forecasting, but also you use it as the last/only date to forecast for the coverage plots
sample_count_predictions <- function(fit_df, fit, forecast_date, nsamp=1000) {
    nloc <- length(unique(fit_df$iloc))
    
    ret_df <- fit_df |>
        filter(date >= forecast_date)
    
    jsamp_fvals <- exp(inla.rjmarginal(nsamp, fit$selection)$samples)
    ex_lam <- ret_df$ex_lam
    
    count_samp <- list_transpose(map(1:nsamp, \(samp) { # invert list so have sampled counts for each row
        lambda <- jsamp_fvals[,samp] * ex_lam
        rpois(length(lambda), lambda) # indep. Poisson for each spacetime series
    }))
    
    ret_df |> 
        mutate(count_samp=count_samp) |> 
        select(where(~!any(is.na(.x))))
}

# A simpler quantile summary not in the FluSight format
summarize_quantiles <- function(pred_samples, nat_samps=NULL, q=c(0.025, 0.25, 0.5, 0.75, 0.975)) {
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