library(tidyverse)
library(lubridate)

# causes error if both forecast_date and pred_idx are NULL, which makes sense since then 
# you wouldn't be making forecasts
forecast_samples <- function(
        fit_df, fit, 
        forecast_date=NULL, pred_idx=NULL, nsamp=1000
) {
    # nloc <- length(unique(fit_df$iloc))
    pred_idx <- parse_number(fit$selection$names)
    
    ret_df <- fit_df[pred_idx,]
    
    jsamp_fvals <- exp(inla.rjmarginal(nsamp, fit$selection)$samples)
    ex_lam <- ret_df$ex_lam
    pred_dim <- length(ex_lam)
    
    count_samp <- list_transpose(map(1:nsamp, \(samp) { # invert list so have sampled counts for each row
        lambda <- jsamp_fvals[,samp] * ex_lam
        rpois(pred_dim, lambda) # indep. Poisson for each spacetime series
    }))
    
    mutate(ret_df, predicted=count_samp)
}

# A simpler quantile summary not in the FluSight format
summarize_quantiles <- function(pred_samples, nat_samps=NULL, q=c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    pred_samples |> 
        # filter(date >= forecast_date) |> 
        bind_rows(nat_samps) |> 
        unnest(predicted) |> 
        group_by(date, location) |> 
        summarize(
            mean=mean(predicted),
            qs=list(value=quantile(predicted, probs=q)), 
            .groups="drop"
        ) |> 
        unnest_wider(qs) |>
        pivot_longer(contains("%"), names_to="quantile") |> 
        mutate(quantile=parse_number(quantile)/100)
}

pred2forecast_quantile <- function(pred_samples, q=c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    summarize_quantiles(q=q_wis) |> 
        left_join(flu, by=c("date", "location")) |> 
        select(date, location, observed=count, predicted=value, quantile_level=quantile) |> 
        as_forecast_quantile()
}