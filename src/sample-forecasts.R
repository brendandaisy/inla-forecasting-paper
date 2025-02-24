library(tidyverse)
library(lubridate)

# response needs to be a tidyselect now, because I'm bad at programming
forecast_samples <- function(fit_df, fit, nsamp=1000, response=count, family="poisson") {
    pred_idx <- parse_number(fit$selection$names)
    
    ret_df <- fit_df[pred_idx,]
    
    t_forecast <- max(filter(fit_df, !is.na({{response}}))$t) # find the t corresponding to the forecast date
    ret_df$horizon <- ret_df$t - t_forecast # works for hindcasting too!
    
    jsamp_fvals <- exp(inla.rjmarginal(nsamp, fit$selection)$samples)
    ex_lam <- ret_df$ex_lam
    pred_dim <- length(ex_lam)
    
    rfun <- switch(
        family,
        poisson = rpois,
        beta = function(n, mu) 
            rbeta(n, mu*fit$summary.hyperpar[1, "mean"], fit$summary.hyperpar[1, "mean"]*(1-mu))
    )
    
    pred_samp <- list_transpose(map(1:nsamp, \(samp) { # invert list so have sampled counts for each row
        lambda <- jsamp_fvals[,samp] * ex_lam
        rfun(pred_dim, lambda) # indep. Poisson for each spacetime series
    }))
    
    mutate(ret_df, predicted=pred_samp)
}

aggregate_forecast <- function(pred_samp, ..., fun=sum, tags=tibble(location="All")) {
    pred_samp |> 
        mutate(predicted=map(predicted, \(p) tibble(sample_id=seq_along(p), predicted=p))) |> 
        unnest(predicted) |> 
        group_by(date, horizon, ..., sample_id) |> 
        summarise(predicted=fun(predicted), .groups="drop") |> 
        bind_cols(tags) |> 
        select(-sample_id) |> 
        nest(predicted=predicted) |> 
        mutate(predicted=map(predicted, ~.$predicted))
}

# ... passed to `fit_inla_model`
forecast_baseline <- function(fit_df, model=baseline_rw1(), loc_sub=NULL, ...) {
    if (is.null(loc_sub))
        loc_sub <- unique(fit_df$location)
    
    base_preds <- map_dfr(loc_sub, \(loc) {
        fit_base <- fit_df |>
            filter(location == loc) |> 
            fit_inla_model(model, ...)
        
        fit_df |>
            filter(location == loc) |> 
            forecast_samples(fit_base, nsamp=5000)
    })
    
    bind_rows(base_preds)
}

# A simpler quantile summary not in the FluSight format
# model: optionally add a column with the model name
# ... additional groups to aggregate by or variables to pass through to summarize
summarize_quantiles <- function(pred_samples, ..., nat_samps=NULL, q=c(0.025, 0.25, 0.5, 0.75, 0.975), model=NULL) {
    pred_samples |> 
        # filter(date >= forecast_date) |> 
        bind_rows(nat_samps) |> 
        unnest(predicted) |> 
        group_by(date, location, horizon, ...) |> 
        summarize(
            mean=mean(predicted),
            qs=list(value=quantile(predicted, probs=q)), 
            .groups="drop"
        ) |> 
        unnest_wider(qs) |>
        pivot_longer(contains("%"), names_to="quantile") |> 
        mutate(quantile=parse_number(quantile)/100, model=model)
}

# pred2forecast_quantile <- function(pred_samples, q=c(0.025, 0.25, 0.5, 0.75, 0.975)) {
#     summarize_quantiles(q=q_wis) |> 
#         left_join(flu, by=c("date", "location")) |> 
#         select(date, location, observed=count, predicted=value, quantile_level=quantile) |> 
#         as_forecast_quantile()
# }