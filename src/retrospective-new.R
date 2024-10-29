library(tidyverse)
library(INLA)
library(lubridate)

forecast_retro <- function(
        fit_df, model, forecast_dates, weeks_ahead=4, loc_sub=NULL, 
        q=c(0.025, 0.25, 0.5, 0.75, 0.975), nsamp=5000,...) {
    
    if (is.null(loc_sub))
        loc_sub <- unique(fit_df$location)
    
    map(forecast_dates, \(fd) {
        fd_exact <- min(filter(fit_df, date >= fd)$date)
        fd_end <- fd_exact + weeks(weeks_ahead-1)
        
        # trim the fit_df to contain the training data and just the exact dates/locations to forecast
        pred_df <- fit_df |> 
            filter(date < fd | (date <= fd_end & location %in% loc_sub))
        
        # replace the dates to forecast with NA
        pred_idx <- which(pred_df$date >= fd_exact)
        pred_df$count[pred_idx] <- NA
        
        fit <- fit_inla_model(
            pred_df, fd, model, 
            pred_idx=pred_idx, joint_forecast=FALSE, ...
        )
        
        pred_summ <- sample_count_predictions(pred_df, fit, fd_exact, nsamp=nsamp) |> 
            summarize_quantiles(q=q) |> 
            pivot_wider(names_from=quantile)
        
        list(pred=pred_summ, fit_df=pred_df, fit=fit)
    })
}

update_pred_quantiles <- function(fr, q, nsamp=5000) {
    map(fr, \(el) {
        fd <- min(el$pred$date)
        
        pred_summ <- sample_count_predictions(el$fit_df, el$fit, fd, nsamp=nsamp) |> 
            summarize_quantiles(q=q) |> 
            pivot_wider(names_from=quantile)
        
        list(pred=pred_summ, fit_df=el$fit_df, fit=el$fit)
    })
}

plot_retro <- function(fr, disease_df, loc_sub=NULL) {
    date_min <- min(map_vec(fr, ~min(.$pred$date))) - weeks(4)
    date_max <- max(map_vec(fr, ~max(.$pred$date))) + weeks(4)
    
    if (is.null(loc_sub)) {
        loc_sub <- unique(fr[[1]]$pred$location)
    }
    else {
        loc_sub <- intersect(unique(fr[[1]]$pred$location), loc_sub)
    }
    
    df_sub <- filter(
        disease_df, 
        date >= date_min, 
        date <= date_max,
        location %in% loc_sub
    )
    
    gg <- ggplot(df_sub, aes(date)) +
        facet_wrap(~location, scales="free_y")
        # coord_cartesian_panels(panel_limits=panel_limits) +
        # scale_x_date(breaks=recent_dates[seq(1, length(recent_dates), 2)], date_labels="%d %b", guide=guide_axis(angle=45)) +
    
    for (el in fr) {
        pred_sub <- filter(el$pred, location %in% loc_sub)
        
        gg <- gg +
            geom_ribbon(
                aes(ymin=`0.025`, ymax=`0.975`), 
                pred_sub,
                fill="skyblue", alpha=0.5, col=NA
            ) +
            geom_ribbon(
                aes(ymin=`0.25`, ymax=`0.75`), 
                pred_sub,
                fill="blue1", alpha=0.5, col=NA
            ) +
            geom_line(aes(y=mean), pred_sub, col="blue4")
    }
    gg + geom_point(aes(y=count), size=0.79, alpha=0.8)
}

baseline_epi_df <- function(disease_df) {
    disease_df |> 
        select(geo_value=location, time_value=date, count) |> 
        epiprocess::as_epi_df()
}



#' coverage_timeseries
#' @description
#' Produce a continuous sequence of forecast quantiles within an interval for a specific forecast
#' horizon. Useful for coverage plots.
#'
#' @param fit_df a dataframe of data as produced by `prep_fit_data_<disease>`
#' @param model an INLA model formula
#' @param start start date to begin forecasting. This does not have to be an exact date
#' within `fit_df` - truncation of `fit_df` begins for dates >= `start`
#' @param end end date to stop forecasting. This does not have to be an exact date
#' within `fit_df`
#' @param horizon number of weeks ahead to do each forecast
#' @param loc_sub optionally, limit forecasts to just a vector of locations. Note that model
#' fitting/training still uses all locations in this case
#' @param nsamp number of samples to draw with `sample_count_predictions`
#' @param ... additional arguments to passed to `fit_inla_model`
#'
#' @return a tibble containing forecast quantiles for each date between `start` and `end`
coverage_timeseries <- function(fit_df, model, start, end, horizon=2, loc_sub=NULL, nsamp=2000,...) {
    forecast_dates <- unique(filter(fit_df, date >= start, date <= end)$date)
    if (is.null(loc_sub))
        loc_sub <- unique(fit_df$location)
    
    map_dfr(forecast_dates, \(fd) {
        obs_end_date <- fd - weeks(horizon)
        
        # trim the fit_df to contain the training data and just the exact date to forecast
        pred_df <- fit_df |> 
            filter(date <= obs_end_date | (date == fd & location %in% loc_sub))
        
        # replace the date to forecast with NA
        pred_idx <- which(pred_df$date == fd)
        pred_df$count[pred_idx] <- NA
        
        fit <- fit_inla_model(
            pred_df, fd, model, 
            pred_idx=pred_idx, joint_forecast=FALSE, ...
        )
        
        # TODO: change to a list including fit object
        sample_count_predictions(pred_df, fit, fd, nsamp=nsamp) |> 
            summarize_quantiles() |> 
            pivot_wider(names_from=quantile) |> 
            mutate(last_date=obs_end_date, .after=date)
    })
}
