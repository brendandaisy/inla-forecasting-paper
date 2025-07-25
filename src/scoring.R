library(INLA)
library(tidyverse)
library(lubridate)
library(scoringutils)

score_pred_quantiles <-  function(disease_df, pred_quantiles, ...) {
    fc_quant <- pred_quantiles |>
        left_join(disease_df, by=c("date", "location")) |> # add truth data back in
        select(date, location, horizon, model, observed=count, predicted=value, quantile_level=quantile, ...) |>
        as_forecast_quantile()

    interval_coverage_95 <- purrr::partial(interval_coverage, interval_range=95)
    score(fc_quant, metrics=list(wis=wis, pic_50=interval_coverage, pic_95=interval_coverage_95))
}
