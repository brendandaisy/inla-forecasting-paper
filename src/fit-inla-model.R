

# forecast_date = the final IN of sample date, i.e. most recent date with available count data. Ignored if joint_forecast=FALSE
# pc_prior_u = vector of length 2 controlling the scale of the seasonal and short-term effects, in that order
# pred_idx = specific indices for which to do joint posterior prediction. Ignored if joint_forecast=TRUE
#
# graph = state adjacency graph, to be used as specified by model
fit_inla_model <- function(
        fit_df, model,
        pc_prior_u=c(1, 1),
        family="poisson",
        response=count,
        pred_idx=NULL,
        forecast_date=NULL,
        q=c(0.025, 0.25, 0.5, 0.75, 0.975),
        graph=NULL, dic=FALSE, config=FALSE
) {
    # the PC priors c(u, a) give the probability a that the standard deviation between weeks exceeds u
    # increasing u increases prior beliefs that there will be large jumps between weeks
    hyper_epwk <- list(prec=list(prior="pc.prec", param=c(pc_prior_u[1], 0.01)))
    hyper_wk <- list(prec=list(prior="pc.prec", param=c(pc_prior_u[2], 0.01)))
    
    mod <- as.formula(model)
    
    if (is.null(pred_idx)) {
        if (is.null(forecast_date)) {
            pred_idx <- which(is.na(pull(fit_df, {{response}})))
            print("forecast_date and pred_idx both NULL. This may have unexpected effects if there are NAs in fit_df prior to when forecasting should begin, but is otherwise fine.")
        } else {
            pred_idx <- which(is.na(pull(fit_df, {{response}})) & fit_df$date > forecast_date)
        }
    }
    
    fit <- inla(
        mod, family=family, data=fit_df, # poisson regression link
        E=fit_df$ex_lam,
        quantiles=q,
        selection=if (length(pred_idx) == 0) NULL else list(Predictor=pred_idx),
        control.fixed=list(prec=1, expand.factor.strategy="inla"),
        control.compute=list(dic=dic, mlik=FALSE, return.marginals.predictor=TRUE, config=config),
        control.predictor=list(link=1) # produce marginal fitted values with default (log) link function
    )
    
    return(fit)
}

# fit_baseline <- function(
#         fit_df, model=baseline_rw1(),
#         forecast_date=NULL,
#         pred_idx=NULL,
#         q=c(0.025, 0.25, 0.5, 0.75, 0.975),
#         dic=FALSE, config=FALSE
# ) {
#     u_emp <- sd(fit_df$count, na.rm=TRUE) # set u to empirical sd
#     hyper_wk <- list(prec=list(prior="pc.prec", param=c(1, 0.01)))
#     
#     mod <- as.formula(model)
#     
#     if (!is.null(forecast_date))
#         pred_idx <- which(fit_df$date > forecast_date)
#     
#     fit <- inla(
#         mod, family="poisson", data=fit_df, # poisson regression link
#         E=fit_df$ex_lam,
#         quantiles=q,
#         selection=if (is.null(pred_idx)) NULL else list(Predictor=pred_idx),
#         control.compute=list(dic=dic, mlik=FALSE, return.marginals.predictor=TRUE, config=config),
#         # control.family=list(hyper=list(prec=list(initial=Inf, fixed=TRUE))),
#         control.predictor=list(link=1)
#     )
#     
#     return(fit)
# }
