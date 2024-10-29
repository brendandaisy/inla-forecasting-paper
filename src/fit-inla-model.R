# forecast_date = the final IN of sample date, i.e. most recent date with available count data. Ignored if joint_forecast=FALSE
# pc_prior_u = vector of length 2 controlling the scale of the seasonal and short-term effects, in that order
# pred_idx = specific indices for which to do joint posterior prediction. Ignored if joint_forecast=TRUE
#
# graph = state adjacency graph, to be used as specified by model
fit_inla_model <- function(
        fit_df, model,
        forecast_date=NULL,
        pc_prior_u=c(1, 1),
        pred_idx=NULL,
        q=c(0.025, 0.25, 0.5, 0.75, 0.975),
        graph=NULL, dic=FALSE, config=FALSE
) {
    # the PC priors c(u, a) give the probability a that the standard deviation between weeks exceeds u
    # increasing u increases prior beliefs that there will be large jumps between weeks
    hyper_epwk <- list(prec=list(prior="pc.prec", param=c(pc_prior_u[1], 0.01)))
    hyper_wk <- list(prec=list(prior="pc.prec", param=c(pc_prior_u[2], 0.01)))
    
    mod <- as.formula(model)
    
    if (!is.null(forecast_date))
        pred_idx <- which(fit_df$date > forecast_date)
    
    fit <- inla(
        mod, family="poisson", data=fit_df, # poisson regression link
        E=fit_df$ex_lam,
        quantiles=q,
        selection=if (is.null(pred_idx)) NULL else list(Predictor=pred_idx),
        control.fixed=list(prec=1, expand.factor.strategy="inla"),
        control.compute=list(dic=dic, mlik=FALSE, return.marginals.predictor=TRUE, config=config),
        control.predictor=list(link=1) # produce marginal fitted values with default (log) link function
    )
    
    return(fit)
}
