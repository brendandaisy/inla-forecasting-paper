# forecast_date = the first OUT of sample date, i.e. date to start predictions
# pc_prior_u = vector of length 2 controlling the scale of the seasonal and short-term effects, in that order
# pred_idx = specific indices for which to do joint posterior prediction. Ignored if joint_forecast=TRUE
#
# graph = state adjacency graph, to be used as specified by model
fit_inla_model <- function(
        fit_df, forecast_date, model,
        pc_prior_u=c(1, 1),
        joint_forecast=TRUE,
        pred_idx=NULL,
        q=c(0.025, 0.25, 0.5, 0.75, 0.975),
        graph=NULL, dic=FALSE
) {
    # the PC priors c(u, a) give the probability that the standard deviation between weeks u is greater than a
    # increasing u increases prior beliefs that there will be large jumps between weeks
    hyper_epwk <- list(prec=list(prior="pc.prec", param=c(pc_prior_u[1], 0.01)))
    hyper_wk <- list(prec=list(prior="pc.prec", param=c(pc_prior_u[2], 0.01)))
    
    mod <- as.formula(model)
    
    if (joint_forecast)
        pred_idx <- which(fit_df$date >= forecast_date)
    
    fit <- inla(
        mod, family="poisson", data=fit_df, # poisson regression link
        E=fit_df$ex_lam,
        quantiles=q,
        selection=if (is.null(pred_idx)) NULL else list(Predictor=pred_idx),
        control.compute=list(dic=dic, mlik=FALSE, return.marginals.predictor=TRUE),
        control.predictor=list(link=1) # produce marginal fitted values with default (log) link function
    )
    
    return(fit)
}
