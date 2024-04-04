library(tidyverse)
library(INLA)
library(lubridate)

inla.setOption(inla.mode="classic")

load_us_graph <- function(rsv, f="data/us-state-boundaries.shp") {
  us0 <- read_sf(f)
  
  us <- us0 |> 
    filter(name %in% unique(rsv$location)) |> 
    select(state=name, division, region)
  
  state_order <- fct_inorder(unique(rsv$location))
  
  # sort order of states to match their order of appearance in rsv data
  us |> 
    mutate(state=factor(state, levels(state_order))) |> 
    arrange(state)
}

us_adj_mat <- function(us) {
  us_adj <- us |> 
    poly2nb() |> 
    nb2mat(style="W", zero.policy=TRUE)
  
  us_adj[which(us$state == "Florida"), which(us$state == "Puerto Rico")] <- 1
  us_adj[which(us$state == "Puerto Rico"), which(us$state == "Florida")] <- 1
  # us_adj[which(us$state == "Georgia"), which(us$state == "Puerto Rico")] <- 1
  # us_adj[which(us$state == "Puerto Rico"), which(us$state == "Georgia")] <- 1
  us_adj[which(us$state == "Washington"), which(us$state == "Alaska")] <- 1
  us_adj[which(us$state == "Alaska"), which(us$state == "Washington")] <- 1
  # us_adj[which(us$state == "California"), which(us$state == "Alaska")] <- 1
  # us_adj[which(us$state == "Alaska"), which(us$state == "California")] <- 1
  us_adj[which(us$state == "California"), which(us$state == "Hawaii")] <- 1
  us_adj[which(us$state == "Hawaii"), which(us$state == "California")] <- 1
  return(us_adj)
}

# TODO: try to make this robust to a much wider range of scenarios
prep_fit_data <- function(flu, weeks_ahead=4) {
    ret <- flu |> 
        filter(location != "US") |> # make sure US is not in training data
        group_by(date) |> 
        mutate(t=cur_group_id(), epiweek=epiweek(date), .after=date) |> # add a time counter starting from 1 for earliest week
        ungroup() |> 
        mutate(
            snum=as.numeric(fct_inorder(location)), # INLA needs groups as ints starting from 1, so add numeric state code
            ex_lam=population
        )
    
    # make a dataframe to hold group info for forecasting
    pred_df <- expand_grid(
        tibble(
            date=duration(1:weeks_ahead, "week") + max(ret$date),
            t=1:weeks_ahead + max(ret$t),
            epiweek=epiweek(date)
        ),
        distinct(ret, location, snum) # makes pairs of new times X each state
    )
    
    # go and find most recent population/ex_lam values for each state
    state_data <- ret |> 
        slice_max(date) |> 
        distinct(location, population, ex_lam)
    
    # add on most recent pop/ex_lam values, with protection against adding/deleting unintended rows
    pred_df <- left_join(
        pred_df, state_data, 
        by=c("location"), unmatched="error", relationship="many-to-one"
    )
    
    bind_rows(ret, pred_df) |> # add to data for counts to be NAs
        arrange(t)
}

# forecast_date = the first OUT of sample date, i.e. date to start predictions
# graph = state adjacency graph, to be used as specified by model
fit_inla_model <- function(
        fit_df, forecast_date, model=current_rsv_model(),
        q=c(0.025, 0.5, 0.975), pc_prior_u=c(0.5, 1),
        graph=NULL, joint=TRUE, diagnostics=FALSE
    ) {
    # the PC priors c(u, a) give the probability that the standard deviation between weeks u is greater than a
    # increasing u increases prior beliefs that there will be large jumps between weeks
    hyper_epwk <- list(prec=list(prior="pc.prec", param=c(pc_prior_u[1], 0.01)))
    hyper_wk <- list(prec=list(prior="pc.prec", param=c(pc_prior_u[2], 0.01)))
    
    mod <- as.formula(model)
    
    pred_idx <- which(fit_df$date >= forecast_date)
    fit <- inla(
        mod, family="poisson", data=fit_df, # poisson regression link
        E=fit_df$ex_lam,
        quantiles=q,
        selection=if (joint) list(Predictor=pred_idx) else NULL,
        control.compute=c(
            if (diagnostics) list(dic=TRUE, mlik=TRUE) else list(mlik=FALSE),
            list(return.marginals.predictor=TRUE)
        ),
        control.predictor=list(link=1) # produce marginal fitted values with default (log) link function
    )
    return(fit)
}

rsv_model_new <- function(epiweek=TRUE, week=TRUE) {
    #        f(epiweek, model="ar", hyper=hyper_epwk, order=4) +
    'count ~ 1 + location +
        f(epiweek, model="rw2", hyper=hyper_epwk, group=snum, control.group=list(model="iid"), scale.model=TRUE) +
        f(t, model="ar1", hyper=hyper_wk, group=snum, control.group=list(model="exchangeable"))'
}

rsv_model_exchangeable <- function(epiweek=TRUE, week=TRUE) {
    paste0(
        'count ~ 1 + location + ',
        if (epiweek) 'f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, scale.model=TRUE) + ' else "",
        if (week) 'f(t, model="ar1", hyper=hyper_wk, group=snum, control.group=list(model="exchangeable"))' else ""
    )
}

# assumes the following variables are available in the environment:
# epiweek, snum, week, graph, and lists storing hyperpriors for the seasonal and weekly
# random effects, `hyper_epwk` and `hyper_wk` respectively
rsv_model_pcar <- function(epiweek=TRUE, week=TRUE) {
    paste0(
        'count ~ 1 + location + ',
        if (epiweek) 'f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, scale.model=TRUE) + ' else "",
        if (week) 'f(snum, model="besagproper", graph=graph, hyper=hyper_wk, group=t, control.group=list(model="ar1"))' else ""
    )
}

rsv_model_2 <- function() {
    'count ~ 1 + location + f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, scale.model=TRUE) +
    f(t, model="ar1", hyper=hyper_wk, group=snum, control.group=list(model="exchangeable")) +
    f(tt, model="rw1", hyper=hyper_wk, group=snum, control.group=list(model="iid"))'
}
