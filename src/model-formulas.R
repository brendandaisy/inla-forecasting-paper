library(glue)

# note that 
model_formula <- function(
        seasonal=c("none", "shared", "iid"),
        temporal=c("ar1", "rw2"),
        spatial=c("iid", "exchangeable", "besag")
) {
    seasonal <- case_match(
        seasonal,
        "none" ~ "",
        "shared" ~ 'f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, scale.model=TRUE)',
        "iid" ~ 'f(epiweek, model="rw2", scale.model=TRUE, cyclic=TRUE, 
            hyper=hyper_epwk, group=iloc, control.group=list(model="iid"))'
    )
    
    scaling <- case_match(
        temporal,
        c("rw1", "rw2") ~ "scale.model=TRUE, ",
        "ar1" ~ ""
    )
    
    weekly_main <- glue('f(t, model="{temporal}", {scaling}hyper=hyper_wk)')
    
    # TODO: yes the besagproper model makes the precision difficult to compare to other versions,
    # but we can argue it still should be in the same ballbark which is the best I think we can 
    # do in a reasonable amount of time anyway. Also, we can fix d->1 to eliminate a hyperparam
    weekly_interaction <- case_match(
        spatial,
        "iid" ~ glue('f(t2, model="{temporal}", hyper=hyper_wk, 
                     group=iloc, {scaling}control.group=list(model="iid"))'),
        "exchangeable" ~ glue('f(t2, model="ar1", hyper=hyper_wk, 
                     group=iloc, control.group=list(model="exchangeable"))'),
        "besag" ~ glue('f(iloc, model="besag", hyper=hyper_wk, graph=graph, scale.model=TRUE,
                     group=t2, control.group=list(model="{temporal}"))'), # note the group model automatically scaled
    )
    
    glue("count ~ 1 + location + {seasonal} + {weekly_main} + {weekly_interaction}")
}
