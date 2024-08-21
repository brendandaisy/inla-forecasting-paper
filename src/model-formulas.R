library(glue)

# note that 
model_formula <- function(
        seasonal=c("none", "shared", "iid"),
        temporal=c("none", "ar1", "rw2"),
        spatial=c("none", "iid", "exchangeable", "besag", "besagproper") # 7/26 you noticed besagproper was indeed faster and also had slighly smaller prediction intervals for covid :)
) {
    seasonal <- case_match(
        seasonal,
        "none" ~ "",
        "shared" ~ '+ f(epiweek, model="rw2", cyclic=TRUE, hyper=hyper_epwk, scale.model=TRUE)',
        "iid" ~ '+ f(epiweek, model="rw2", scale.model=TRUE, cyclic=TRUE, 
            hyper=hyper_epwk, group=iloc, control.group=list(model="iid"))'
    )
    
    scaling <- case_match(
        temporal,
        c("rw1", "rw2") ~ "scale.model=TRUE, ",
        "ar1" ~ ""
    )
    
    weekly_main <- case_match(
        temporal,
        "none" ~ "",
        c("ar1", "rw2") ~ glue('+ f(t, model="{temporal}", {scaling}hyper=hyper_wk)')
    )
    
    weekly_interaction <- case_match(
        spatial,
        "none" ~ "",
        "iid" ~ glue('+ f(t2, model="{temporal}", hyper=hyper_wk, 
                     group=iloc, {scaling}control.group=list(model="iid"))'),
        "exchangeable" ~ glue('+ f(t2, model="{temporal}", hyper=hyper_wk, 
                     group=iloc, {scaling}control.group=list(model="exchangeable"))'),
        "besag" ~ glue('+ f(iloc, model="besag", hyper=hyper_wk, graph=graph, scale.model=TRUE,
                     group=t2, control.group=list(model="{temporal}"))'), # note the group model automatically scaled
        "besagproper" ~ glue('+ f(iloc, model="besagproper", hyper=hyper_wk, graph=graph,
                     group=t2, control.group=list(model="{temporal}"))')
    )
    
    glue("count ~ 1 + location {seasonal} {weekly_main} {weekly_interaction}")
}
