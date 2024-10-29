library(tidyverse)
library(INLA)
library(lubridate)
library(cowplot)
library(sf)

source("src/prep-fit-data.R")
source("src/model-formulas.R")
source("src/fit-inla-model.R")
source("src/sample-forecasts.R")

inla.setOption(inla.mode="classic")

rsv <- read_csv("data/weekly-rsv-us.csv") |> 
    filter(!(date %within% interval(ymd("2020-03-01"), ymd("2021-06-01")))) |> 
    mutate(season_week=(epiweek-40) %% 52 + 1) # epiweek 40 is start of RSV season, going off provided season column

flu <- read_csv("data/weekly-flu-us.csv") |> 
    filter(
        location != "US", # won't need nat'l for these figures
        date > "2021-09-01" # No real point in looking at data before around here
    ) |> 
    mutate(
        season=case_when(
            date < "2021-09-01" ~ "2020-21",
            date >= "2021-09-01" & date < "2022-09-01" ~ "2021-22",
            date >= "2022-09-01" & date < "2023-09-01" ~ "2022-23",
            date >= "2023-09-01" ~ "2023-24"
        ),
        season_week=(epiweek-35) %% 52 + 1
    )

covid <- read_csv("data/weekly-covid-us.csv") |> 
    filter(location != "US") |> 
    mutate(
        season=case_when( # just define as same as flu season
            date < "2021-09-01" ~ "2020-21",
            date >= "2021-09-01" & date < "2022-09-01" ~ "2021-22",
            date >= "2022-09-01" & date < "2023-09-01" ~ "2022-23",
            date >= "2023-09-01" ~ "2023-24"
        ),
        season_week=(epiweek-35) %% 52 + 1
    )

graph <- load_us_graph(flu) |> sf2mat()
model <- model_formula(seasonal="shared", temporal="ar1", spatial="besagproper")

fit_df_flu <- prep_fit_data_flu_covid(flu, weeks_ahead=0, ex_lam=population) |> 
    mutate(location=factor(location))
fit_flu <- fit_inla_model(fit_df_flu, model, graph=graph, config=TRUE) # config was for trying conditional sampling

fit_df_rsv <- prep_fit_data_rsv(rsv, weeks_ahead=0, ex_lam=pop_served) |> 
    mutate(location=factor(location))
fit_rsv <- fit_inla_model(fit_df_rsv, model_formula(seasonal="shared", temporal="ar1", spatial="exchangeable"))

fit_df_covid <- prep_fit_data_flu_covid(covid, weeks_ahead=0, ex_lam=population) |> 
    mutate(location=factor(location))
fit_covid <- fit_inla_model(fit_df_covid, model, graph=graph)


#  option 1 for spatial main effects: US map--------------------------------------
plot_loc_effect <- function(fit_df, fit, us, us_pts=NULL) {
    loc_risk <- tibble(
        full=unique(fit_df$location),
        risk=fit$summary.fixed$mean[-1]
    )
    
    ret <- us |> 
        filter(!(full %in% iso_states)) |> 
        left_join(loc_risk, by=c("full")) |> 
        ggplot() +
        geom_sf(aes(fill=risk), col="black") +
        scale_fill_gradient2() +
        theme_map()
    
    if (!is.null(us_pts)) {
        us_pts <- left_join(us_pts, loc_risk, by=c("full"))
        
        ret <- ret +
            geom_sf(aes(fill=risk), us_pts, shape=21, size=2.4, col="white") +
            geom_sf_text(aes(label=abbr), us_pts, nudge_x=150000)
    }
    return(ret)
}

us <- us_map(regions="states")
iso_states <- c("Alaska", "Hawaii", "Puerto Rico")
small_states <- c("District of Columbia", "Rhode Island")

us_pts <- tribble(
    ~abbr, ~full, ~x, ~y,
    "AK", "Alaska", 2350000, -500000,
    "DC", "District of Columbia", 2350000, -700000,
    "RI", "Rhode Island", 2700000, -700000,
    "HI", "Hawaii", 2350000, -900000,
    "PR", "Puerto Rico", 2700000, -900000
) |> 
    st_as_sf(crs=st_crs(us), coords=c("x", "y"))

p1 <- plot_loc_effect(fit_df_flu, fit_flu, us, us_pts)
p2 <- plot_loc_effect(fit_df_rsv, fit_rsv, us)
p3 <- plot_loc_effect(fit_df_covid, fit_covid, us, us_pts)

#  option two for spatial main effects: simple mean + CIs-------------------------
summ_fx <- function(fx_marg) {
    imap_dfr(fx_marg, ~{
        ci5 <- inla.hpdmarginal(0.5, .x)
        ci95 <- inla.hpdmarginal(0.95, .x)
        tibble(
            var = .y,
            m = inla.zmarginal(.x, TRUE)$mean,
            lo5 = ci5[,1], hi5 = ci5[,2],
            lo95 = ci95[,1], hi95 = ci95[,2]
        )
    })
}

post_loc_flu <- fit_flu$marginals.fixed |> 
    summ_fx() |> 
    slice(2:n())

post_loc_rsv <- fit_rsv$marginals.fixed |> 
    summ_fx() |> 
    slice(2:n())

post_loc_covid <- fit_covid$marginals.fixed |> 
    summ_fx() |> 
    slice(2:n())

p1 <- bind_rows(
    mutate(post_loc_flu, disease="flu"),
    mutate(post_loc_rsv, disease="rsv"),
    mutate(post_loc_covid, disease="covid")
) |> 
    mutate(state=str_remove(var, "location"), disease=fct_inorder(disease)) |> 
    ggplot(aes(m, fct_rev(state), col=disease)) +
    geom_vline(xintercept=0, col='gray70', linetype="dashed", size=1.5) +
    geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth=1, alpha=0.6) +
    # geom_linerange(aes(xmin = lo5, xmax = hi5), linewidth=2) +
    geom_point(size=1.5) +
    facet_wrap(~disease, nrow=1) +
    labs(x="log relative risk", y=NULL) +
    theme_bw() +
    theme(legend.position="none")

seasonal_summ <- imap_dfr(list(flu=fit_flu, rsv=fit_rsv, covid=fit_covid), \(ft, d) {
    ft$summary.random$epiweek |> 
        as_tibble() |> 
        select(epiweek=ID, mean, contains("quant")) |> 
        mutate(disease=d)
})

p4 <- seasonal_summ |> 
    mutate(disease=fct_inorder(disease)) |> 
    ggplot(aes(epiweek, mean, col=disease, fill=disease)) +
    annotate(
        "segment", 
        x=0, y=0, xend=53, yend=0, 
        col="gray40", linetype="13", linewidth=1.1
    ) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col=NA, alpha=0.3) +
    geom_ribbon(aes(ymin=`0.25quant`, ymax=`0.75quant`), col=NA, alpha=0.5) +
    geom_line(alpha=0.8) +
    scale_x_continuous(n.breaks=4) +
    labs(x="epiweek", y="seasonal effect", col=NULL, fill=NULL) +
    theme_mod_output()

get_main_effect <- function(fit_df, ft) {
    ft$summary.random$t |> 
        as_tibble() |> 
        mutate(date=unique(fit_df$date)) |> 
        select(date, mean, contains("quant"))
}

get_interaction_effect <- function(fit_df, ft) {
    fit_df |> 
        arrange(location, date) |> 
        bind_cols(as_tibble(ft$summary.random$iloc)) |> 
        select(date, location, mean, contains("quant"))
}

short_main_summ <- bind_rows(
    mutate(get_main_effect(fit_df_flu, fit_flu), disease="flu"),
    mutate(get_main_effect(fit_df_rsv, fit_rsv), disease="rsv"),
    mutate(get_main_effect(fit_df_covid, fit_covid), disease="covid")
)

p5 <- short_main_summ |> 
    mutate(disease=fct_inorder(disease)) |> 
    filter(date > "2021-09-01", date <= max(flu$date)) |> 
    ggplot(aes(date, mean, col=disease, fill=disease)) +
    geom_hline(yintercept=0, col="gray40", linetype="13", linewidth=1.1) +
    # annotate(
    #     "segment", 
    #     min(short_main_summ$date), y=0, xend=max(short_main_summ$date), yend=0,
    #     col="gray40", linetype="13", linewidth=1.1
    # ) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col=NA, alpha=0.3) +
    geom_ribbon(aes(ymin=`0.25quant`, ymax=`0.75quant`), col=NA, alpha=0.5) +
    geom_line(alpha=0.8) +
    scale_x_date(date_breaks="3 months", date_labels="%b %y", guide=guide_axis(angle=45)) +
    labs(y="short-term main effect", col=NULL, fill=NULL) +
    theme_mod_output()

# put everything together - version 1
plot_grid(
    plot_grid(p1, p2, p3, nrow=1, labels="AUTO"),
    plot_grid(p4, p5, rel_widths=c(0.8, 1), labels=c("D", "E")),
    nrow=2
)

ggsave("figs/fig4-draft1.pdf", width=13, height=7.5)

# version 2

plot_grid(
    plot_grid(
        p1, 
        plot_grid(NULL, p4, NULL, rel_heights=c(0.15, 1, 0.15), nrow=3), 
        rel_widths=c(1, 0.9), labels=c("A", "B")
    ),
    plot_grid(p5, NULL, rel_widths=c(1, 0.2), nrow=1, labels=c("C")), 
    nrow=2, rel_heights=c(1, 0.95)
)

ggsave("figs/fig4-draft2.pdf", width=8, height=8)

###
inter_summ <- bind_rows(
    mutate(get_interaction_effect(fit_df_flu, fit_flu), disease="flu"),
    mutate(get_main_effect(fit_df_rsv, fit_rsv), disease="rsv"),
    mutate(get_main_effect(fit_df_covid, fit_covid), disease="covid")
)

###
fit_flu$misc$configs$contents

###
get_post_sd <- function(tau) {
    inla.tmarginal(\(x) sqrt(1/x), tau)
}

m1 <- get_post_sd(fit$marginals.hyperpar$`Precision for epiweek`)

m2 <- inla.tmarginal(
    \(x) sqrt(1/x),
    fit$marginals.hyperpar$`Precision for t`
)

max(m2[,1])

p1 <- ggplot() +
    stat_function(
        fun=~inla.dmarginal(., m1), 
        geom="area", fill="lightblue", col="gray40", alpha=0.7
    ) +
    xlim(0, 1.3) +
    labs(x="value", y="density", title="seasonal standard deviation") +
    theme_mod_output()

p2 <- ggplot() +
    stat_function(
        fun=~inla.dmarginal(., m2), 
        geom="area", fill="lightblue", col="gray40", alpha=0.7
    ) +
    xlim(0, 0.8) +
    labs(x="value", y="density", title="short-term standard deviation") +
    theme_mod_output()

p3 <- ggplot() +
    stat_function(
        fun=~inla.dmarginal(., fit$marginals.hyperpar$`Rho for t`), 
        geom="area", fill="lightblue", col="gray40", alpha=0.7
    ) +
    xlim(0, 1) +
    labs(x="value", y="density", title="correlation between days") +
    theme_mod_output()

p4 <- ggplot() +
    stat_function(
        fun=~inla.dmarginal(., fit$marginals.hyperpar$`GroupRho for t`), 
        geom="area", fill="lightblue", col="gray40", alpha=0.7
    ) +
    xlim(0, 1) +
    labs(x="value", y="density", title="correlation between age groups") +
    theme_mod_output()

plot_grid(p1, p2, p3, p4, nrow=2)