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

max(flu$date)

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

fit_df_flu <- prep_fit_data(flu, weeks_ahead=0, ex_lam=population) |> 
    mutate(location=factor(location)) # necessary to include all 52 state intercepts!
fit_flu <- fit_inla_model(fit_df_flu, model, graph=graph, config=TRUE) # config was for trying conditional sampling

fit_df_rsv <- prep_fit_data(rsv, weeks_ahead=0, ex_lam=pop_served) |> 
    mutate(location=factor(location))

fit_rsv <- fit_inla_model(
    fit_df_rsv, 
    model_formula(seasonal="shared", temporal="ar1", spatial="exchangeable"),
    config=TRUE
)

fit_df_covid <- prep_fit_data(covid, weeks_ahead=0, ex_lam=population) |> 
    mutate(location=factor(location))
fit_covid <- fit_inla_model(fit_df_covid, model, graph=graph, config=TRUE)

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

us <- load_us_graph(flu) |> # for HHS region codes
    st_drop_geometry() |> 
    mutate(region=ifelse(region == 9, 3, region)) |> 
    mutate(region=factor(region, labels=c("Northeast", "Midwest", "South", "West/Pacific")))

post_loc <- bind_rows(
    mutate(post_loc_covid, disease="COVID-19"),
    mutate(post_loc_flu, disease="Influenza"),
    mutate(post_loc_rsv, disease="RSV")
) |> 
    mutate(state=str_remove(var, "location"), disease=fct_inorder(disease)) |> 
    left_join(us, by=c("state")) |> 
    # arrange(disease, m) |> 
    arrange(region) |> 
    mutate(state=fct_inorder(state))

# how many states were significantly different from zero for each disease?
post_loc |> 
    filter(disease == "COVID-19") |> 
    filter(lo95 > 0 | hi95 < 0)

library(grid)
.pt <- 1 / 0.352777778
len0_null <- function(x) {
    if (length(x) == 0)  NULL
    else                 x
}

theme_right_border <- function(colour = "black", size = 1, linetype = 1) {
    # use with e.g.: ggplot(...) + theme( panel.border=theme_right_border() ) + ...
    structure(
        list(colour = colour, size = size, linetype = linetype),
        class = c("theme_right_border", "element_blank", "element")
    )
}
element_grob.theme_right_border <- function(
        element, x = 0, y = 0, width = 1, height = 1,
        colour = NULL, size = NULL, linetype = NULL,
        ...) {
    gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
    element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
    polylineGrob(
        x = c(x+width, x+width), y = c(y,y+height), ..., default.units = "npc",
        gp = modifyList(element_gp, gp),
    )
}

p1 <- ggplot(post_loc, aes(m, fct_rev(state), col=m > 0)) +
    geom_vline(xintercept=0, col='gray70', size=0.93, linetype="11") +
    geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth=1, alpha=0.6) +
    # geom_linerange(aes(xmin = lo5, xmax = hi5), linewidth=2) +
    geom_point(size=1.5) +
    facet_grid(region~disease, scales="free_y", space="free_y", switch="y") +
    scale_color_manual(values=c("#AB76DE", "salmon2")) +
    # scale_color_gradient2(low="#ba6cb1", mid="gray70", high="#f0b924") +
    labs(x="Location effect", y=NULL) +
    theme_half_open() +
    # background_grid(major="y") +
    theme(
        legend.position="none",
        strip.background.y=theme_right_border(),
        strip.placement="outside",
        # strip.text.y=element_blank(),
        axis.text.y=element_text(size=rel(0.8))
        # plot.margin=unit(c(0, 0, 0, 0), "cm")
    )

#  Supp fig: association between coefficients for flu and covid-------------------
# library(gghighlight)
# 
# hl_states <- c("Vermont", "Puerto Rico", "Oklahoma", "Kentucky", "Washington", "Utah")
# 
# post_loc |>
#     filter(disease != "rsv") |>
#     pivot_wider(id_cols=c(state), names_from=disease, values_from=c(m, lo95, hi95)) |>
#     ggplot(aes(m_flu, m_covid, col=state)) +
#     geom_hline(yintercept=0, col="gray70", linetype="dashed") +
#     geom_vline(xintercept=0, col="gray70", linetype="dashed") +
#     geom_linerange(aes(xmin=lo95_flu, xmax=hi95_flu)) +
#     geom_linerange(aes(ymin=lo95_covid, ymax=hi95_covid)) +
#     geom_point() +
#     labs(x="Flu", y="COVID-19") +
#     theme_half_open() +
#     gghighlight(state %in% hl_states, max_highlight=Inf)
#     # theme(legend.position="bottom")
# 
# ggsave("figs/loc-coeffs-flu-covid.pdf", width=4.1, height=4)

#  version of above supplemental plot for spencer---------------------------------
post_loc_flu_c19 <- post_loc |>
    filter(disease != "RSV") |>
    pivot_wider(id_cols=c(state), names_from=disease, values_from=c(m, lo95, hi95))

post_loc_flu_rsv <- post_loc |>
    filter(disease != "COVID-19") |>
    pivot_wider(id_cols=c(state), names_from=disease, values_from=c(m, lo95, hi95)) |> 
    filter(!is.na(m_RSV))

post_loc_c19_rsv <- post_loc |>
    filter(disease != "Influenza") |>
    pivot_wider(id_cols=c(state), names_from=disease, values_from=c(m, lo95, hi95)) |> 
    filter(!is.na(m_RSV))

cor(post_loc_flu_c19$`m_COVID-19`, post_loc_flu_c19$m_Influenza)
cor(post_loc_flu_rsv$m_Influenza, post_loc_flu_rsv$m_RSV)
cor(post_loc_c19_rsv$`m_COVID-19`, post_loc_c19_rsv$m_RSV)

p1 <- post_loc |>
    filter(disease != "RSV") |>
    pivot_wider(id_cols=c(state), names_from=disease, values_from=c(m, lo95, hi95)) |>
    ggplot(aes(m_Influenza, `m_COVID-19`)) +
    geom_hline(yintercept=0, col="gray70", linetype="dashed") +
    geom_vline(xintercept=0, col="gray70", linetype="dashed") +
    geom_linerange(aes(xmin=lo95_Influenza, xmax=hi95_Influenza), col="gray50", alpha=0.5) +
    geom_linerange(aes(ymin=`lo95_COVID-19`, ymax=`hi95_COVID-19`), col="gray50", alpha=0.5) +
    geom_point(col="gray50") +
    geom_smooth(method="lm", se=FALSE, col="tomato", fullrange=TRUE) +
    annotate("text", label=as.expression(substitute(italic(r)^2~"="~0.174)), x=-0.88, y=0.72, col="tomato") +
    # coord_cartesian(xlim=c(-1.2, 1), ylim=c(-1, 1)) +
    labs(x="Influenza", y="COVID-19") +
    theme_half_open()

p2 <- post_loc_flu_rsv |>
    ggplot(aes(m_Influenza, m_RSV)) +
    geom_hline(yintercept=0, col="gray70", linetype="dashed") +
    geom_vline(xintercept=0, col="gray70", linetype="dashed") +
    geom_linerange(aes(xmin=lo95_Influenza, xmax=hi95_Influenza), col="gray50", alpha=0.5) +
    geom_linerange(aes(ymin=lo95_RSV, ymax=hi95_RSV), col="gray50", alpha=0.5) +
    geom_point(col="gray50") +
    geom_smooth(method="lm", se=FALSE, col="tomato", fullrange=TRUE) +
    annotate("text", label=as.expression(substitute(italic(r)^2~"="~0.37)), x=-1.08, y=1.09, col="tomato") +
    # coord_cartesian(xlim=c(-1.2, 1), ylim=c(-1, 1)) +
    labs(x="Influenza", y="RSV") +
    theme_half_open()

p3 <- post_loc_c19_rsv |>
    ggplot(aes(`m_COVID-19`, m_RSV)) +
    geom_hline(yintercept=0, col="gray70", linetype="dashed") +
    geom_vline(xintercept=0, col="gray70", linetype="dashed") +
    geom_linerange(aes(xmin=`lo95_COVID-19`, xmax=`hi95_COVID-19`), col="gray50", alpha=0.5) +
    geom_linerange(aes(ymin=lo95_RSV, ymax=hi95_RSV), col="gray50", alpha=0.5) +
    geom_point(col="gray50") +
    geom_smooth(method="lm", se=FALSE, col="tomato", fullrange=TRUE) +
    annotate("text", label=as.expression(substitute(italic(r)^2~"="~0.495)), x=-0.48, y=1.07, col="tomato") +
    # coord_cartesian(xlim=c(-1.2, 1), ylim=c(-1, 1)) +
    labs(x="COVID-19", y="RSV") +
    theme_half_open()

plot_grid(p1, p2, p3, nrow=1, rel_widths=c(1, 1, 1))
ggsave("figs/loc-coeff-corr.pdf", width=8.9, height=3.1)


#  Figure 4 drafts----------------------------------------------------------------
# ggplot(post_loc, aes(state, m, col=disease)) +
#     geom_hline(yintercept=0, col='gray70', size=1.3) +
#     geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth=1.08, alpha=0.6, position=position_dodge(width=1)) +
#     # geom_linerange(aes(xmin = lo5, xmax = hi5), linewidth=2) +
#     geom_point(size=2.1, position=position_dodge(width=1)) +
#     facet_grid(.~region, scales="free_x", space="free_x") +
#     labs(y="log relative risk", x=NULL) +
#     scale_x_discrete(guide=guide_axis(angle=40)) +
#     theme_half_open() +
#     background_grid(major="x") +
#     theme(
#         legend.position="none",
#         strip.background=element_blank(),
#         strip.text=element_blank(),
#         axis.text.x=element_text(size=rel(0.96)),
#         axis.title.y=element_text(size=rel(1.05))
#         # panel.spacing=unit(1, "pt")
#     )
# 
# ggsave("figs/fig4A-draft1.pdf", width=11.4, height=4.2)

# p1 <- ggplot(post_loc, aes(m, fct_rev(state), col=disease)) +
#     geom_vline(xintercept=0, col='gray70', linetype="dashed", size=1.5) +
#     geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth=1, alpha=0.6) +
#     # geom_linerange(aes(xmin = lo5, xmax = hi5), linewidth=2) +
#     geom_point(size=1.5) +
#     facet_wrap(~disease, nrow=1) +
#     labs(x="log relative risk", y=NULL) +
#     theme_bw() +
#     theme(legend.position="none")

seasonal_summ <- imap_dfr(list(covid=fit_covid, flu=fit_flu, rsv=fit_rsv), \(ft, d) {
    ft$summary.random$epiweek |> 
        as_tibble() |> 
        select(epiweek=ID, mean, contains("quant")) |> 
        mutate(disease=d)
})

seasonal_summ |> 
    group_by(disease) |> 
    filter(`0.025quant` > 0, epiweek >= 30) |> 
    slice_min(epiweek)

p4 <- seasonal_summ |> 
    mutate(disease=factor(disease, labels=c("COVID-19", "Influenza", "RSV"))) |> 
    ggplot(aes(epiweek, mean, col=disease, fill=disease)) +
    # annotate(
    #     "segment", 
    #     x=0, y=0, xend=53, yend=0, 
    #     col="gray40", linetype="13", linewidth=1.1
    # ) +
    geom_hline(yintercept=0, col="gray70", linetype="dashed", linewidth=1.03) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col=NA, alpha=0.3) +
    geom_ribbon(aes(ymin=`0.25quant`, ymax=`0.75quant`), col=NA, alpha=0.5) +
    geom_line(alpha=0.8, linewidth=1.03) +
    scale_x_continuous(breaks=seq(0, 50, 10)) +
    scale_color_manual(values=c("#7DC0A6", "#919FC7", "#DA8EC0")) +
    scale_fill_manual(values=c("#7DC0A6", "#919FC7", "#DA8EC0")) +
    labs(x="Epiweek", y="Seasonal effect", col=NULL, fill=NULL) +
    theme_half_open() +
    theme(legend.position="none")

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
    mutate(get_main_effect(fit_df_covid, fit_covid), disease="COVID-19"),
    mutate(get_main_effect(fit_df_flu, fit_flu), disease="Influenza"),
    mutate(get_main_effect(fit_df_rsv, fit_rsv), disease="RSV")
)

p5 <- short_main_summ |> 
    mutate(disease=fct_inorder(disease)) |> 
    filter(date > "2021-09-01", date <= max(flu$date)) |> 
    ggplot(aes(date, mean, col=disease, fill=disease)) +
    geom_hline(yintercept=0, col="gray70", linetype="dashed", linewidth=1.03) +
    # annotate(
    #     "segment", 
    #     min(short_main_summ$date), y=0, xend=max(short_main_summ$date), yend=0,
    #     col="gray40", linetype="13", linewidth=1.1
    # ) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), col=NA, alpha=0.3) +
    geom_ribbon(aes(ymin=`0.25quant`, ymax=`0.75quant`), col=NA, alpha=0.5) +
    geom_line(alpha=0.8, linewidth=1.03) +
    scale_x_date(date_breaks="3 months", date_labels="%b %y", guide=guide_axis(angle=45)) +
    scale_color_manual(values=c("#7DC0A6", "#919FC7", "#DA8EC0")) +
    scale_fill_manual(values=c("#7DC0A6", "#919FC7", "#DA8EC0")) +
    labs(y="Short-term main effect", col=NULL, fill=NULL, x=NULL) +
    theme_half_open() +
    theme(legend.position="none")

legend <- get_plot_component(
    p5 + theme(legend.position="right", legend.justification="center", legend.direction="vertical"),
    "guide-box-right",
    return_all=TRUE
)

# put everything together - version 1
# plot_grid(
#     plot_grid(p1, p2, p3, nrow=1, labels="AUTO"),
#     plot_grid(p4, p5, rel_widths=c(0.8, 1), labels=c("D", "E")),
#     nrow=2
# )
# 
# ggsave("figs/fig4-draft1.pdf", width=13, height=7.5)

# version 2

# plot_grid(
#     plot_grid(
#         p1, 
#         plot_grid(NULL, p4, NULL, rel_heights=c(0.15, 1, 0.15), nrow=3), 
#         rel_widths=c(1, 0.9), labels=c("A", "B")
#     ),
#     plot_grid(p5, NULL, rel_widths=c(1, 0.2), nrow=1, labels=c("C")), 
#     nrow=2, rel_heights=c(1, 0.95)
# )

plot_grid(
    plot_grid(p4, p5, nrow=2, labels=c("A", "B")),
    legend,
    p1,
    nrow=1, rel_widths=c(0.67, 0.23, 1), labels=c("", "", "C")
)

# ggdraw(plot_grid(pfinal, legend, nrow=1, rel_widths=c(1, 0.15)))

ggsave("figs/fig4-draft5.pdf", width=10.5, height=6.8)

# getting credible intervals for seasonal min and max-----------------------------
post_seasonal_ci <- function(post, tag) {
    phi <- inla.posterior.sample.eval("epiweek", post)
    phi_df <- post_samp_seasonal2df(phi)
    
    peak <- phi_df |> 
        group_by(name) |> 
        slice_max(value, n=1) |> 
        pull(epiweek)
    
    start <- phi_df |> 
        group_by(name) |> 
        slice_min(value, n=1) |> 
        pull(epiweek)
    
    cross <- phi_df |> 
        filter(epiweek > mean(start), value > 0) |> 
        group_by(name) |> 
        slice_min(value) |> 
        pull(epiweek)
    
    labs <- str_c(tag, "_", c("start", "crossover", "peak"))
    
    ret <- list(
        quantile(start, c(0.025, 0.5, 0.975)),
        quantile(cross, c(0.025, 0.5, 0.975)),
        (quantile((peak - 30) %% 52, c(0.025, 0.5, 0.975)) + 29) %% 52 + 1
    )
    names(ret) <- labs
    
    return(ret)
}

post_samp_seasonal2df <- function(phi) {
    phi[1:52,] |> 
        as_tibble() |> 
        mutate(epiweek=1:n()) |> 
        pivot_longer(contains("sample")) |> 
        mutate(name=fct_inorder(str_extract(name, "\\d+")))
}

post_flu <- inla.posterior.sample(5000, fit_flu)
post_rsv <- inla.posterior.sample(5000, fit_rsv)
post_covid <- inla.posterior.sample(5000, fit_covid)

post_seasonal_ci(post_flu, "flu")
post_seasonal_ci(post_rsv, "rsv")
post_seasonal_ci(post_covid, "covid")

alpha <- inla.posterior.sample.eval("t", post_flu)
delta <- inla.posterior.sample.eval("iloc", post_flu)

ggplot(phi_df, aes(epiweek, value, col=name, group=name)) +
    geom_line()

# what is the average location effect across US regions for each disease?
locations <- str_c("location", unique(flu$location))

regions <- load_us_graph(flu) |> # for HHS region codes
    st_drop_geometry() |> 
    mutate(region=ifelse(region == 9, 3, region))
# , idx=1:n()) |> 
#     group_by(region) |> 
#     group_map(~pull(., idx))

beta_flu <- inla.posterior.sample.eval(locations, post_flu)
beta_covid <- inla.posterior.sample.eval(locations, post_covid)

beta_reg_flu <- beta_flu |> 
    as_tibble() |> 
    bind_cols(regions) |> 
    pivot_longer(contains("sample")) |> 
    group_by(region, name) |> 
    summarize(beta_region=mean(value)) |> 
    group_by(region) |> 
    reframe(
        q=quantile(beta_region, c(0.025, 0.25, 0.5, 0.75, 0.975)), 
        level=c("lo95", "lo50", "median", "hi50", "hi95")
    ) |> 
    pivot_wider(names_from=level, values_from=q) |> 
    mutate(disease="Influenza")

beta_reg_covid <- beta_covid |> 
    as_tibble() |> 
    bind_cols(regions) |> 
    pivot_longer(contains("sample")) |> 
    group_by(region, name) |> 
    summarize(beta_region=mean(value)) |> 
    group_by(region) |> 
    reframe(
        q=quantile(beta_region, c(0.025, 0.25, 0.5, 0.75, 0.975)), 
        level=c("lo95", "lo50", "median", "hi50", "hi95")
    ) |> 
    pivot_wider(names_from=level, values_from=q) |> 
    mutate(disease="COVID-19")

beta_reg <- bind_rows(beta_reg_covid, beta_reg_flu) |> 
    mutate(region=factor(region, labels=c("Northeast", "Midwest", "South", "West/Pacific")))

ggplot(beta_reg, aes(median, region, col=disease)) +
    geom_vline(xintercept=0, col='gray70', size=0.93, linetype="11") +
    geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth=1, alpha=0.6) +
    geom_linerange(aes(xmin = lo50, xmax = hi50), linewidth=1.5, alpha=0.8) +
    geom_point(size=2) +
    facet_wrap(~disease) +
    scale_color_manual(values=c("#7DC0A6", "#919FC7")) +
    scale_x_continuous(breaks=c(-0.5, -0.25, 0, 0.25, 0.5)) +
    labs(x="Location effect regional average", y=NULL) +
    theme_half_open() +
    # background_grid(major="y") +
    theme(
        legend.position="none"
        # strip.background.y=element_blank(),
        # strip.text.y=element_blank(),
        # axis.text.y=element_text(size=rel(0.8))
        #   6plot.margin=unit(c(0, 0, 0, 0), "cm")
    )

ggsave("figs/loc-coeff-regional.pdf", width=5.6, height=4)

# compare to a simple "mean of means" approach
post_loc |> 
    group_by(region, disease) |> 
    summarize(mean_reg=mean(m)) |> 
    arrange(disease)

# graveyard-----------------------------------------------------------------------
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