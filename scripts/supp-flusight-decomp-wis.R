library(tidyverse)
library(INLA)
library(sf)
library(spdep)
library(lubridate)
library(covidHubUtils)
library(scoringutils)
library (ggplot2)
library(readr)
library(DT)
library(cowplot)

scores_horizon <- read.csv("results/wis_season_by_model_horizon.csv")

scores_h_small <- scores_horizon |>
    filter(model %in% c("FluSight-ensemble", "UGA_flucast-INFLAenza")) %>%
    mutate(model = recode(model,
                          "UGA_flucast-INFLAenza" = "INFLAenza",
                          "FluSight-ensemble" = "Ensemble"
    ))

scores_h_larger <- scores_horizon %>%
    filter(
        horizon == 3,
        model %in% c("FluSight-ensemble", "UGA_flucast-INFLAenza", "UMass-flusion", "MIGHTE-Nsemble", "PSI-PROF_beta")
    ) %>%
    mutate(model = recode(model,
                          "UMass-flusion"     = "E1",
                          "MIGHTE-Nsemble"    = "E3",
                          "PSI-PROF_beta"     = "I2",
                          "UGA_flucast-INFLAenza" = "INFLAenza",
                          "FluSight-ensemble" = "Ensemble"
    ))

df_long <- scores_h_small %>%
    pivot_longer(
        cols = c("underprediction", "dispersion", "overprediction"),
        names_to = "component",
        values_to = "value"
    )

df_long <- df_long |>
    mutate(horizon = horizon + 1)

all_horizon <- ggplot(df_long, aes(x = model, y = value, fill = component)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ horizon, nrow = 1, labeller = labeller(horizon = function(x) paste0("H=", x))) +
    labs(
        x = NULL,
        y = "WIS",
        fill = "Prediction error"#,
        #title = "Decomposition of WIS by Model and Horizon"
    ) +
    scale_fill_manual(
        values = c(
            "underprediction" = "#377eb8",
            "dispersion" = "#984ea3",
            "overprediction" = "#ff7f00"
        )
    ) +
    #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_minimal() +
    theme(
        legend.position  = "right",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.spacing.x = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.text = element_text(size = 8)
    )

print(all_horizon)

##
scores_h_larger <- scores_h_larger %>%
    mutate(model = factor(model, levels = c("E1", "Ensemble", "E3", "INFLAenza", "I2")))


df_long <- scores_h_larger %>%
    pivot_longer(
        cols = c("underprediction", "dispersion", "overprediction"),
        names_to = "component",
        values_to = "value"
    )

df_long <- df_long |>
    mutate(horizon = horizon + 1)

top_model <- ggplot(df_long, aes(x = model, y = value, fill = component)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ horizon, nrow = 1, labeller = labeller(horizon = function(x) paste0("H=", x))) +
    labs(
        x = NULL,
        y = "WIS",
        fill = "Prediction error"#,
        #title = "Decomposition of WIS by Model and Horizon"
    ) +
    scale_fill_manual(
        values = c(
            "underprediction" = "#377eb8",
            "dispersion" = "#984ea3",
            "overprediction" = "#ff7f00"
        )
    ) +
    theme_minimal() +
    theme(
        legend.position  = "right",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.spacing.x = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.text = element_text(size = 8)
    )

print(top_model)

plot_grid(all_horizon, top_model, ncol = 1, align = "v", labels = "auto", rel_heights = c(1, 1))

# Save the stacked plot as a PNG file

ggsave(
    "figs/WIS_components_plots_changes.png",
    width = 8.8, 
    height = 5.5, 
    units = "in", 
    bg = "white")
