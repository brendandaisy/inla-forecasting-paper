library(tidyverse)
library(lubridate)
library(cowplot)

rsv <- read_csv("data/weekly-rsv-us.csv") |> 
    filter(!(date %within% interval(ymd("2020-03-01"), ymd("2021-06-01")))) |> 
    mutate(season_week=(epiweek-40) %% 52 + 1) # epiweek 40 is start of RSV season, going off provided season column

flu <- read_csv("data/weekly-flu-us.csv") |> 
    filter(
        location != "US", # won't need nat'l for these figures
        date > "2021-09-01" # No real point in looking at data before around here
    ) |>
    filter(date >= "2021-10-07") |> # also go ahead and remove pre season dates for this fig
    mutate(
        season=case_when(
            date <= "2022-10-01" ~ "2021-22",
            date > "2022-10-01" & date <= "2023-09-30" ~ "2022-23",
            date > "2023-09-30" ~ "2023-24"
        ),
        season_week=(epiweek-40) %% 52 + 1
    )

flu_nat <- flu |> 
    group_by(epiweek, season_week) |> 
    summarise(mean=mean(weekly_rate), .groups="drop")

flu_pr <- flu |>
    filter(location == "Puerto Rico", season %in% c("2022-23", "2023-24"))
    # group_by(epiweek, season_week) |> 
    # summarise(mean=mean(weekly_rate), .groups="drop")

p1 <- ggplot(flu_nat, aes(season_week, mean)) +
    geom_line(
        aes(season_week, weekly_rate, group=interaction(season, location)), 
        filter(flu, location != "Puerto Rico", season %in% c("2022-23", "2023-24")), 
        col="gray50", alpha=0.1
    ) +
    geom_line(col="black", linewidth=1.1) +
    geom_line(aes(season_week, weekly_rate, group=interaction(season, location)), rsv_ga, col="#AB76DE", linewidth=1.1) +
    annotate("text", label="Puerto Rico", x=42, y=2, col="#AB76DE", size=5) +
    labs(x="Respiratory season week", y="Flu Admits per 100k") +
    scale_x_continuous(breaks=seq(0, 50 , 5)) +
    coord_cartesian(ylim=c(0, 8)) +
    theme_half_open()

rsv_nat <- rsv |> 
    group_by(epiweek, season_week) |> 
    summarise(mean=mean(weekly_rate), .groups="drop")

rsv_ga <- rsv |>
    filter(location == "Georgia", season %in% c("2022-23", "2023-24"))
    # group_by(epiweek, season_week) |> 
    # summarise(mean=mean(weekly_rate), .groups="drop")

p2 <- ggplot(rsv_nat, aes(season_week, mean)) +
    geom_line(
        aes(season_week, weekly_rate, group=interaction(season, location)), 
        filter(rsv, location != "Georgia", season %in% c("2022-23", "2023-24")), 
        col="gray50", alpha=0.1
    ) +
    geom_line(col="black", linewidth=1.1) +
    geom_line(aes(season_week, weekly_rate, group=interaction(season, location)), rsv_ga, col="salmon2", linewidth=1.1) +
    annotate("text", label="Georgia", x=44, y=2, col="salmon2", size=5) +
    labs(x="Respiratory season week", y="RSV Admits per 100k") +
    scale_x_continuous(breaks=seq(0, 50 , 5)) +
    coord_cartesian(ylim=c(0, 5)) +
    theme_half_open()

plot_grid(p1, p2, labels="AUTO")
ggsave("figs/fig-S5-v3.pdf", width=9.2, height=4.5)
