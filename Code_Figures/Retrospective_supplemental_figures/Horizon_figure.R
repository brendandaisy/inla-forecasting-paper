library(covidHubUtils)
library(scoringutils)
library(gt)
library(tidyverse)
library(patchwork)

#########
#read in data
#########
overall_df_flu <- read_csv("../data/overall_df_flu.csv")

overall_df_covid <- read_csv("../data/overall_df_covid.csv")

overall_df_rsv <- read_csv("../data/overall_df_rsv.csv")

simple <- read_csv("../data/simple_scored_horizon.csv")

simple_ <- simple |>
  select(model, horizon, r_wis, improvement, percent)

score_all <- simple_ %>%
  mutate(model = recode(model,
                        "COVID-19_forecast_model" = "COVID-19",
                        "Influenza_forecast_model" = "Influenza",
                        "RSV_forecast_model" = "RSV"))


###########
#Panel A
###########

#######
#function
########
make_rwis <- function(df, baseline_model) {
  df <- df %>%
    group_by(horizon, location, date, quantile) %>%
    mutate(
      baseline_wis = interval_score[model == baseline_model],  # extract baseline WIS
      baseline_wis = ifelse(is.na(baseline_wis) | baseline_wis == 0, 0.001, baseline_wis),  # replace 0 or NA baseline WIS with 0.001
      r_wis = ifelse(model == baseline_model, 1, interval_score / baseline_wis),  # calculate relative WIS
      improvement = ifelse(model == baseline_model, 0, (1 - r_wis) * 100)  # calculate improvement
    ) %>%
    ungroup()
  
  # Calculate the proportion of forecasts better than or equal to the baseline by horizon and disease
  df_summary <- df %>%
    group_by(horizon, model) %>%
    summarise(
      proportion_better_than_baseline = mean(r_wis <= 1, na.rm = TRUE)  # proportion of forecasts better than or equal to the baseline
    ) %>%
    ungroup()
  
  return(df_summary)
}


#########
#plot
########
library(cowplot)


custom_colors <- c(
  "COVID-19" = "#66C2A5",
  "Influenza" = "#8DA0CB",
  "RSV" = "#E78AC3"
)

# Create the plot
slide12 <- ggplot(score_all, aes(x = horizon, y = percent, color = model)) +
  geom_point(size = 3) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "azure4", size = 1) +
  labs(
    x = "Horizon",
    y = "Relative improvement (Baseline)") +
  theme_minimal() +
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  scale_color_manual(values = custom_colors) +  
  scale_x_continuous(breaks = unique(score_all$horizon)) +
  scale_y_continuous(labels = scales::percent) 

slide12 <- slide12 + background_grid(major = 'xy', minor = 'y')


print(slide12)

# Save the plot with white background using save_plot
#save_plot("../Figures/sup_hor_simple.png", slide12, base_width = 8, base_asp = 1.7, bg = "white")


#############
#Panel B
############
score_flu_all <- score(overall_df_flu) 
score_flu_all <- make_rwis(score_flu_all, "baseline_flu")

score_covid_all <- score(overall_df_covid) 
score_covid_all <- make_rwis(score_covid_all, "baseline_covid")

score_rsv_all <- score(overall_df_rsv) 
score_rsv_all <- make_rwis(score_rsv_all, "baseline_rsv")

scoredt2 <- rbind(score_covid_all, score_flu_all)

score_all_ht2 <- rbind(scoredt2, score_rsv_all)

score_all_ht2 <- score_all_ht2 |>
  filter(model %in% c("covid", "flu", "rsv_2018"))

score_all_ht2 <- score_all_ht2 %>%
  mutate(model = recode(model,
                        "covid" = "COVID-19",
                        "flu" = "Influenza",
                        "rsv_2018" = "RSV"))

#############
#plot
#############
library(cowplot)

# Define the custom colors for the models
custom_colors <- c(
  "COVID-19" = "#66C2A5",
  "Influenza" = "#8DA0CB",
  "RSV" = "#E78AC3"
)

# Create the plot with error bars
s12 <- ggplot(score_all_ht2, aes(x = horizon, y = proportion_better_than_baseline, color = model)) +
  geom_point(size = 3) +  # Plot the middle point
  labs(
    x = "Horizon",
    y = "Proportion Better or Equal (Baseline)") +
  theme_minimal() +
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()   
  ) +
  scale_color_manual(values = custom_colors) +  
  scale_x_continuous(breaks = unique(score_all_ht2$horizon)) #+
#scale_y_continuous(labels = scales::percent)

# Add grid lines for major and minor axes
s12 <- s12 + background_grid(major = 'xy', minor = 'y')

# Print the plot
print(s12)

# Save the plot with a white background using save_plot
#save_plot("../Figures/sup_proportion.png", s12, base_width = 10,  base_height = 6, bg = "white")

combined_plot <- slide12 + s12 + plot_layout(ncol = 2, widths = c(1, 1))  # Make plot2 wider
print(combined_plot)
#ggsave("../Figures/combined_plot_horizon.png", combined_plot, width = 14, height = 7)