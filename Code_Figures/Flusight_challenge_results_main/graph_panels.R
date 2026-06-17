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


#################
#Panel A
#################

#data prep

date_df <- read_csv("../data/wis_season_by_model_date.csv")

#function
make_improv <- function(df, model1, model2) {
  df <- df %>%
    group_by(reference_date) %>%
    mutate(
      rel_wis_model1 = ifelse(model == model1, rel_wis, NA),  # Extract rel_wis for model1
      rel_wis_model2 = ifelse(model == model2, rel_wis, NA)   # Extract rel_wis for model2
    ) %>%
    fill(rel_wis_model1, rel_wis_model2, .direction = "downup") %>%  # Ensure both values are available in all rows
    mutate(
      improvement = (1 - rel_wis_model2 / rel_wis_model1) * 100  # Calculate relative improvement between the two models
    ) %>%
    ungroup()
  
  return(df)
}

# Step 2: Filter for UGA_flucast_INFLAenza model
loc_df_baseline <- date_df |>
  filter(model %in% c("UGA_flucast-INFLAenza", "FluSight-baseline"))

loc_df_ensemble <-date_df |>
  filter(model %in% c("UGA_flucast-INFLAenza", "FluSight-ensemble"))

baseline <- make_improv(loc_df_baseline, "FluSight-baseline", "UGA_flucast-INFLAenza")

baseline <- baseline |>
  mutate(percent_baseline = improvement / 100)

#ensemble prep
ensemble <- make_improv(loc_df_ensemble, "FluSight-ensemble", "UGA_flucast-INFLAenza")

ensemble <- ensemble |>
  mutate(percent_ensemble = improvement / 100)


baseline <- baseline |>
  filter(model == "UGA_flucast-INFLAenza")

for_merge <- ensemble |>
  filter(model == "UGA_flucast-INFLAenza") |>
  select (model, reference_date, percent_ensemble)

uga_only_date <- baseline |>
  left_join(for_merge, by = c("model", "reference_date"))

truth_df <- read_csv("../data/truth_df_flu.csv")

truth_US <- truth_df |>
  filter(location_name == "US") |>
  select(date, count) |>
  rename(reference_date = date)


####
#plot
####
library(ggplot2)
library(dplyr)

custom_colors <- c(
  "percent_ensemble" = "coral2",
  "percent_baseline" = "#2b8cbe"
)

# Reshape the data to long format for easier plotting
df_long <- uga_only_date %>%
  pivot_longer(cols = c(percent_baseline, percent_ensemble), 
               names_to = "type", 
               values_to = "percent")

df_long <- df_long |>
  left_join(truth_US, by = "reference_date")

#merge the truth_df to get column count  
# Define the range of rel_wis and count
rel_wis_range <- range(df_long$percent)
count_range <- range(df_long$count)

# Function to rescale the count values to match the range of rel_wis
rescale_count <- function(x) {
  (x - count_range[1]) / (count_range[2] - count_range[1]) * (rel_wis_range[2] - rel_wis_range[1]) + rel_wis_range[1]
}

# Create the plot
plot <- ggplot(df_long, aes(x = reference_date, y = percent, color = type)) +
  
  # Add the rescaled count values with a secondary y-axis, plotting behind the dots and segments
  geom_line(aes(y = rescale_count(count)), color = "azure3", size = 0.7) +
  
  # Add a horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "azure4", size = 1) +
  
  # Add points for percent_baseline and percent_ensemble
  geom_point(size = 2) +
  
  # Add a line connecting the points for each location
  geom_line(aes(group = reference_date), color = "gray45") +
  
  # Format the y-axis to show percentages
  scale_y_continuous(labels = scales::percent) +
  
  # Apply custom colors
  scale_color_manual(values = custom_colors) +
  
  labs(x = NULL, y = "Realative Improvement", color = "Type") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) + 
  
  # Use a clean minimal theme
  theme_minimal() +
  theme_cowplot() +
  
  theme(
    legend.position = "none",
    legend.title = element_blank(),   # Remove the legend title
    #axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  )

# Adding background grids to the top panel plot
plot <- plot + background_grid(major = 'xy')

# Display the plot
plot

#ggsave("../Figures/panelB.3_wide.png", plot = plot, width = 10, height = 6, dpi = 300, bg = "white") 


#############
#Panels B
#############
wis_season_by_location_name <- read_csv("../data/wis_season_by_model_location_name.csv")

#function
make_improv <- function(df, model1, model2) {
  df <- df %>%
    group_by(location_name) %>%
    mutate(
      rel_wis_model1 = ifelse(model == model1, rel_wis, NA),  # Extract rel_wis for model1
      rel_wis_model2 = ifelse(model == model2, rel_wis, NA)   # Extract rel_wis for model2
    ) %>%
    fill(rel_wis_model1, rel_wis_model2, .direction = "downup") %>%  # Ensure both values are available in all rows
    mutate(
      improvement = (1 - rel_wis_model2 / rel_wis_model1) * 100  # Calculate relative improvement between the two models
    ) %>%
    ungroup()
  
  return(df)
}

# Step 2: Filter for UGA_flucast_INFLAenza model
loc_df_baseline <- wis_season_by_location_name |>
  filter(model %in% c("UGA_flucast-INFLAenza", "FluSight-baseline"))

loc_df_ensemble <- wis_season_by_location_name |>
  filter(model %in% c("UGA_flucast-INFLAenza", "FluSight-ensemble"))

baseline <- make_improv(loc_df_baseline, "FluSight-baseline", "UGA_flucast-INFLAenza")

baseline <- baseline |>
  mutate(percent_baseline = improvement / 100)

ensemble <- make_improv(loc_df_ensemble, "FluSight-ensemble", "UGA_flucast-INFLAenza")

ensemble <- ensemble |>
  mutate(percent_ensemble = improvement / 100)


baseline <- baseline |>
  filter(model == "UGA_flucast-INFLAenza")

for_merge <- ensemble |>
  filter(model == "UGA_flucast-INFLAenza") |>
  select (model, location_name, percent_ensemble)

uga_only <- baseline |>
  left_join(for_merge, by = c("model", "location_name"))

######
#plot
######
library(ggplot2)
library(dplyr)

custom_colors <- c(
  "percent_ensemble" = "coral2",
  "percent_baseline" = "#2b8cbe"
)

# Reshape the data to long format for easier plotting
df_long <- uga_only %>%
  pivot_longer(cols = c(percent_baseline, percent_ensemble), 
               names_to = "type", 
               values_to = "percent")

# Reorder the location_name factor levels based on the percent_baseline values (descending)
df_long <- df_long %>%
  mutate(location_name = factor(location_name, 
                                levels = uga_only %>%
                                  arrange(desc(percent_baseline)) %>%
                                  pull(location_name)))

# Create the plot
plot <- ggplot(df_long, aes(x = location_name, y = percent, color = type)) +
  
  # Add a horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "azure4", size = 1) +
  
  # Add points for percent_baseline and percent_ensemble
  geom_point(size = 2) +
  
  # Add a line connecting the points for each location
  geom_line(aes(group = location_name), color = "gray45") +
  
  # Format the y-axis to show percentages
  scale_y_continuous(labels = scales::percent) +
  
  # Apply custom colors
  scale_color_manual(values = custom_colors) +
  
  labs(x = NULL, y = "Realative Improvement", color = "Type") +
  
  # Use a clean minimal theme
  theme_minimal() +
  theme_cowplot() +
  
  theme(
    legend.position = "none",
    legend.title = element_blank(),   # Remove the legend title
    #axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  )

# Adding background grids to the top panel plot
plot <- plot + background_grid(major = 'xy')

# Display the plot
plot

#ggsave("../Figures/panelB.3_average_wide.png", plot = plot, width = 10, height = 6, dpi = 300, bg = "white") 