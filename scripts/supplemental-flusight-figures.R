# --------------------------------------------------------------------------------
# supplemental figures evaluating performance during the 2023-24 challenge--------
# --------------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(DT)
library(cowplot)

###############
#PIC Plots
###############

##################################### #read in and clean data 
date_df <- read_csv("results/wis_season_by_model_date.csv")

truth_df <- read_csv("results/truth_df_flu.csv")

truth_US <- truth_df |>
  filter(location_name == "US") |>
  select(date, count) |>
  rename(reference_date = date)

PIC_date <- date_df |>
  filter(model %in% c("UGA_flucast-INFLAenza", "FluSight-baseline", "FluSight-ensemble"))

PIC_date <- PIC_date |>
  mutate(coverage_50 = (cov_50/100),
         coverage_95 = (cov_95/100))

#add the truth column 
truth_US_clean <- truth_US |>
  filter(reference_date >= "2023-10-14")

################################### #Plotting
# Custom colors for each model
model_colors <- c("FluSight-baseline" = "#2b8cbe", 
                  "FluSight-ensemble" = "coral2", 
                  "UGA_flucast-INFLAenza" = "darkseagreen")

custom_labels <- c(
  "FluSight-ensemble" = "Ensemble",
  "FluSight-baseline" = "Baseline",
  "UGA_flucast-INFLAenza" = "INFLAenza"
)

plot_counts <- ggplot(data = truth_US_clean, aes(x = reference_date, y = count)) +
  geom_line(color = "black") +
  labs(x = NULL, y = "Hospital admits") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
  theme_minimal() +
  theme_cowplot() +
  scale_y_continuous(labels = scales::comma)

plot_counts <- plot_counts + background_grid(major = 'xy')

# Generate individual plots
plot_cov_50 <- ggplot(data = PIC_date, aes(x = reference_date, y = coverage_50, color = model)) +
  # Add a horizontal line at y = 0
  geom_hline(yintercept = .50, linetype = "dashed", color = "azure4", size = 1) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors, labels = custom_labels) +
  labs(x = NULL, y = "50% Coverage", color = "Model") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
  #labs(y = "cov_50") +
  theme_minimal() +
  theme_cowplot() +
  scale_y_continuous(labels = scales::percent) +
  
  theme(
    legend.position = "right",
    legend.title = element_blank(),   # Remove the legend title
    #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  )

plot_cov_50 <- plot_cov_50 + background_grid(major = 'xy')

plot_cov_95 <- ggplot(data = PIC_date, aes(x = reference_date, y = coverage_95, color = model)) +
  # Add a horizontal line at y = 0
  geom_hline(yintercept = .95, linetype = "dashed", color = "azure4", size = 1) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors, labels = custom_labels) +
  labs(x = NULL, y = "95% Coverage", color = "Model") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
  #labs(y = "cov_95") +
  theme_minimal()+
  theme_cowplot() +
  scale_y_continuous(labels = scales::percent) +
  
  theme(
    legend.position = "right",
    legend.title = element_blank(),   # Remove the legend title
    #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  )

plot_cov_95 <- plot_cov_95 + background_grid(major = 'xy')

# Stack the plots
stacked_plots <- plot_grid(plot_counts, plot_cov_50, plot_cov_95, ncol = 1, align = "v", labels = "AUTO")

# Display the stacked plots
print(stacked_plots)

#ggsave("PIC_date_Flusight_counts.png", plot = stacked_plots, width = 10, height = 7, dpi = 300, bg = "white")


##################
#Location Plot
##################

########################## # Read in data and prep
wis_season_by_location_name <- read_csv("results/wis_season_by_model_location_name.csv")

PIC_location <- wis_season_by_location_name |>
  filter(model %in% c("UGA_flucast-INFLAenza", "FluSight-baseline", "FluSight-ensemble"))

PIC_location <- PIC_location |>
  mutate(coverage_50 = (cov_50/100),
         coverage_95 = (cov_95/100))

############################## # graph 

# Custom colors for each model
model_colors <- c("FluSight-baseline" = "#2b8cbe", 
                  "FluSight-ensemble" = "coral2", 
                  "UGA_flucast-INFLAenza" = "darkseagreen")

# Custom labels for key
custom_labels <- c(
  "FluSight-ensemble" = "Ensemble",
  "FluSight-baseline" = "Baseline",
  "UGA_flucast-INFLAenza" = "INFLAenza"
)

# Compute mean coverage for sorting
PIC_location <- PIC_location %>%
  group_by(location_name, model) %>%
  summarise(mean_cov_50 = mean(coverage_50, na.rm = TRUE),
            mean_cov_95 = mean(coverage_95, na.rm = TRUE), .groups = "drop")

# Get ordering based on coverage_95 for UGA_flucast-INFLAenza
uga_order <- PIC_location %>%
  filter(model == "UGA_flucast-INFLAenza") %>%
  arrange(mean_cov_95) %>%
  pull(location_name)

# Apply this ordering to all models
PIC_location <- PIC_location %>%
  mutate(location_name = factor(location_name, levels = uga_order))


# Compute the mean coverage for each model (for reference lines)
mean_coverage <- PIC_location %>%
  group_by(model) %>%
  summarise(mean_cov_50 = mean(mean_cov_50, na.rm = TRUE),
            mean_cov_95 = mean(mean_cov_95, na.rm = TRUE), .groups = "drop")

# Generate plot for Coverage 50%
plot_cov_50 <- ggplot(data = PIC_location, aes(x = location_name, y = mean_cov_50, color = model)) +
  geom_hline(yintercept = .50, linetype = "dashed", color = "azure4", size = 1) +  # Reference line
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors, labels = custom_labels) +
  labs(x = NULL, y = "50% Coverage", color = "Model") +
  theme_minimal() +
  theme_cowplot() +
  scale_y_continuous(labels=scales::label_percent()) +
  theme(
    #legend.position = NULL,  # Ensure the legend is included
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  background_grid(major = 'xy')

# Add mean coverage lines per model
for (model in names(model_colors)) {
  avg_cov_50 <- mean_coverage %>% filter(model == !!model) %>% pull(mean_cov_50)
  plot_cov_50 <- plot_cov_50 + 
    geom_hline(yintercept = avg_cov_50, color = model_colors[[model]], linetype = "solid", size = 0.8)
}

# Generate plot for Coverage 95%
plot_cov_95 <- ggplot(data = PIC_location, aes(x = location_name, y = mean_cov_95, color = model)) +
  geom_hline(yintercept = .95, linetype = "dashed", color = "azure4", size = 1) +  # Reference line
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors, labels = custom_labels) +
  labs(x = NULL, y = " 95% Coverage", color = "Model") +
  theme_minimal() +
  theme_cowplot() +
  scale_y_continuous(labels=scales::label_percent()) +
  theme(
    #legend.position = c(-5, 0.9),
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  background_grid(major = 'xy')

# Add mean coverage lines per model
for (model in names(model_colors)) {
  avg_cov_95 <- mean_coverage %>% filter(model == !!model) %>% pull(mean_cov_95)
  plot_cov_95 <- plot_cov_95 + 
    geom_hline(yintercept = avg_cov_95, color = model_colors[[model]], linetype = "solid", size = 0.8)
}

# Extract the shared legend from one plot
shared_legend <- get_legend(
  plot_cov_50 + theme(
    legend.position = "right",
    legend.justification = c(0.5, 0.55)  # Moves legend closer to center
  )
)

# Remove individual legends from plots
plot_cov_50 <- plot_cov_50 + theme(legend.position = "none")
plot_cov_95 <- plot_cov_95 + theme(legend.position = "none")

# Stack the two plots vertically
stacked_plots <- plot_grid(
  plot_cov_50, plot_cov_95, 
  ncol = 1, align = "v", labels = "AUTO"
)

final_plot <- plot_grid(
  stacked_plots, shared_legend, 
  ncol = 2, rel_widths = c(1.5, 0.2)  # Adjust width to give space for legend
)

# Display the final plot
print(final_plot)


##################
#Horizon Plot
##################

########################## # Read in data and prep
horizon_df <- read.csv("results/wis_season_by_model_horizon.csv")

PIC_horizon <- horizon_df |>
  filter(model %in% c("UGA_flucast-INFLAenza", "FluSight-baseline", "FluSight-ensemble"))

PIC_horizon <- PIC_horizon |>
  mutate(coverage_50 = (cov_50/100),
         coverage_95 = (cov_95/100))

PIC_horizon <- PIC_horizon %>%
  mutate(horizon = horizon + 1)

############################## # graph 
# Custom colors for each model
model_colors <- c("FluSight-baseline" = "#2b8cbe", 
                  "FluSight-ensemble" = "coral2", 
                  "UGA_flucast-INFLAenza" = "darkseagreen")

# Custom labels for key
custom_labels <- c(
  "FluSight-ensemble" = "Ensemble",
  "FluSight-baseline" = "Baseline",
  "UGA_flucast-INFLAenza" = "INFLAenza"
)


# Generate individual plots
plot_cov_50 <- ggplot(data = PIC_horizon, aes(x = horizon, y = coverage_50, color = model)) +
  # Add a horizontal line at y = 0
  geom_hline(yintercept = .50, linetype = "dashed", color = "azure4", size = 1) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors, labels = custom_labels) +
  labs(x = "Horizon", y = "50% Coverage", color = "Model") +
  #scale_x_date(date_breaks = "1 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
  #labs(y = "cov_50") +
  theme_minimal() +
  theme_cowplot() +
  scale_y_continuous(labels=scales::label_percent()) +
  
  theme(
    legend.position = "right",
    legend.title = element_blank(),   # Remove the legend title
    axis.text.x = element_text(angle = 00, hjust = 1, vjust = 0.5) 
  )

plot_cov_50 <- plot_cov_50 + background_grid(major = 'xy')

plot_cov_95 <- ggplot(data = PIC_horizon, aes(x = horizon, y = coverage_95, color = model)) +
  # Add a horizontal line at y = 0
  geom_hline(yintercept = .95, linetype = "dashed", color = "azure4", size = 1) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors, labels = custom_labels) +
  labs(x = "Horizon", y = "95% Coverage", color = "Model") +
  #scale_x_date(date_breaks = "1 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
  #labs(y = "cov_95") +
  theme_minimal()+
  theme_cowplot() +
  scale_y_continuous(labels = scales::label_percent()) +
  
  theme(
    legend.position = "right",
    legend.title = element_blank(),   # Remove the legend title
    axis.text.x = element_text(angle = 00, hjust = 1, vjust = 0.5) 
  )

plot_cov_95 <- plot_cov_95 + background_grid(major = 'xy')

# Stack the plots
#stacked_plots <- plot_grid(plot_cov_50, plot_cov_95, ncol = 1, align = "v", labels = "AUTO")

# Display the stacked plots
#print(stacked_plots)

# Extract the shared legend from one plot
shared_legend <- get_legend(
  plot_cov_50 + theme(
    legend.position = "right",
    legend.justification = c(0.5, 0.55)  # Moves legend closer to center
  )
)

# Remove individual legends from plots
plot_cov_50 <- plot_cov_50 + theme(legend.position = "none")
plot_cov_95 <- plot_cov_95 + theme(legend.position = "none")

# Stack the two plots vertically
stacked_plots <- plot_grid(
  plot_cov_50, plot_cov_95, 
  ncol = 1, align = "v", labels = "AUTO"
)

final_plot <- plot_grid(
  stacked_plots, shared_legend, 
  ncol = 2, rel_widths = c(1.5, 0.2)  # Adjust width to give space for legend
)

# Display the final plot
print(final_plot)

##############
#Horizon Plot 2
##############

############################# # function needed
make_improv <- function(df, model1, model2) {
  df <- df %>%
    group_by(horizon) %>%
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

############################# #read in and clean data
truth_US <- truth_df |>
  filter(location_name == "US") |>
  select(date, count) |>
  rename(reference_date = date)

#################################

# Step 2: Filter for UGA_flucast_INFLAenza model
loc_df_baseline <- horizon_df |>
  filter(model %in% c("UGA_flucast-INFLAenza", "FluSight-baseline"))

loc_df_ensemble <- horizon_df |>
  filter(model %in% c("UGA_flucast-INFLAenza", "FluSight-ensemble"))

#Use function 
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
  select (model, horizon, percent_ensemble)

uga_only <- baseline |>
  left_join(for_merge, by = c("model", "horizon"))

#################################### # plot 
custom_colors <- c(
  "percent_ensemble" = "coral2",
  "percent_baseline" = "#2b8cbe"
)

custom_labels <- c(
  "percent_ensemble" = "Ensemble",
  "percent_baseline" = "Baseline"
)

# Reshape the data to long format for easier plotting
df_long <- uga_only %>%
  pivot_longer(cols = c(percent_baseline, percent_ensemble), 
               names_to = "type", 
               values_to = "percent")


# Create the plot
plot <- ggplot(df_long, aes(x = horizon, y = percent, color = type, group = type)) +
  
  # Add a horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "azure4", size = 1) +
  
  # Add a line connecting the points for each type
  geom_line(size = 1.2) +  # Add colored lines
  
  # Add points for percent_baseline and percent_ensemble
  geom_point(size = 2.5) +
  
  # Format the y-axis to show percentages
  scale_y_continuous(labels = scales::percent) +
  
  # Apply custom colors
  scale_color_manual(values = custom_colors, labels = custom_labels) +
  
  labs(x = "Horizon", y = "Relative Improvement", color = "Type") +
  
  # Use a clean minimal theme
  theme_minimal() +
  theme_cowplot() +
  
  theme(
    legend.position = "right",
    legend.title = element_blank(),   # Remove the legend title
    #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.x = element_text(angle = 00, hjust = 1, vjust = 0.5, size = 14), # Keep x-axis labels small
    axis.text.y = element_text(size = 14),     # Larger y-axis text
    axis.title.y = element_text(size = 16)#,    # Larger y-axis title
  )

# Adding background grids to the top panel plot
plot <- plot + background_grid(major = 'xy')

# Display the plot
plot
