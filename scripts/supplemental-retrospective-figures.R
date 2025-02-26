# --------------------------------------------------------------------------------
# supplemental figures summarizing performance in the retrospective analysis------
# --------------------------------------------------------------------------------
library(tidyverse)

###########################
#Relative Performance Plots
############################

######################## #read in data and basic prep
score_all_d <- read_csv("results/date_score.csv")

score_all_d <- score_all_d |>
  mutate(average_percent = mean(percent, na.rm = TRUE))

truth_df_all <- read_csv("results/truth_df_fig3.csv")

truth_df <- truth_df_all |>
  filter(location == "US") |>
  filter(date >= ymd("2022-08-30"),
         date <= ymd("2024-04-27"))

truth_df <- truth_df %>%
  mutate(disease = recode(disease,
                          "covid" = "COVID-19 Historical Data",
                          "flu" = "Influenza Historical Data",
                          "rsv" = "RSV Historical Data"))

truth_df2 <- truth_df |>
  filter(date >= ymd("2022-09-17"))

######################### #Plot

####top panel
library(cowplot)

# Define the custom colors
custom_colors <- c(
  "COVID-19 Historical Data" = "#66C2A5",
  "Influenza Historical Data" = "#8DA0CB",
  "RSV Historical Data" = "#E78AC3"
)

# Create an empty list to store the plots
plot_list_pt <- list()

# Loop through each unique model and generate a separate plot for each
disease <- unique(truth_df2$disease)

for (dis in disease) {
  # Filter the data for each model
  score_model <- subset(truth_df2, disease == dis)
  
  # Create the plot for each model
  plot <- ggplot(score_model, aes(x = date, y = count, color = disease, group = disease)) +
    geom_line(size = 1.2) +  
    #geom_hline(yintercept = 0, linetype = "dashed", color = "azure4", size = 1) +
    labs(
      #title = dis,  # Use model name as title
      x = NULL,
      y = "Hospital admits") +
    theme_minimal() +
    theme_cowplot() +
    theme(
      #plot.title = element_text(hjust = 0.5),
      legend.position = "none",  # No need for legend in individual plots
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
      # Add background grid for major x and y axes
      #panel.grid.major.x = element_line(color = "grey80"),
      #panel.grid.major.y = element_line(color = "grey80"),
      #panel.grid.minor = element_blank()  # Remove minor grid lines
    ) +
    scale_color_manual(values = custom_colors) +  
    #scale_x_date(date_breaks = "3 months", date_labels = "%b '%y", guide = guide_axis(angle = 0), expand = c(0.02, 0.02)) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
    scale_y_continuous(labels = scales::comma)
  
  # Adding background grids to the top panel plot
  plot <- plot + background_grid(major = 'xy')
  
  # Save the plot in the list using the model name as the key
  plot_list_pt[[dis]] <- plot
}


print(plot_list_pt[["COVID-19 Historical Data"]])
print(plot_list_pt[["Influenza Historical Data"]])
print(plot_list_pt[["RSV Historical Data"]])

######## #Bottom panel

# Define the custom colors
custom_colors <- c(
  "COVID-19_forecast_model" = "#66C2A5",
  "Influenza_forecast_model" = "#8DA0CB",
  "RSV_forecast_model" = "#E78AC3"
  #"COVID-19 Model" = "#66C2A5",
  #"Influenza Model" = "#8DA0CB",
  #"RSV Model" = "#E78AC3"
)


plot_list_pb <- list()

models <- unique(score_all_d$model)

for (mod in models) {
  # Filter the data for each model
  score_model <- subset(score_all_d, model == mod)
  
  plot <- ggplot(score_model, aes(x = date, y = percent, color = model, group = model)) +
    geom_point(size = 2) + 
    geom_line(aes(y = average_percent, color = model, group = model), size = 0.75) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "azure4", size = 1) +
    labs(
      #title = mod,  # Use model name as title
      x = NULL ,
      y = "Relative improvement") +
    theme_minimal() +
    theme_cowplot() +
    theme(
      # plot.title = element_text(hjust = 0.5),
      legend.position = "none",  
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
      #panel.grid.major.x = element_line(color = "grey80"),
      #panel.grid.major.y = element_line(color = "grey80"),
      #panel.grid.minor = element_blank()  
    ) +
    scale_color_manual(values = custom_colors) +  
    #scale_x_date(date_breaks = "3 months", date_labels = "%b '%y", guide = guide_axis(angle = 0), expand = expansion(mult = c(0.01, 0.01))) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
    scale_y_continuous(labels = scales::percent)
  
  plot <- plot + background_grid(major = 'xy')
  
  plot_list_pb[[mod]] <- plot
}

print(plot_list_pb[["COVID-19_forecast_model"]])
print(plot_list_pb[["Influenza_forecast_model"]])
print(plot_list_pb[["RSV_forecast_model"]])

############################# #stack and save

# Create a vector of matching diseases/models
diseases_truth <- c("COVID-19 Historical Data", "Influenza Historical Data", "RSV Historical Data")
#models_forecast <- c("COVID-19 Model", "Influenza Model", "RSV Model")
models_forecast <- c("COVID-19_forecast_model", "Influenza_forecast_model", "RSV_forecast_model")

# Loop through the diseases and corresponding models
for (i in seq_along(diseases_truth)) {
  
  # Get the truth plot and forecast plot for the current disease/model
  truth_plot <- plot_list_pt[[diseases_truth[i]]]
  forecast_plot <- plot_list_pb[[models_forecast[i]]]
  
  # Stack the truth plot on top of the forecast plot
  stacked_plot <- plot_grid(truth_plot, forecast_plot, ncol = 1, align = "v", labels = "AUTO", rel_heights = c(1, 1))
  print(stacked_plot)
  
  # Save the stacked plot as a PNG file
  # file_name <- paste0(i, "improvement_retro_date.png")
  # ggsave(filename = file_name, plot = stacked_plot, path = "/Users/maya/Documents/Figures/Code_Figures/Retrospective_supplemental_figures", width = 8.8, height = 5.5, units = "in", bg = "white")
}


###################
#PIC Plot date
###################

##################### # Top Panel graph

# Define the custom colors
custom_colors <- c(
  "COVID-19 Historical Data" = "#66C2A5",
  "Influenza Historical Data" = "#8DA0CB",
  "RSV Historical Data" = "#E78AC3"
)

# Create an empty list to store the plots
plot_list_pt <- list()

# Loop through each unique model and generate a separate plot for each
disease <- unique(truth_df2$disease)

for (dis in disease) {
  # Filter the data for each model
  score_model <- subset(truth_df2, disease == dis)
  
  # Create the plot for each model
  plot <- ggplot(score_model, aes(x = date, y = count, color = disease, group = disease)) +
    geom_line(size = 1.2) +  
    labs(
      x = NULL,
      y = "Hospital admits") +
    theme_minimal() +
    theme_cowplot() +
    theme(
      legend.position = "none",  # No need for legend in individual plots
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    scale_color_manual(values = custom_colors) +  
    scale_x_date(date_breaks = "3 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
    scale_y_continuous(labels = scales::comma)
  
  # Adding background grids to the top panel plot
  plot <- plot + background_grid(major = 'xy')
  
  # Save the plot in the list using the model name as the key
  plot_list_pt[[dis]] <- plot
}


print(plot_list_pt[["COVID-19 Historical Data"]])
print(plot_list_pt[["Influenza Historical Data"]])
print(plot_list_pt[["RSV Historical Data"]])


################################# #bottom panel data prep and graph
score_all_d <- score_all_d |>
  filter(date <= "2024-04-27")

score_all_d <- score_all_d %>%
  group_by(model) %>%
  complete(date = seq(min(date), max(date), by = "1 week")) %>%
  fill(model, .direction = "down") %>% # Ensure model column is filled in each group
  ungroup()

# Custom colors for each model
model_colors <- c("COVID-19_forecast_model" = "#66C2A5", 
                  "Influenza_forecast_model" = "#8DA0CB", 
                  "RSV_forecast_model" = "#E78AC3")

# Initialize lists to store the individual plots
plot_list_cov_50 <- list()
plot_list_cov_95 <- list()
plot_list_truth <- list()

# Define the diseases/models for plotting
diseases_truth <- c("COVID-19 Historical Data", "Influenza Historical Data", "RSV Historical Data")
models_forecast <- c("COVID-19_forecast_model", "Influenza_forecast_model", "RSV_forecast_model")

# Loop through each model to create individual plots for cov_50, cov_95, and truth
for (i in seq_along(models_forecast)) {
  
  # Filter data for each model
  model_data <- score_all_d %>% filter(model == models_forecast[i])
  
  # Plot cov_50
  plot_list_cov_50[[models_forecast[i]]] <- ggplot(data = model_data, aes(x = date, y = coverage_50, color = model)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "azure4", size = 1) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = model_colors) +
    labs(x = NULL, y = "50% Coverage") +
    scale_x_date(date_breaks = "3 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
    theme_minimal() +
    theme_cowplot() +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "none", legend.title = element_blank()) +
    background_grid(major = 'xy')
  
  # Plot cov_95
  plot_list_cov_95[[models_forecast[i]]] <- ggplot(data = model_data, aes(x = date, y = coverage_95, color = model)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "azure4", size = 1) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = model_colors) +
    labs(x = NULL, y = "95% Coverage") +
    scale_x_date(date_breaks = "3 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
    theme_minimal() +
    theme_cowplot() +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "none", legend.title = element_blank()) +
    background_grid(major = 'xy')
  
}

########################################## #stack and save
# output_dir <- "/Users/maya/Documents/Figures/Code_Figures/Retrospective_supplemental_figures"

# Loop through the diseases and corresponding models to create and save stacked plots as PNGs
for (i in seq_along(diseases_truth)) {
  
  # Get the truth plot and forecast plots for the current disease/model
  truth_plot <- plot_list_pt[[diseases_truth[i]]]
  forecast_plot <- plot_list_cov_50[[models_forecast[i]]]
  improve_plot <- plot_list_cov_95[[models_forecast[i]]]
  
  # Stack the plots
  stacked_plot <- plot_grid(truth_plot, forecast_plot, improve_plot, ncol = 1, align = "v", labels = "AUTO", rel_heights = c(1, 1, 1))
  
  print(stacked_plot)
  # Define file name for the PNG output
  # file_name <- paste0(i, "_retro_coverage.png")
  # 
  # # Save the stacked plot as a PNG file
  # ggsave(filename = file_name, plot = stacked_plot, path = output_dir, width = 10, height = 7, units = "in", bg = "white")
}

###################
#PIC Plot location
###################

################################### #read in data and do some cleaning
score_all_l <- read_csv("results/location_score.csv")

score_all_l <- score_all_l |>
  mutate(cov_50 = (coverage_50*100),
         cov_95 = (coverage_95*100))

##################################### #graph

# use this one
library(scales)

# Custom colors for each model
model_colors <- c("COVID-19_forecast_model" = "#66C2A5", 
                  "Influenza_forecast_model" = "#8DA0CB", 
                  "RSV_forecast_model" = "#E78AC3")

custom_labels <- c(
  "COVID-19_forecast_model" = "COVID-19",
  "Influenza_forecast_model" = "Influenza",
  "RSV_forecast_model" = "RSV"
)

PIC_location <- score_all_l %>%
  group_by(location, model) %>%
  summarise(mean_cov_50 = mean(coverage_50, na.rm = TRUE),
            mean_cov_95 = mean(coverage_95, na.rm = TRUE), .groups = "drop")


# Compute the mean coverage for each model
mean_coverage <- score_all_l %>%
  group_by(model) %>%
  summarise(mean_cov_50 = mean(coverage_50, na.rm = TRUE),
            mean_cov_95 = mean(coverage_95, na.rm = TRUE)) 

# Get ordering based on coverage_95 for UGA_flucast-INFLAenza
uga_order <- PIC_location %>%
  filter(model == "COVID-19_forecast_model") %>%
  arrange(mean_cov_95) %>%
  pull(location)

# Apply this ordering to all models
PIC_location <- PIC_location %>%
  mutate(location = factor(location, levels = uga_order))


# Generate plot for Coverage 50%
plot_cov_50 <- ggplot(data = PIC_location, aes(x = location, y = mean_cov_50, color = model)) +
  geom_hline(yintercept = .50, linetype = "dashed", color = "azure4", size = 1) +  # Reference line
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors, labels = custom_labels) +
  labs(x = NULL, y = "50% Coverage", color = "Model") +
  theme_minimal() +
  theme_cowplot() +
  scale_y_continuous(labels = percent) +
  theme(
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
plot_cov_95 <- ggplot(data = PIC_location, aes(x = location, y = mean_cov_95, color = model)) +
  geom_hline(yintercept = .95, linetype = "dashed", color = "azure4", size = 1) +  # Reference line
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors, labels = custom_labels) +
  labs(x = NULL, y = "95% Coverage", color = "Model") +
  theme_minimal() +
  theme_cowplot() +
  scale_y_continuous(labels = percent) +
  theme(
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

#ggsave("PIC_Retro_location_lines2.png", plot = final_plot, width = 11, height = 10, dpi = 300, bg = "white") 


##################
#horizon PIC
##################

####################### #read data and clean 
score_all_h <- read_csv("results/horizon_score.csv")

score_all_h <- score_all_h |>
  mutate(cov_50 = (coverage_50*100),
         cov_95 = (coverage_95*100))

############################## # graph it out 
# Custom colors for each model
model_colors <- c("COVID-19_forecast_model" = "#66C2A5", 
                  "Influenza_forecast_model" = "#8DA0CB", 
                  "RSV_forecast_model" = "#E78AC3")

custom_labels <- c(
  "COVID-19_forecast_model" = "COVID-19",
  "Influenza_forecast_model" = "Influenza",
  "RSV_forecast_model" = "RSV"
)

# Generate individual plots
plot_cov_50 <- ggplot(data = score_all_h , aes(x = horizon, y = coverage_50, color = model)) +
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
  scale_y_continuous(labels = scales::percent) +
  
  theme(
    legend.position = "none",
    legend.title = element_blank(),   # Remove the legend title
    axis.text.x = element_text(angle = 00, hjust = 1, vjust = 0.5) 
  )

plot_cov_50 <- plot_cov_50 + background_grid(major = 'xy')

plot_cov_95 <- ggplot(data = score_all_h, aes(x = horizon, y = coverage_95, color = model)) +
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
  scale_y_continuous(labels = scales::percent) +
  
  theme(
    legend.position = "none",
    legend.title = element_blank(),   # Remove the legend title
    axis.text.x = element_text(angle = 00, hjust = 1, vjust = 0.5) 
  )

plot_cov_95 <- plot_cov_95 + background_grid(major = 'xy')

# Stack the plots
#stacked_plots <- plot_grid(plot_cov_50, plot_cov_95, ncol = 1, align = "v", labels = "AUTO")
#
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

#ggsave("PIC_Retro_horizon2.png", plot = final_plot, width = 10, height = 6, dpi = 300, bg = "white") 


################
#Horizon figure
################

##################### # read in and clean data 
simple <- score_all_h

simple_ <- simple |>
  select(model, horizon, r_wis, improvement, percent)

score_all <- simple_ %>%
  mutate(model = recode(model,
                        "COVID-19_forecast_model" = "COVID-19",
                        "Influenza_forecast_model" = "Influenza",
                        "RSV_forecast_model" = "RSV"))

############################### #graph 

#date_evi$model <- factor(date_evi$model, levels = c("UGA_flucast-INFLAenza", "UGA_flucast-Copycat", "UGA_flucast-SIRSea"))

custom_colors <- c(
  "COVID-19" = "#66C2A5",
  "Influenza" = "#8DA0CB",
  "RSV" = "#E78AC3"
)

# Create the plot
slide12 <- ggplot(score_all, aes(x = horizon, y = percent, color = model)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "azure4", size = 1) +
  labs(
    x = "Horizon",
    y = "Relative improvement (Baseline)") +
  theme_minimal() +
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  scale_color_manual(values = custom_colors) +  
  scale_x_continuous(breaks = unique(score_all$horizon)) +
  scale_y_continuous(labels = scales::percent) 

slide12 <- slide12 + background_grid(major = 'xy', minor = 'y')


print(slide12)

# Save the plot with white background using save_plot
# save_plot("sup_hor_simple.png", slide12, base_width = 8, base_asp = 1.7, bg = "white")

