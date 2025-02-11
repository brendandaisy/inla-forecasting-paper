library(covidHubUtils)
library(scoringutils)
library(gt)
library(tidyverse)
library(cowplot)
library (ggplot2)

##############
#date
##############

#truth data
truth_df_all <- read_csv("../data/truth_df.csv")

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

library(cowplot)

## plot top panel 

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
      y = "US Counts") +
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

# You can also save each plot individually
#for (mod in names(plot_list)) {
#save_plot(paste0("../Figures/sup_date_", mod, ".png"), 
#            plot_list[[mod]], base_width = 8, base_asp = 1.7, bg = "white")
#}

## date info
score_all_d <- read_csv("../data/simple_scored_date.csv")

score_all_d <- score_all_d |>
  filter(date <= "2024-04-27")

score_all_d <- score_all_d %>%
  group_by(model) %>%
  complete(date = seq(min(date), max(date), by = "1 week")) %>%
  fill(model, .direction = "down") %>% # Ensure model column is filled in each group
  ungroup()

#plot

library(ggplot2)
library(cowplot)
library(dplyr)

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
    labs(x = NULL, y = "Coverage 50%") +
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
    labs(x = NULL, y = "Coverage 95%") +
    scale_x_date(date_breaks = "3 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
    theme_minimal() +
    theme_cowplot() +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.position = "none", legend.title = element_blank()) +
    background_grid(major = 'xy')
  
}


#######
#plot date
#######
# Define the directory to save PNG files
#output_dir <- "../Figures/"

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
  #file_name <- paste0(i, "_retro_coverage.png")
  
  # Save the stacked plot as a PNG file
  #ggsave(filename = file_name, plot = stacked_plot, path = output_dir, width = 10, height = 6, units = "in", bg = "white")
}

# Close the PDF device
#dev.off()



########
#Location
##########
score_all_l <- read_csv("../data/simple_scored_location.csv")

score_all_l <- score_all_l |>
  mutate(cov_50 = (coverage_50*100),
         cov_95 = (coverage_95*100))


model_colors <- c("COVID-19_forecast_model" = "#66C2A5", 
                  "Influenza_forecast_model" = "#8DA0CB", 
                  "RSV_forecast_model" = "#E78AC3")

#####
#plot location
#######

# Generate individual plots
plot_cov_50 <- ggplot(data = score_all_l, aes(x = location, y = coverage_50, color = model)) +
  # Add a horizontal line at y = 0
  geom_hline(yintercept = .50, linetype = "dashed", color = "azure4", size = 1) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors) +
  labs(x = NULL, y = "Coverage 50%") +
  #scale_x_date(date_breaks = "1 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
  #labs(y = "cov_50") +
  theme_minimal() +
  theme_cowplot() +
  scale_y_continuous(labels = scales::percent) +
  
  theme(
    legend.position = "none",
    legend.title = element_blank(),   # Remove the legend title
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  )

plot_cov_50 <- plot_cov_50 + background_grid(major = 'xy')

plot_cov_95 <- ggplot(data = score_all_l, aes(x = location, y = coverage_95, color = model)) +
  # Add a horizontal line at y = 0
  geom_hline(yintercept = .95, linetype = "dashed", color = "azure4", size = 1) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors) +
  labs(x = NULL, y = "Coverage 95%") +
  #scale_x_date(date_breaks = "1 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
  #labs(y = "cov_95") +
  theme_minimal()+
  theme_cowplot() +
  scale_y_continuous(labels = scales::percent) +
  
  theme(
    legend.position = "none",
    legend.title = element_blank(),   # Remove the legend title
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  )

plot_cov_95 <- plot_cov_95 + background_grid(major = 'xy')

# Stack the plots
stacked_plots <- plot_grid(plot_cov_50, plot_cov_95, ncol = 1, align = "v", labels = "AUTO")

# Display the stacked plots
print(stacked_plots)

#ggsave("../Figures/PIC_Retro_location.png", plot = stacked_plots, width = 8, height = 10, dpi = 300, bg = "white") 

#############
#Horizon
############
score_all_h <- read_csv("../data/simple_scored_horizon.csv")

score_all_h <- score_all_h |>
  mutate(cov_50 = (coverage_50*100),
         cov_95 = (coverage_95*100))

#######
#Horizon plot
######
model_colors <- c("COVID-19_forecast_model" = "#66C2A5", 
                  "Influenza_forecast_model" = "#8DA0CB", 
                  "RSV_forecast_model" = "#E78AC3")

# Generate individual plots
plot_cov_50 <- ggplot(data = score_all_h , aes(x = horizon, y = coverage_50, color = model)) +
  # Add a horizontal line at y = 0
  geom_hline(yintercept = .50, linetype = "dashed", color = "azure4", size = 1) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors) +
  labs(x = NULL, y = "Coverage 50%") +
  #scale_x_date(date_breaks = "1 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
  #labs(y = "cov_50") +
  theme_minimal() +
  theme_cowplot() +
  scale_y_continuous(labels = scales::percent) +
  
  theme(
    legend.position = "none",
    legend.title = element_blank(),   # Remove the legend title
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  )

plot_cov_50 <- plot_cov_50 + background_grid(major = 'xy')

plot_cov_95 <- ggplot(data = score_all_h, aes(x = horizon, y = coverage_95, color = model)) +
  # Add a horizontal line at y = 0
  geom_hline(yintercept = .95, linetype = "dashed", color = "azure4", size = 1) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = model_colors) +
  labs(x = NULL, y = "Coverage 95%") +
  #scale_x_date(date_breaks = "1 months", date_labels = "%b '%y", guide = guide_axis(angle = 0)) +
  #labs(y = "cov_95") +
  theme_minimal()+
  theme_cowplot() +
  scale_y_continuous(labels = scales::percent) +
  
  theme(
    legend.position = "none",
    legend.title = element_blank(),   # Remove the legend title
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  )

plot_cov_95 <- plot_cov_95 + background_grid(major = 'xy')

# Stack the plots
stacked_plots <- plot_grid(plot_cov_50, plot_cov_95, ncol = 1, align = "v", labels = "AUTO")

# Display the stacked plots
print(stacked_plots)

#ggsave("../Figures/PIC_Retro_horizon.png", plot = stacked_plots, width = 10, height = 6, dpi = 300, bg = "white") 