library(covidHubUtils)
library(scoringutils)
library(gt)
library(tidyverse)


####################
#date 
##################


###############
#read in data 
###############
score_all_d <- read_csv("../data/simple_scored_date.csv")

score_all_d <- score_all_d |>
  mutate(average_percent = mean(percent, na.rm = TRUE))

truth_df_all <- read_csv("../data/truth_df.csv")

truth_df <- truth_df_all |>
  filter(location == "US") |>
  filter(date >= ymd("2022-08-30"),
         date <= ymd("2024-04-27"))

#TODO: FINISH FIXING 
truth_df <- truth_df %>%
  mutate(disease = recode(disease,
                          "covid" = "COVID-19 Historical Data",
                          "flu" = "Influenza Historical Data",
                          "rsv" = "RSV Historical Data"))

truth_df2 <- truth_df |>
  filter(date >= ymd("2022-09-17"))


###########
#date top panels 
###########
#use this 
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
      y = "Hospitlization Counts") +
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

# You can also save each plot individually
#for (mod in names(plot_list)) {
#save_plot(paste0("../Figures/sup_date_", mod, ".png"), 
#            plot_list[[mod]], base_width = 8, base_asp = 1.7, bg = "white")
#}


###########
#Bottom Panel date 
##########
library(cowplot)

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

#models <- unique(score_all_j$model)
models <- unique(score_all_d$model)

for (mod in models) {
  # Filter the data for each model
  #score_model <- subset(score_all_j, model == mod)
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


#############
#date graph together 
############
library(cowplot)
library(ggplot2)

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
  
  # Save the stacked plot as a PNG file
  #file_name <- paste0(i, "improvement_retro_date.png")
  #ggsave(filename = file_name, plot = stacked_plot, path = "../Figures/", width = 8.8, height = 5.5, units = "in", bg = "white")
  
  library(grid)
  library(cowplot)
  
  # Loop through the diseases and corresponding models
  for (i in seq_along(diseases_truth)) {
    
    # Get the truth plot and forecast plot for the current disease/model
    truth_plot <- plot_list_pt[[diseases_truth[i]]]
    forecast_plot <- plot_list_pb[[models_forecast[i]]]
    
    # Align the plots manually using axis_canvas() and plot_grid
    aligned_plot <- plot_grid(
      truth_plot, 
      forecast_plot, 
      ncol = 1, 
      align = "v", 
      axis = "tb"  # Align both top and bottom axes
    )
    
    # Print the aligned plot
    print(aligned_plot)
  }
  
  
}