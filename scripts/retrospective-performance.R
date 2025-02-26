# --------------------------------------------------------------------------------
# display performance of the retrospective analysis (Figure 3)--------------------
# --------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(cowplot)
library(covidHubUtils)
library(scoringutils)
library(usmap)

#############
#Left hand panel
#############

##############################  #read in data 
truth_df <- read.csv("/Users/maya/Desktop/figures_paper/data/truth_df_fig3.csv")
plot_df_all <- read.csv("/Users/maya/Desktop/figures_paper/data/plot_df_fig3.csv")

plot_df_all <- plot_df_all |>
  mutate(date = as.Date(date))

######################################  # Function for plotting 
generate_plots <- function(plot_df, truth_df) {
  
  # Process the forecast data to keep only specified quantiles
  forecast_df <- plot_df %>%
    spread(quantile, prediction) 
  
  truth_df <- truth_df %>%
    filter(date >= ymd("2022-08-30"))
  
  # Create a list to store the plots for each location
  plots_list <- list()
  
  for(curr_location in unique(forecast_df$location)) {
    
    location_forecast <- forecast_df %>%
      filter(location == curr_location)
    
    location_truth <- truth_df %>%
      filter(location == curr_location)
    
    # Suppress warnings about NAs in the plot - NAs are supposed to be in there 
    suppressWarnings({
      
      # Create individual plots for each disease with specific y-axis limits
      covid_plot <- ggplot() +
        geom_point(data = location_truth %>% filter(disease == "covid"), aes(x = date, y = count), color = "black", size = 0.5) +
        geom_ribbon(data = location_forecast %>% filter(disease == "covid"), aes(x = date, ymin = `0.025`, ymax = `0.975`), fill = "#7DC0A6", alpha = 0.3) +
        geom_line(data = location_forecast %>% filter(disease == "covid"), aes(x = date, y = `0.5`), color = "#7DC0A6") +
        scale_y_continuous(limits = c(0, 65000), oob = scales::rescale_none, labels = scales::comma) +  # Custom limits for COVID
        theme_minimal() +
        theme_cowplot() +
        labs(x = NULL, y = "Hospital admits") +
        scale_x_date(date_breaks = "3 months", date_labels = "%b '%y", guide = guide_axis(angle = 0), expand = c(0, 0)) +
        background_grid(major = 'xy', minor = 'y') +
        theme(legend.position = "none",
              axis.title.x = element_blank(), # Remove x-axis labels from here
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 10),
              #plot.title = element_text(size = 10, face = "plain", hjust = 0.5),
              plot.title = element_blank(),  # Remove the title
              strip.background = element_blank(),
              strip.text = element_text(face = "bold", size = 10)) +
        ggtitle("COVID-19")
      
      flu_plot <- ggplot() +
        geom_point(data = location_truth %>% filter(disease == "flu"), aes(x = date, y = count), color = "black", size = 0.5) +
        geom_ribbon(data = location_forecast %>% filter(disease == "flu"), aes(x = date, ymin = `0.025`, ymax = `0.975`), fill = "#919FC7", alpha = 0.3) +
        geom_line(data = location_forecast %>% filter(disease == "flu"), aes(x = date, y = `0.5`), color = "#919FC7") +
        scale_y_continuous(limits = c(0, 35000), oob = scales::rescale_none, labels = scales::comma) +  # Custom limits for Flu
        theme_minimal() +
        theme_cowplot() +
        labs(x = NULL, y = "Hospital admits") +
        scale_x_date(date_breaks = "3 months", date_labels = "%b '%y", guide = guide_axis(angle = 0), expand = c(0, 0)) +
        background_grid(major = 'xy', minor = 'y') +
        theme(legend.position = "none",
              axis.title.x = element_blank(), # Remove x-axis labels from here
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 10),
              #plot.title = element_text(size = 10, face = "plain", hjust = 0.5),
              plot.title = element_blank(),  # Remove the title
              strip.background = element_blank(),
              strip.text = element_text(face = "bold", size = 10)) +
        ggtitle("Influenza")
      
      rsv_plot <- ggplot() +
        geom_point(data = location_truth %>% filter(disease == "rsv"), aes(x = date, y = count), color = "black", size = 0.5) +
        geom_ribbon(data = location_forecast %>% filter(disease == "rsv"), aes(x = date, ymin = `0.025`, ymax = `0.975`), fill = "#DA8EC0", alpha = 0.3) +
        geom_line(data = location_forecast %>% filter(disease == "rsv"), aes(x = date, y = `0.5`), color = "#DA8EC0") +
        scale_y_continuous(limits = c(0, 20000), oob = scales::rescale_none, labels = scales::comma) +  # Custom limits for COVID
        theme_minimal() +
        theme_cowplot() +
        labs(x = NULL, y = "Hospital admits") +
        scale_x_date(date_breaks = "3 months", date_labels = "%b '%y", guide = guide_axis(angle = 0), expand = c(0, 0)) +
        background_grid(major = 'xy', minor = 'y') +
        theme(legend.position = "none",
              #axis.title.x = element_blank(), # Remove x-axis labels from here
              #axis.text.x = element_blank(),
              #axis.ticks.x = element_blank(),
              strip.background = element_blank(),
              axis.text.x = element_text(size = 10),  # Change x-axis text size
              axis.text.y = element_text(size = 10),  # Change y-axis text size
              #plot.title = element_text(size = 10, face = "plain", hjust = 0.5),
              plot.title = element_blank(),  # Remove the title
              strip.text = element_text(face = "bold", size = 10)) +
        ggtitle("RSV")
      
      # Combine the plots with shared x-axis using plot_grid
      combined_plot <- plot_grid(covid_plot, flu_plot, rsv_plot, ncol = 1, align = "v", rel_heights = c(1, 1, 1.2))
      # Store both plots in the list with the location name as the key
      #plots_list[[curr_location]] <- list(covid_plot = covid_plot, flu_plot = flu_plot)
      plots_list[[curr_location]] <- combined_plot
    })
  }
  
  # Return the list of plots for each location
  return(plots_list)
  
}

###################################### #Create gaps for Flu and RSV (exclude summer months)
# Define the date range and diseases of interest
start_date <- as.Date("2023-04-22")
end_date <- as.Date("2023-09-09")
target_diseases <- c("flu", "rsv")

# Modify the dataframe
plot_df_all_gaps <- plot_df_all %>%
  mutate(
    prediction = if_else(
      date >= start_date & date <= end_date & disease %in% target_diseases,
      NA_real_,  # Set to NA for numeric columns
      prediction
    )
  )

plot_df_all_gaps <- plot_df_all_gaps %>% mutate(date = as.Date(date))
truth_df <- truth_df %>% mutate(date = as.Date(date))
##############################  #generate and save plots (PLEASE CHANGE path)

plots_all <- generate_plots(plot_df_all_gaps, truth_df) #changed to flu so green

combined_plot <- plot_grid(plotlist = plots_all, ncol = 1, align = "v", rel_heights = c(1, 1, 5))
print(combined_plot)


################
#Right panel (maps)
#################

##############################  #needed functions
make_rwis <- function(df, baseline_model) {
  df <- df %>%
    #group_by(location, date) %>%
    group_by(location) %>%
    mutate(
      baseline_wis = interval_score[model == baseline_model],
      r_wis = ifelse(model == baseline_model, 1, interval_score / baseline_wis)
    ) %>%
    ungroup() %>%
    select(-baseline_wis) # remove temp column
  
  return(df)
}

##############################  #Read in data
score_graph_location_covid <- read_csv("/Users/maya/Desktop/figures_paper/data/score_graph_location_covid_fig3.csv")
score_graph_location_flu <- read_csv("/Users/maya/Desktop/figures_paper/data/score_graph_location_flu_fig3.csv")
score_graph_location_rsv <- read_csv("/Users/maya/Desktop/figures_paper/data/score_graph_location_rsv_fig3.csv")

############################################ #
# Prepare the US map data with FIPS codes
us <- us_map(regions = "states") %>% 
  st_as_sf()  # Convert to sf object

iso_states <- c("Alaska", "Hawaii", "Puerto Rico", "District of Columbia", "Rhode Island")

# Define custom positions for small or isolated states
us_pts <- tribble(
  ~abbr, ~full, ~x, ~y,
  "AK", "Alaska", 2350000, -500000,
  "DC", "District of Columbia", 2350000, -700000,
  "RI", "Rhode Island", 2700000, -700000,
  "HI", "Hawaii", 2350000, -900000,
  "PR", "Puerto Rico", 2700000, -900000
) |> 
  st_as_sf(crs = st_crs(us), coords = c("x", "y"))

# Define the plot function
plot_loc_effect <- function(map_data, us, us_pts = NULL) {
  ret <- map_data %>%
    ggplot() +
    geom_sf(aes(fill = r_wis), col = "black", size = 0.1) +
    scale_fill_gradient2(
      name = "rWIS",
      #low = "#BE94E6",
      low = "#253494" ,
      #low = "#AB76DE", 
      mid = "white",
      high = "salmon2",
      #high = "red4" ,
      midpoint = 1,  # Set the midpoint for scaling
      limits = c(0.4, 1.65),  # Define the limits
      breaks = c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6),  # Breaks for the legend
      labels = c("0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6"),  # Labels for the legend
      na.value = "grey70"  # Set NA values to grey
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      #legend.position = "none"
      legend.text = element_text(size = 8),      # Adjust legend text size
      legend.title = element_text(size = 12)     # Adjust legend title size
    )
  
  if (!is.null(us_pts)) {
    # Verify if `r_wis` or equivalent column exists in `us_pts`
    if ("r_wis" %in% names(us_pts)) {
      ret <- ret +
        geom_sf(data = us_pts, aes(fill = r_wis), shape = 21, size = 5, col = "white") +
        geom_sf_text(data = us_pts, aes(label = abbr), nudge_x = 150000, size = 4)
    } else {
      stop("Error: `r_wis` column not found in `us_pts`.")
    }
  }
  return(ret)
}

############################################ #Covid plot
map_df <- score_graph_location_covid %>%
  filter(model == "covid") %>%
  mutate(fips = fips(location)) # Convert state names to FIPS codes

map_data <- us_map(regions = "states") %>%
  left_join(map_df, by = "fips") |>
  filter(!full %in% c("Alaska", "Hawaii"))

r_wis_iso <- map_data |>
  filter(full %in% iso_states) |>
  select(abbr, full, r_wis)

df_covid <- data.frame(
  #abbr = c("AK", "DC", "HI", "RI", "PR"),
  full = c("Alaska", "District of Columbia", "Hawaii", "Rhode Island", "Puerto Rico"),
  r_wis = c(0.8394664, 0.5309340, 1.0024976, 0.5516186, 0.5730461))

us_pts_covid <- left_join(us_pts, df_covid, by = "full") 

# Generate the map
p_covid <- plot_loc_effect(map_data, us, us_pts_covid)
print(p_covid)

# Save the plot as a PNG
ggsave(
  filename = "covid_map_p.png",   # File name
  plot = p_covid,              # Plot object
  width = 10,                # Width in inches
  height = 8,                # Height in inches
  dpi = 600,                  # Resolution in dots per inch
  bg = "white"
)

############################################ #Flu plot
map_df <- score_graph_location_flu %>%
  filter(model == "flu") %>%
  mutate(fips = fips(location)) # Convert state names to FIPS codes

map_data <- us_map(regions = "states") %>%
  left_join(map_df, by = "fips") |>
  filter(!full %in% c("Alaska", "Hawaii"))

r_wis_iso <- map_data |>
  filter(full %in% iso_states) |>
  select(abbr, full, r_wis)

df_flu <- data.frame(
  #abbr = c("AK", "DC", "HI", "RI", "PR"),
  full = c("Alaska", "District of Columbia", "Hawaii", "Rhode Island", "Puerto Rico"),
  r_wis = c(0.8188897, 0.6378093, 0.7721396, 0.5542517, 1.6152894))

us_pts_flu <- left_join(us_pts, df_flu, by = "full") 

# Generate the map
p_flu <- plot_loc_effect(map_data, us, us_pts_flu)
print(p_flu)

# Save the plot as a PNG
ggsave(
  filename = "flu_map_p.png",   # File name
  plot = p_flu,              # Plot object
  width = 10,                # Width in inches
  height = 8,                # Height in inches
  dpi = 300                  # Resolution in dots per inch
)

############################################ #RSV plot
map_df <- score_graph_location_rsv %>%
  filter(model == "rsv_2018") %>%
  mutate(fips = fips(location)) # Convert state names to FIPS codes

map_data <- us_map(regions = "states") %>%
  left_join(map_df, by = "fips") |>
  filter(!full %in% c("Alaska", "Hawaii"))

map_data <- map_data %>%
  mutate(r_wis = as.numeric(r_wis))

df_RSV <- data.frame(
  #abbr = c("AK", "DC", "HI", "RI", "PR"),
  full = c("Alaska", "District of Columbia", "Hawaii", "Rhode Island", "Puerto Rico"),
  r_wis = c(NA, NA, NA, NA, NA))

us_pts_RSV <- left_join(us_pts, df_RSV, by = "full") 

us_pts_RSV <- us_pts_RSV %>%
  mutate(r_wis = as.numeric(r_wis))

# Generate the map
p_RSV <- plot_loc_effect(map_data, us, us_pts_RSV)
print(p_RSV)

# Save the plot as a PNG
ggsave(
  filename = "RSV_map_P.png",   # File name
  plot = p_RSV,              # Plot object
  width = 10,                # Width in inches
  height = 8,                # Height in inches
  dpi = 300                  # Resolution in dots per inch
)

############################################ #combine and print
object <- plot_grid(combined_plot, plot_grid(p_covid, p_flu, p_RSV, ncol = 1))

save_plot(
  "fig2_test_new2.png", 
  object, 
  base_width = 10,  # Try increasing width
  base_height = 6,  # Adjust for the height aspect
  base_asp = 1.8,  # Adjust aspect ratio to better fit content
  bg = "white"
)

