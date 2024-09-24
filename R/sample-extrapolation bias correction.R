# Load necessary libraries
library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("data/PCB_Puyallup_monitoring_data.csv")

# Log-transform the 'pcb_uggwet' column using natural log
data <- data %>% mutate(log_pcb_uggwet = log(pcb_uggwet))  

# Calculate the overall geometric mean and variance of log-transformed concentrations
overall_geom_mean <- exp(mean(data$log_pcb_uggwet, na.rm = TRUE))  # Geometric mean
overall_variance <- var(data$log_pcb_uggwet, na.rm = TRUE)  # Variance of log-transformed values

# Define the function to simulate individual concentrations
simulate_individual_concentrations <- function(geom_mean, variance, n) {
  # Simulate from log-normal distribution
  rlnorm(n, meanlog = log(geom_mean), sdlog = sqrt(variance))
}

# Apply the correct method based on Composite_Count to account for bias when pooling samples
data <- data %>%
  rowwise() %>%
  mutate(
    corrected_variance = ifelse(Composite_Count > 1, log(1 + (overall_variance / Composite_Count)), 0),  # Apply correction if Composite_Count > 1;  variance shrinks when more samples are composited
    geom_mean_pooled = ifelse(Composite_Count > 1, exp(log_pcb_uggwet - (corrected_variance / 2)), pcb_uggwet)  # Caudill's correction  is applied as pooled result may overestimate the average concentration across individuals. Correcting for this ensures that the pooled geometric mean more accurately represents the individual-level mean.
    #adjusting the pooled geometric mean to more closely represent what the individual sample concentrations would have been without the bias introduced by compositing.
    #note: Should we be calculating arithmetic mean rather than geometric mean? (current equation is calculating geometric mean)
  ) %>%
  ungroup()

# Generate individual concentrations
set.seed(123)  # For reproducibility

data <- data %>%
  rowwise() %>%
  mutate(
    individual_concentrations = list(simulate_individual_concentrations(geom_mean_pooled, corrected_variance, Composite_Count))
  ) %>%
  ungroup()

# Unnest the 'individual_concentrations' column to create a new dataframe with only the concentrations
individual_concentrations_df <- data %>%
  unnest(individual_concentrations) %>%
  select(individual_concentration = individual_concentrations)

# View the resulting dataframe
head(individual_concentrations_df)
