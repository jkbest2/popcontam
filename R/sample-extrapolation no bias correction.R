# Load necessary libraries
library(dplyr)
library(tidyr)
library(writexl)

# Read the data
data <- read.csv("data/PCB_Puyallup_monitoring_data.csv")

# Log-transform the 'pcb_uggwet' column using natural log
data <- data %>% mutate(log_pcb_uggwet = log(pcb_uggwet))  

# Calculate the overall variance of log-transformed concentrations
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
    corrected_variance = ifelse(Composite_Count > 1, overall_variance, 0)  # Using overall variance (no correction), no variance for individual samples
  ) %>%
  ungroup()


# Generate individual concentrations
set.seed(123)  # For reproducibility

data <- data %>%
  rowwise() %>%
  mutate(
    individual_concentrations = list(simulate_individual_concentrations(pcb_uggwet, corrected_variance, Composite_Count))
  ) %>%
  ungroup()

expanded_data <- data %>%
  unnest(individual_concentrations)

# Unnest the 'individual_concentrations' column to create a new dataframe with only the concentrations
individual_concentrations_df <- data %>%
  unnest(individual_concentrations) %>%
  select(individual_concentration = individual_concentrations)

# View the resulting dataframe
head(individual_concentrations_df)

