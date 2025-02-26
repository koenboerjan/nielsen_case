# Install required packages
# -----------------------


install_requirements <- function() {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(plotly)
  library(data.table)
  library(future.apply)
  library(duckdb)
  library(DBI)
  source("likelihood_functions.R")
  source("simulate_data.R")
  source("read_data.R")
  source("obtain_insights_data.R")
  source("simulation_study.R")
  source("likelihood_functions_frequency.R")
  source("p_z_given_s_matrices.R")
}