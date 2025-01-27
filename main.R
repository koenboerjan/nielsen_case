#Clean environment
rm(list = ls())

#Load libraries and functions
prepare_setup <- function() {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  source("likelihood_functions.R")
  source("simulate_data.R")
  source("read_data.R")
}


prepare_setup()

n <- 1000000
dataset <- simulate_dataset_male_female(n, 0.1, 0.1, 0.6, 0.7)

real_dataset <- read_exposures()
head(real_dataset)

x <- optimize_loglikelihood(dataset, segmentation = 'gender')

x2 <- optimize_loglikelihood(real_dataset, segmentation = 'age')
compute_reach_estimated_segments(real_dataset)
