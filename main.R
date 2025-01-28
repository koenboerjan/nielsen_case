#Clean environment
rm(list = ls())

#Load libraries and functions
prepare_setup <- function() {
  source("requirements.R")
  install_requirements()
}


prepare_setup()

n <- 1000000
dataset <- simulate_dataset_male_female(n, fraction_male = 0.7, exposure_male = 0.03, exposure_female = 0.05,
                                        estimate_correct_male = 0.6, estimate_correct_female = 0.7)

real_dataset <- read_exposures()
head(real_dataset)

x <- optimize_loglikelihood(dataset, segmentation = 'gender')

x2 <- optimize_loglikelihood(real_dataset, segmentation = 'age')
compute_reach_estimated_segments(real_dataset)
