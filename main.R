#Clean environment
rm(list = ls())

#Load libraries and functions
prepare_setup <- function() {
  source("requirements.R")
  install_requirements()
}


prepare_setup()

n <- 1000000
dataset <- simulate_dataset_male_female(n, fraction_male = 0.484, exposure_male = 0.00769, exposure_female = 0.0131 ,
                                        estimate_correct_male = 0.603 , estimate_correct_female = 0.688)

real_dataset <- read_exposures(site_id_input = "1438679", platform_input = "MBL")
head(real_dataset)

x <- optimize_loglikelihood(dataset, segmentation = 'gender', with_prior = FALSE)

x_2 <- optimize_loglikelihood(real_dataset, segmentation = 'gender', with_prior = FALSE)
x_2 <- optimize_loglikelihood(real_dataset, segmentation = 'gender', with_prior = TRUE)

compute_reach_estimated_segments(real_dataset)

create_summary_estimates(real_dataset, c('gender'))
create_summary_true(real_dataset, c('age'))

create_bubble_chart(real_dataset)
insight_exposed_individuals(real_dataset, c('demo'))

evaluation(real_dataset, x_2$par, 'gender')
