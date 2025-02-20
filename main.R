#Clean environment
rm(list = ls())

#Load libraries and functions
prepare_setup <- function() {
  source("requirements.R")
  install_requirements()
}


prepare_setup()

real_dataset<-read_exposures()
print(real_dataset)
n <- 1000000
dataset <- simulate_dataset(n, c(0.516, 0.484), c(0.02,0.06), matrix(c(0.7,0.3,0.2,0.8), ncol = 2),
                                       "gender", seed_number = 0, print_simulation = TRUE) 

real_dataset <- read_exposures(site_id_input = "1438679", platform_input = "MBL")
head(real_dataset)
 
x <- optimize_loglikelihood(dataset, segmentation = 'gender', with_prior = FALSE)

x_2 <- optimize_loglikelihood(real_dataset, segmentation = 'gender', with_prior = FALSE)
x_2 <- optimize_loglikelihood(real_dataset, segmentation = 'gender', with_prior = TRUE)
x_3 <- optimize_loglikelihood(real_dataset, segmentation = 'age', with_prior = TRUE)
x_4 <- optimize_loglikelihood(real_dataset, segmentation = 'demo', with_prior = TRUE)
compute_reach_estimated_segments(real_dataset)

create_summary_estimates(real_dataset, c('gender'))
create_summary_true(real_dataset, c('age'))

create_bubble_chart(real_dataset)
insight_exposed_individuals(real_dataset, c('demo'))

simulate_evaluate_fraction_3_segments()
simulate_evaluate_exposure_2_segments()
simulate_evaluate_accuracy_2_segments()

evaluation(real_dataset, x_2$par, 'gender')
