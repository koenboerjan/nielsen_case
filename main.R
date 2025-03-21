#Clean environment
rm(list = ls())

#Load libraries and functions
prepare_setup <- function() {
  source("requirements.R")
  install_requirements()
}


prepare_setup()

n <- 1000000
dataset <- simulate_dataset(n, c(0.516, 0.484), c(0.02,0.06), matrix(c(0.7,0.3,0.2,0.8), ncol = 2),
                                       "gender", seed_number = 0, print_simulation = TRUE) 

dataset_test <- read_only_exposures()
dataset_demos <- combine_panel_id_graph() 

compute_p_z_given_s_including_prior(dataset_test, 'gender')
compute_p_z_given_s(dataset_test, 'gender')

real_dataset <- read_exposures()
real_dataset <- read_exposures_frequency()
head(real_dataset)

x <- optimize_loglikelihood(dataset, segmentation = 'gender', with_prior = FALSE)

x_2 <- optimize_loglikelihood(real_dataset, segmentation = 'demo', with_prior = FALSE)
x_2 <- optimize_loglikelihood(real_dataset, segmentation = 'gender', with_prior = TRUE)

compute_reach_estimated_segments(real_dataset)

create_summary_estimates(dataset_demos, c('gender'))
create_summary_true(dataset_demos, c('gender'))

create_bubble_chart(dataset_demos)
insight_exposed_individuals(real_dataset, c('age'))

simulate_evaluate_fraction_2_segments()
simulate_evaluate_fraction_3_segments()
simulate_evaluate_exposure_2_segments()
simulate_evaluate_accuracy_2_segments()

evaluation(real_dataset, x_2$par, 'gender')
