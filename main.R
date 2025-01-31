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

x <- optimize_loglikelihood(dataset, segmentation = 'gender',with_prior = TRUE)
x2 <- optimize_loglikelihood(real_dataset, segmentation = 'age')
compute_reach_estimated_segments(real_dataset)

conditional_probs<-compute_p_z_given_s_including_prior(dataset = real_dataset,segmentation = 'gender')
segment_responses<-compute_segment_sizes(dataset = real_dataset,segmentation = 'gender')
true_segments<-compute_true_segment_sizes(dataset = real_dataset,segmentation = 'gender')

segment_responses[1,1]

initial_betas <- c(0, 0)
maximize_likelihood <- optim(
  par = initial_betas,
  fn = loglikelihood_segments_based,
  p_z_given_s = conditional_probs,
  segment_responses = segment_responses,
  true_segments = true_segments,
  method = 'L-BFGS-B',
)

beta <- maximize_likelihood$par
print(beta)

# Compute probabilities using plogis
pr_1_given_female <- plogis(beta[1])
pr_1_given_male <- plogis(beta[2])

print(c(pr_1_given_female, pr_1_given_male))

true_data_fractions <- real_dataset %>%
  filter(complete.cases(.)) %>%  # Corrected complete.cases usage
  group_by(true_gender, response) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(true_gender) %>%
  mutate(fraction = count / sum(count))

print(true_data_fractions)
