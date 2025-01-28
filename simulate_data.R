simulate_dataset_male_female <- function(n_test, fraction_male, exposure_male, exposure_female, estimate_correct_male, 
                                         estimate_correct_female,  seed_number = 0) {
  set.seed(seed_number)
  male_response <- rbinom(n_test*fraction_male, 1, exposure_male)
  female_response <- rbinom(n_test*(1-fraction_male), 1, exposure_female)
  responses <- c(male_response, female_response)
  
  male_estimates <- rbinom(n_test*fraction_male, 1, estimate_correct_male)
  female_estimates <- rbinom(n_test*(1-fraction_male), 1, (1 - estimate_correct_female))
  estimates <- c(male_estimates, female_estimates)
  
  true_gender <- c(rep(1,n_test*fraction_male), rep(0,n_test*(1-fraction_male)))
  
  sample_df <- data.frame(estimates, true_gender, responses)
  colnames(sample_df) <- c("estimated_gender",  "true_gender", "response")
  
  true_male_reach <- sum(male_response)/(n_test*fraction_male)
  true_female_reach <- sum(female_response)/(n_test*(1-fraction_male))
  
  cat(paste('Dataset for male and female devision has been generated, input parameters are: \n sample size: ', n_test, 
              '\n exposure male:', exposure_male, '\n exposure female: ', exposure_female, 
              '\n correct estimation male: ', estimate_correct_male, '\n correct estimation female: ', 
              estimate_correct_female, '\nTrue sampled reach are, \n for male: ', true_male_reach, '\n for female: ',
              true_female_reach))
  
  return(sample_df)
}

