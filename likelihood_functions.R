compute_reach_estimated_segments <- function(dataset) {
  responses_per_segment <- dataset %>%
    group_by(graph_gender) %>%
    summarise(count = sum(response), .groups = 'drop') %>%
    group_by(graph_gender) %>%
    mutate(probability = count / sum(count))
  
  female <- ifelse(test_df$Gender_est == 0, 1, 0)
  male <- ifelse(test_df$Gender_est == 1, 1, 0)
  est_seg <- bind_cols(female,male)
  
  count_female <- sum(female)
  count_male <- sum(male)
  true_female <- count_female*p_is[1, 1] + count_male*p_is[1,2]
  true_male <- count_male*p_is[2, 2] + count_male*p_is[2,1]
  reach_female <- test$count[1]/count_female
  reach_male <- test$count[2]/count_male
  print(paste("Reach male", reach_male,". Reach female:", reach_female))
  
}

compute_p_z_given_s <- function(dataset){
  # Not used yet, can maybe be used as prior
  # p_s <- universe_estimates %>%
  #   group_by(gender_bucket) %>%
  #   summarise(count = sum(num_persons)) %>%
  #   group_by(gender_bucket) %>%
  #   mutate(probability = count / (universe_estimates$tot_persons[1]))
  
  
  conditional_probs <- dataset %>%
    group_by(estimated_gender, true_gender) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(estimated_gender) %>%
    mutate(probability = count / sum(count))
  conditional_probs_matrix <- matrix(conditional_probs$probability,2,2)
  
  # female <- ifelse(dataset$graph_gender == 0, 1, 0)
  # male <- ifelse(dataset$graph_gender == 1, 1, 0)
  # est_seg <- bind_cols(female,male)
  
  p_z_given_s <- t(conditional_probs_matrix)
  
  return(p_z_given_s)
}

loglikelihood <- function(beta) {
  print(beta)
  log_value <- 0
  for (i in 1:n_test) {
    # if (male[i] == 1) {
    #   test_exposure_data[i] <- 0
    # }
    log_value_between <- 0
    for (s in 1:2) {
      e_z_beta <- exp(beta[s])
      log_value_between <- log_value_between + expect_indicator[i, s] * 
        (((e_z_beta/(1 + e_z_beta))**test_exposure_data[i])*(1/(1 + e_z_beta))**(1 - test_exposure_data[i]))
    }
    log_value <- log_value + log(log_value_between)
  }
  print(log_value)
  return(-1*log_value)
}

compute_segment_sizes <- function(dataset) {
  segments_response <- dataset %>%
    group_by(estimated_gender) %>%
    summarise(response_count = sum(response),no_response_count = sum(1-response), .groups = 'drop')
  
  mat_segments_response <- matrix(c(segments_response$response_count, segments_response$no_response_count),2,2)
  
  return(mat_segments_response)
}

loglikelihood_segments_based <- function(beta, p_z_given_s, segment_responses) {
  print(beta)
  log_value <- 0
  for (i in 1:2) {
    for (j in 1:2) {
      response_j <- ifelse(j == 1, 1, 0)
      value_segment <- 0
      for (s in 1:2) {
        e_z_beta <- exp(beta[s])
        value_segment <- value_segment + p_z_given_s[s,i] * (((e_z_beta/(1 + e_z_beta))**response_j)*(1/(1 + e_z_beta))**(1 - response_j))
      }
      log_value <- log_value + segment_responses[i,j] * log(value_segment)
    }
  }
  print(log_value)
  return(-1*log_value)
}

optimize_loglikelihood <- function(dataset) {
  p_z_given_s <- compute_p_z_given_s(dataset)
  print(p_z_given_s)
  response_segments <- compute_segment_sizes(dataset)
  best_result <- optim(par=c(0,0), fn=loglikelihood_segments_based, p_z_given_s = p_z_given_s, segment_responses = response_segments)
  
  reach_female <- exp(best_result$par[1])/(1+exp(best_result$par[1]))
  reach_male <- exp(best_result$par[2])/(1+exp(best_result$par[2]))
  cat(paste('Reach of true segmants are ,\n for male: ', reach_male, '\n for female: ', reach_female, '\n with beta:', 
            best_result$par[1], ' ', best_result$par[2]))
  return(best_result)
}
