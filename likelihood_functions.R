compute_reach_estimated_segments <- function(dataset) {
  responses_per_segment <- dataset %>%
    group_by(estimated_gender) %>%
    summarise(count = sum(response), .groups = 'drop') %>%
    group_by(estimated_gender) %>%
    mutate(probability = count / sum(count))
  
  female <- ifelse(dataset$estimated_gender == 'female', 1, 0)
  male <- ifelse(dataset$estimated_gender == 'male', 1, 0)
  est_seg <- bind_cols(female,male)
  
  count_female <- sum(female)
  count_male <- sum(male)
  # true_female <- count_female*p_is[1, 1] + count_male*p_is[1,2]
  # true_male <- count_male*p_is[2, 2] + count_male*p_is[2,1]
  reach_female <- responses_per_segment$count[1]/count_female
  reach_male <- responses_per_segment$count[2]/count_male
  print(paste("Reach male", reach_male,". Reach female:", reach_female))
  
}

compute_p_z_given_s <- function(dataset, segmentation){
  # Not used yet, can maybe be used as prior
  # p_s <- universe_estimates %>%
  #   group_by(gender_bucket) %>%
  #   summarise(count = sum(num_persons)) %>%
  #   group_by(gender_bucket) %>%
  #   mutate(probability = count / (universe_estimates$tot_persons[1]))
  
  # Clean missing true data
  dataset <- dataset %>%
    filter(complete.cases(.))
  
  if (segmentation == "gender") {
    conditional_probs <- dataset %>%
      group_by(estimated_gender, true_gender) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(estimated_gender) %>%
      mutate(probability = count / sum(count))
  }
  
  if (segmentation == "age") {
    conditional_probs <- dataset %>%
      group_by(estimated_age, true_age) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(estimated_age) %>%
      mutate(probability = count / sum(count))
  }
  
  dim_matrix <- sqrt(length(conditional_probs$probability))
  
  print(conditional_probs)
  
  conditional_probs_matrix <- matrix(conditional_probs$probability, dim_matrix, dim_matrix)
  
  # female <- ifelse(dataset$graph_gender == 0, 1, 0)
  # male <- ifelse(dataset$graph_gender == 1, 1, 0)
  # est_seg <- bind_cols(female,male)
  
  p_z_given_s <- t(conditional_probs_matrix)
  
  print(p_z_given_s)
  
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

compute_segment_sizes_gender <- function(dataset, segmentation) {
  head(dataset)
  
  if (segmentation == "gender") {
    segments_response <- dataset %>%
      group_by(estimated_gender) %>%
      summarise(response_count = sum(response),no_response_count = sum(1-response), .groups = 'drop')
  }
  
  if (segmentation == "age") {
    segments_response <- dataset %>%
      group_by(estimated_age) %>%
      summarise(response_count = sum(response),no_response_count = sum(1-response), .groups = 'drop')
  }
  
  print(segments_response)
  mat_segments_response <- matrix(c(segments_response$response_count, segments_response$no_response_count), 
                                  ncol = 2)

  print(mat_segments_response)
  return(mat_segments_response)
}

loglikelihood_segments_based <- function(beta, p_z_given_s, segment_responses) {
  log_value <- 0
  segment_count <- dim(p_z_given_s)[1]
  print(segment_count)
  for (i in 1:segment_count) {
    for (j in 1:2) {
      response_j <- ifelse(j == 1, 1, 0)
      value_segment <- 0
      for (s in 1:segment_count) {
        e_z_beta <- exp(beta[s])
        value_segment <- value_segment + p_z_given_s[i,s] * (((e_z_beta/(1 + e_z_beta))**response_j)*(1/(1 + e_z_beta))**(1 - response_j))
      }
      log_value <- log_value + segment_responses[i,j] * log(value_segment)
    }
  }
  return(-1*log_value)
}

optimize_loglikelihood <- function(dataset, segmentation) {
  p_z_given_s <- compute_p_z_given_s(dataset, segmentation = segmentation)
  response_segments <- compute_segment_sizes_gender(dataset, segmentation = segmentation)
  
  segment_count <- dim(p_z_given_s)[1]
  initial_par <- rep(0, segment_count)
  best_result <- optim(par=initial_par, fn=loglikelihood_segments_based, p_z_given_s = p_z_given_s, segment_responses = response_segments)
  
  if (segmentation == 'gender') {
    reach_female <- exp(best_result$par[1])/(1+exp(best_result$par[1]))
    reach_male <- exp(best_result$par[2])/(1+exp(best_result$par[2]))
    cat(paste('Reach of true segmants are ,\n for male: ', reach_male, '\n for female: ', reach_female, '\n with beta:', 
              best_result$par[1], ' ', best_result$par[2]))
  }
  if (segmentation == 'age') {
    reach_1 <- exp(best_result$par[1])/(1+exp(best_result$par[1]))
    reach_2 <- exp(best_result$par[2])/(1+exp(best_result$par[2]))
    reach_3 <- exp(best_result$par[3])/(1+exp(best_result$par[3]))
    reach_4 <- exp(best_result$par[4])/(1+exp(best_result$par[4]))
    cat(paste(reach_1, '\n', reach_2, '\n', reach_3, '\n', reach_4))
  }
  return(best_result)
}
