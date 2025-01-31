# ----------------------------- Compute Reach for Estimated Segments -----------------------------
compute_reach_estimated_segments <- function(dataset) {
  responses_per_segment <- dataset %>%
    group_by(estimated_gender) %>%
    summarise(count = sum(response), .groups = 'drop') %>%
    group_by(estimated_gender) %>%
    mutate(probability = count / sum(count))
  
  # Calculate gender flags
  female <- ifelse(dataset$estimated_gender == 'female', 1, 0)
  male <- ifelse(dataset$estimated_gender == 'male', 1, 0)
  est_seg <- bind_cols(female, male)
  
  # Compute counts for each gender
  count_female <- sum(female)
  count_male <- sum(male)
  
  # Compute reach
  # true_female <- count_female * p_is[1, 1] + count_male * p_is[1, 2]
  # true_male <- count_male * p_is[2, 2] + count_male * p_is[2, 1]
  reach_female <- responses_per_segment$count[1] / count_female
  reach_male <- responses_per_segment$count[2] / count_male
  
  print(paste("Reach male:", reach_male, ". Reach female:", reach_female))
}

# ----------------------------- Compute Conditional Probabilities including prior (Z|S) -----------------------------
compute_p_z_given_s_including_prior <- function(dataset, segmentation) {
  # Collect universe true priors
  universe_estimates <- read_universe_estimates()
  
  # Clean missing true data
  dataset <- dataset %>% filter(complete.cases(.))
  
  # Calculate conditional probabilities
  if (segmentation == "gender") {
    p_z <- universe_estimates %>%
      group_by(gender_bucket) %>%
      summarise(count = sum(num_persons)) %>%
      mutate(probability = count / (universe_estimates$tot_persons[1]))
    p_s_given_z <- dataset %>%
      group_by(estimated_gender, true_gender) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(true_gender) %>%
      mutate(probability = count / sum(count))
  } else if (segmentation == "age") {
    p_z <- universe_estimates %>%
      group_by(age_bucket) %>%
      summarise(count = sum(num_persons)) %>%
      mutate(probability = count / (universe_estimates$tot_persons[1]))
    p_s_given_z <- dataset %>%
      group_by(estimated_age, true_age) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(true_age) %>%
      mutate(probability = count / sum(count))
  } else if (segmentation == "demo") {
    p_z <- universe_estimates %>%
      group_by(demo3_bucket) %>%
      summarise(count = sum(num_persons)) %>%
      mutate(probability = count / (universe_estimates$tot_persons[1]))
    p_s_given_z <- dataset %>%
      group_by(estimated_demo, true_demo) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(true_demo) %>%
      mutate(probability = count / sum(count))
  }
  
  
  dim_matrix <- sqrt(length(p_s_given_z$probability))
  p_s_given_z_matrix <- matrix(p_s_given_z$probability, nrow = dim_matrix, ncol = dim_matrix)
  
  p_s_and_z <- p_s_given_z_matrix * p_z$probability
  p_s <- colSums(p_s_and_z)
  
  p_s_matrix <- t(matrix(rep((p_s**-1),dim_matrix), ncol = dim_matrix))
  p_z_given_s <- t(p_s_and_z*p_s_matrix)
  
  print(p_z_given_s)
  return(p_z_given_s)
}


# ----------------------------- Compute Conditional Probabilities P(Z|S) -----------------------------
compute_p_z_given_s <- function(dataset, segmentation) {
  # Clean missing true data
  dataset <- dataset %>% filter(complete.cases(.))
  
  # Calculate conditional probabilities
  if (segmentation == "gender") {
    conditional_probs <- dataset %>%
      group_by(estimated_gender, true_gender) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(estimated_gender) %>%
      mutate(probability = count / sum(count))
  } else if (segmentation == "age") {
    conditional_probs <- dataset %>%
      group_by(estimated_age, true_age) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(estimated_age) %>%
      mutate(probability = count / sum(count))
  } else if (segmentation == "demo") {
    conditional_probs <- dataset %>%
      group_by(estimated_demo, true_demo) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(estimated_demo) %>%
      mutate(probability = count / sum(count))
  }
  
  dim_matrix <- sqrt(length(conditional_probs$probability))
  conditional_probs_matrix <- matrix(conditional_probs$probability, nrow = dim_matrix, ncol = dim_matrix)
  
  p_z_given_s <- t(conditional_probs_matrix)
  
  print(p_z_given_s)
  return(p_z_given_s)
}

# # ----------------------------- Compute Log Likelihood -----------------------------
# loglikelihood <- function(beta) {
#   print(beta)
#   log_value <- 0
#   
#   for (i in 1:n_test) {
#     # if (male[i] == 1) {
#     #   test_exposure_data[i] <- 0
#     # }
#     log_value_between <- 0
#     
#     for (s in 1:2) {
#       e_z_beta <- exp(beta[s])
#       log_value_between <- log_value_between + expect_indicator[i, s] * 
#         (((e_z_beta / (1 + e_z_beta))^test_exposure_data[i]) * 
#            ((1 / (1 + e_z_beta))^(1 - test_exposure_data[i])))
#     }
#     
#     log_value <- log_value + log(log_value_between)
#   }
#   
#   print(log_value)
#   return(-log_value)
# }

# ----------------------------- Compute Segment Sizes -----------------------------
compute_segment_sizes <- function(dataset, segmentation) {
  dataset <- dataset %>% filter(is.na(true_gender) | is.na(true_age) | is.na(true_demo))
  # Summarize response and non-response counts
  if (segmentation == "gender") {
    segments_response <- dataset %>%
      group_by(estimated_gender) %>%
      summarise(
        response_count = sum(response),
        no_response_count = sum(1 - response),
        .groups = 'drop'
      )
  } else if (segmentation == "age") {
    segments_response <- dataset %>%
      group_by(estimated_age) %>%
      summarise(
        response_count = sum(response),
        no_response_count = sum(1 - response),
        .groups = 'drop'
      )
  } else if (segmentation == "demo") {
    segments_response <- dataset %>%
      group_by(estimated_demo) %>%
      summarise(
        response_count = sum(response),
        no_response_count = sum(1 - response),
        .groups = 'drop'
      )
  }
  
  mat_segments_response <- matrix(
    c(segments_response$response_count, segments_response$no_response_count),
    ncol = 2
  )
  
  print(segments_response)
  
  return(mat_segments_response)
}

#--------------------------------Compute True Segment Sizes-------------------------------------
compute_true_segment_sizes<-function(dataset,segmentation){
  #Delete the observations with NA values
  dataset <- dataset %>% filter(complete.cases(.))
  # Summarize response and non-response counts
  if (segmentation == "gender") {
    segments_response <- dataset %>%
      group_by(true_gender) %>%
      summarise(
        response_count = sum(response),
        no_response_count = sum(1 - response),
        .groups = 'drop'
      )
  } else if (segmentation == "age") {
    segments_response <- dataset %>%
      group_by(true_age) %>%
      summarise(
        response_count = sum(response),
        no_response_count = sum(1 - response),
        .groups = 'drop'
      )
  } else if (segmentation == "demo") {
    segments_response <- dataset %>%
      group_by(true_demo) %>%
      summarise(
        response_count = sum(response),
        no_response_count = sum(1 - response),
        .groups = 'drop'
      )
  }
  
  mat_segments_response <- matrix(
    c(segments_response$response_count, segments_response$no_response_count),
    ncol = 2
  )
  
  print(segments_response)
  
  return(mat_segments_response)
}
# 
# ----------------------------- Compute Log Likelihood for Segments -----------------------------
loglikelihood_segments_based <- function(beta, p_z_given_s, segment_responses, true_segments) {
  log_value_est <- 0
  log_value_true <- 0
  segment_count <- dim(p_z_given_s)[1]
  
  # Compute likelihood for estimated demographics
  for (i in 1:segment_count) {
    for (j in 1:2) {  # j=1 (response), j=2 (no response)
      response_j <- ifelse(j == 1, 1, 0)
      value_segment <- 0
      
      for (s in 1:segment_count) {
        e_z_beta <- plogis(beta[s])
        value_segment <- value_segment + p_z_given_s[i, s] *
          (((e_z_beta^response_j) * ((1 - e_z_beta)^(1 - response_j))))
      }
      
      log_value_est <- log_value_est + segment_responses[i, j] * log(max(value_segment, 1e-10))
    }
  }
  
  # Compute likelihood for true demographics
  for (z in 1:segment_count) {
    value_segment <- 0  
    for (j in 1:2) {
      response_j <- ifelse(j == 1, 1, 0)
      e_z_beta <- plogis(beta[z])
      value_segment <- value_segment + (((e_z_beta^response_j) * ((1 - e_z_beta)^(1 - response_j))))
    }
    
    log_value_true <- log_value_true + true_segments[z, j] * log(max(value_segment, 1e-10))
  }
  
  # # 🔹 **Regularization term to prevent extreme beta values**
  # lambda <- 0.1   # Adjust this value to control regularization strength
  # reg_term <- lambda * sum(beta^2)
  
  return(-log_value_true - log_value_est )
}



# ----------- Compute Fraction of Response and Compare to Prediction --------------------
evaluation <- function(dataset, beta, segmentation) {
  # Compute predicted probabilities
  pred_probabilities <- exp(beta) / (1 + exp(beta))
  
  if (segmentation == "gender") {
    # Fraction of response for estimated gender segments
    response_for_estimated_seg <- dataset %>%
      group_by(estimated_gender) %>%
      summarise(
        fraction_response = mean(response, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(predicted_probability = case_when(
        estimated_gender == "female" ~ pred_probabilities[1],
        estimated_gender == "male" ~ pred_probabilities[2],
        TRUE ~ NA_real_
      ))
    
    # Fraction of response for true gender segments
    response_for_true_seg <- dataset %>%
      filter(!is.na(true_gender)) %>%
      group_by(true_gender) %>%
      summarise(
        fraction_response = mean(response, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(predicted_probability = case_when(
        true_gender == "female" ~ pred_probabilities[1],
        true_gender == "male" ~ pred_probabilities[2],
        TRUE ~ NA_real_
      ))
  } else if (segmentation == "age") {
    # Fraction of response for estimated age segments
    response_for_estimated_seg <- dataset %>%
      group_by(estimated_age) %>%
      summarise(
        fraction_response = mean(response, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(predicted_probability = pred_probabilities)
    
    # Fraction of response for true age segments
    response_for_true_seg <- dataset %>%
      filter(!is.na(true_age)) %>%
      group_by(true_age) %>%
      summarise(
        fraction_response = mean(response, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(predicted_probability = pred_probabilities)
  }
  
  return(list(
    response_for_true_seg = response_for_true_seg,
    response_for_estimated_seg = response_for_estimated_seg
  ))
}


