# # ----------------------------- Compute Reach for Estimated Segments -----------------------------
# compute_reach_estimated_segments <- function(dataset) {
#   responses_per_segment <- dataset %>%
#     group_by(estimated_gender) %>%
#     summarise(count = sum(response), .groups = 'drop') %>%
#     group_by(estimated_gender) %>%
#     mutate(probability = count / sum(count))
#   
#   # Calculate gender flags
#   female <- ifelse(dataset$estimated_gender == 'female', 1, 0)
#   male <- ifelse(dataset$estimated_gender == 'male', 1, 0)
#   est_seg <- bind_cols(female, male)
#   
#   # Compute counts for each gender
#   count_female <- sum(female)
#   count_male <- sum(male)
#   
#   # Compute reach
#   # true_female <- count_female * p_is[1, 1] + count_male * p_is[1, 2]
#   # true_male <- count_male * p_is[2, 2] + count_male * p_is[2, 1]
#   reach_female <- responses_per_segment$count[1] / count_female
#   reach_male <- responses_per_segment$count[2] / count_male
#   
#   print(paste("Reach male:", reach_male, ". Reach female:", reach_female))
# }

# ----------------------------- Compute Conditional Probabilities including prior (Z|S) -----------------------------
compute_p_z_given_s_including_prior <- function(dataset, segmentation) {
  
  universe_estimates <- read_universe_estimates()
  print(universe_estimates)
  if (segmentation == "gender") {
    # Clean missing true data
    dataset <- dataset %>% filter(!is.na(true_gender))
    
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
    
    dataset <- dataset %>% filter(!is.na(true_age))
    
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
    
    dataset <- dataset %>% filter(!is.na(true_demo))
    
    p_z <- universe_estimates %>%
      group_by(demo3_bucket) %>%
      summarise(count = sum(num_persons)) %>%
      mutate(probability = count / (universe_estimates$tot_persons[1]))
    p_s_given_z <- dataset %>%
      group_by(estimated_demo, true_demo) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(true_demo) %>%
      mutate(probability = count / sum(count))
  }else if (segmentation == "full") {
    dataset <- dataset %>% filter(!is.na(true_demo) & !is.na(true_age) & !is.na(true_gender))
    
    # Generate all possible (true, estimated) pairs
    all_combinations <- expand.grid(
      estimated_gender = 0:1,
      estimated_age = 1:4,
      estimated_demo = 1:5,
      true_gender = 0:1,
      true_age = 1:4,
      true_demo = 1:5
    )
    
    # Compute p_z only for observed true demographic values
    observed_true_demos <- dataset %>% select(true_gender, true_age, true_demo) %>% distinct()
    
    p_z <- universe_estimates %>%
      mutate(probability = num_persons / tot_persons)
    
    # Compute p_s_given_z based on dataset counts
    p_s_given_z <- dataset %>%
      group_by(estimated_gender, estimated_age, estimated_demo, true_gender, true_age, true_demo) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(true_gender, true_age, true_demo) %>%
      mutate(probability = count / sum(count))
    
    # Ensure all 40×40 combinations exist
    p_s_given_z <- all_combinations %>%
      left_join(p_s_given_z, 
                by = c("estimated_gender", "estimated_age", "estimated_demo", "true_gender", "true_age", "true_demo")) %>%
      replace_na(list(count = 0, probability = NA))  # Set missing probabilities to NA first
    
    # Fill missing true demos with uniform probability across all estimated values
    p_s_given_z <- p_s_given_z %>%
      group_by(true_gender, true_age, true_demo) %>%
      mutate(probability = ifelse(all(is.na(probability)), 1/n(), probability)) %>%
      ungroup()
    print(p_z)
    print(p_s_given_z)
  } 
  
  
  dim_matrix <- sqrt(length(p_s_given_z$probability))
  p_s_given_z_matrix <- matrix(p_s_given_z$probability, nrow = dim_matrix, ncol = dim_matrix)
  
  p_s_and_z <- p_s_given_z_matrix * p_z$probability
  p_s <- colSums(p_s_and_z)
  
  p_s_matrix <- t(matrix(rep((p_s**-1),dim_matrix), ncol = dim_matrix))
  p_z_given_s <- t(p_s_and_z*p_s_matrix)
  
  # print(p_z_given_s)
  return(p_z_given_s)
}


# ----------------------------- Compute Conditional Probabilities P(Z|S) -----------------------------
compute_p_z_given_s <- function(dataset, segmentation) {
  # Calculate conditional probabilities
  if (segmentation == "gender") {
    # Clean missing true data
    dataset <- dataset %>% filter(!is.na(true_gender))
    
    conditional_probs <- dataset %>%
      group_by(estimated_gender, true_gender) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(estimated_gender) %>%
      mutate(probability = count / sum(count))
  } else if (segmentation == "age") {
    # Clean missing true data
    dataset <- dataset %>% filter(!is.na(true_age))
    
    conditional_probs <- dataset %>%
      group_by(estimated_age, true_age) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(estimated_age) %>%
      mutate(probability = count / sum(count))
  } else if (segmentation == "demo") {
    # Clean missing true data
    dataset <- dataset %>% filter(!is.na(true_demo))
    
    conditional_probs <- dataset %>%
      group_by(estimated_demo, true_demo) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(estimated_demo) %>%
      mutate(probability = count / sum(count))
  }
  
  dim_matrix <- sqrt(length(conditional_probs$probability))
  conditional_probs_matrix <- matrix(conditional_probs$probability, nrow = dim_matrix, ncol = dim_matrix)
  
  p_z_given_s <- t(conditional_probs_matrix)
  
  # print(p_z_given_s)
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
compute_segment_sizes <- function(dataset, segmentation, use_true_seperate = FALSE) {
  if (use_true_seperate) {
    if (segmentation == "gender") {
      dataset_true_only <- dataset %>% filter(!is.na(true_gender))
      
      segments_response_true <- dataset_true_only %>%
        group_by(true_gender) %>%
        summarise(
          response_count = sum(response),
          no_response_count = sum(1 - response),
          .groups = 'drop'
        )
      
      dataset <- dataset %>% filter(is.na(true_gender))
    } else if (segmentation == "age") {
      dataset_true_only <- dataset %>% filter(!is.na(true_age))
      
      segments_response_true <- dataset_true_only %>%
        group_by(true_age) %>%
        summarise(
          response_count = sum(response),
          no_response_count = sum(1 - response),
          .groups = 'drop'
        )
      dataset <- dataset %>% filter(is.na(true_age))
    } else if (segmentation == "demo") {
      dataset_true_only <- dataset %>% filter(!is.na(true_demo))
      
      segments_response_true <- dataset_true_only %>%
        group_by(true_demo) %>%
        summarise(
          response_count = sum(response),
          no_response_count = sum(1 - response),
          .groups = 'drop'
        )
      dataset <- dataset %>% filter(is.na(true_demo))
    } else if (segmentation=='full'){
      dataset_true_only <-dataset%>% filter(!is.na(true_age)|!is.na(true_gender)|!is.na(true_demo))
      
      segments_response_true<- dataset_true_only %>%
        group_by(true_gender,true_age,true_demo) %>%
        summarise(
          response_count=sum(response),
          no_response_count=sum(1-response),
          .groups='drop'
        )
    }
    
    mat_segments_true_response <- matrix(
      c(segments_response_true$response_count, segments_response_true$no_response_count),
      ncol = 2
    )
  } else {
    mat_segments_true_response <- c(0,0)
  }
  print(mat_segments_true_response)
  
  dataset<-dataset%>%filter(is.na(true_age) | is.na(true_gender) | is.na(true_demo))
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
  } else if (segmentation=="full"){
    segments_response<- dataset%>%
      group_by(estimated_gender,estimated_age,estimated_demo)%>%
      summarise(
        response_count=sum(response),
        no_response_count= sum(1-response),
        .groups='drop'
      )
    
    segments_response_true
  }
  
  mat_segments_response <- matrix(
    c(segments_response$response_count, segments_response$no_response_count),
    ncol = 2
  )
  
  print(segments_response)
  
  return(list(mat_segments_response, mat_segments_true_response))
}

# 
# ----------------------------- Compute Log Likelihood for Segments -----------------------------
loglikelihood_segments_based <- function(beta, p_z_given_s, segment_responses, true_value_weights = c(1,1)) {
  log_value_est <- 0
  log_value_true <- 0
  segment_count <- dim(p_z_given_s)[1]
  
  for (i in 1:segment_count) {
    for (j in 1:2) {
      response_j <- ifelse(j == 1, 1, 0)
      value_segment <- 0
      
      for (s in 1:segment_count) {
        e_z_beta <- exp(beta[s])
        value_segment <- value_segment + p_z_given_s[i, s] *
          (((e_z_beta / (1 + e_z_beta))^response_j) *
             ((1 / (1 + e_z_beta))^(1 - response_j)))
      }
      
      log_value_est <- log_value_est + segment_responses[[1]][i, j] * log(value_segment)
    }
  }
  if (sum(segment_responses[[2]]) != 0) {
    for (z in 1:segment_count) {
      for (j in 1:2) {
        response_j <- ifelse(j == 1, 1, 0)
        e_z_beta <- exp(beta[z]) # Convert beta to probability
        
        value_segment <- ((e_z_beta / (1 + e_z_beta))^response_j) * ((1 / (1 + e_z_beta))^(1 - response_j))
        
        # Log likelihood contribution for this (y_j, z) pair
        log_value_true <- log_value_true + segment_responses[[2]][z, j] * log(value_segment)
      }
    }
  }
  
  return(-1*(true_value_weights[1]*log_value_est + true_value_weights[2]*log_value_true))
}

# ----------------------------- Optimize Log Likelihood -----------------------------
optimize_loglikelihood <- function(dataset, segmentation, with_prior, use_true_seperate = FALSE, print_result = TRUE) {
  if (with_prior) {
    p_z_given_s <- compute_p_z_given_s_including_prior(dataset, segmentation)
  } else {
    p_z_given_s <- compute_p_z_given_s(dataset, segmentation)
  }
  response_segments <- compute_segment_sizes(dataset, segmentation, use_true_seperate)
  
  segment_count <- dim(p_z_given_s)[1]
  initial_par <- rep(0, segment_count)
  
  best_result <- optim(
    par = initial_par,
    fn = loglikelihood_segments_based,
    p_z_given_s = p_z_given_s,
    segment_responses = response_segments,
    method = "L-BFGS-B",     
    lower = c(-10, -10),     
    upper = c(0, 0)
  )
  if (print_result) {
    if (segmentation == "gender") {
      reach_female <- exp(best_result$par[1]) / (1 + exp(best_result$par[1]))
      reach_male <- exp(best_result$par[2]) / (1 + exp(best_result$par[2]))
      cat(paste(
        "Reach of true segments are:",
        "\nFor male:", reach_male,
        "\nFor female:", reach_female,
        "\nWith beta:", best_result$par[1], best_result$par[2]
      ))
    } else if (segmentation == "age") {
      reach <- sapply(1:4, function(i) exp(best_result$par[i]) / (1 + exp(best_result$par[i])))
      cat(paste("Reach for age segments:\n", paste(reach, collapse = "\n")))
    } else if (segmentation == "demo") {
      reach <- sapply(1:5, function(i) exp(best_result$par[i]) / (1 + exp(best_result$par[i])))
      cat(paste("Reach for demo3 segments:\n", paste(reach, collapse = "\n")))
    }
  }
  
  return(best_result)
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
        estimated_gender == "0" ~ pred_probabilities[1],
        estimated_gender == "1" ~ pred_probabilities[2],
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
        true_gender == "0" ~ pred_probabilities[1],
        true_gender == "1" ~ pred_probabilities[2],
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
  } else if (segmentation == "demo") {
    # Fraction of response for estimated demo segments
    response_for_estimated_seg <- dataset %>%
      group_by(estimated_demo) %>%
      summarise(
        fraction_response = mean(response, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(predicted_probability = pred_probabilities)
    
    # Fraction of response for true age segments
    response_for_true_seg <- dataset %>%
      filter(!is.na(true_demo)) %>%
      group_by(true_demo) %>%
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


#------------ Compute True Frequency --------------------------------
compute_true_frequency <- function(dataset, segmentation) {
  
  # Map user-friendly segmentation input to the correct column names
  column_mapping <- list(
    "gender" = "estimated_gender",
    "age" = "estimated_age",
    "demo" = "estimated_demo"
  )

  # Get the actual column name from mapping
  segmentation_col <- column_mapping[[segmentation]]
  
  # Remove rows without valid exposures
  dataset <- dataset %>%
    filter(!is.na(total_exposures) & total_exposures > 0)
  
  # Sum exposures per segment
  total_exposures_segment <- dataset %>%
    group_by(!!sym(segmentation_col)) %>%
    summarise(total_exposures = sum(total_exposures, na.rm = TRUE), .groups = 'drop')
  
  # Find reach per segment (unique users)
  reach_segment <- dataset %>%
    group_by(!!sym(segmentation_col)) %>%
    summarise(reach = n_distinct(person_id), .groups = 'drop')
  
  # Merge and calculate frequency
  frequency_data <- merge(total_exposures_segment, reach_segment, by = segmentation_col)
  frequency_data <- frequency_data %>%
    mutate(frequency = total_exposures / reach)
  
  # Print and return results
  print(frequency_data)
  return(frequency_data)
}


