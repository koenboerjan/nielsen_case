
#----------------------------- Compute Conditional Probabilities including prior (Z|S) -----------------------------
compute_p_z_given_s_including_prior <- function(dataset, segmentation) {
  # Collect universe true priors
  universe_estimates <- read_universe_estimates()

  # Calculate conditional probabilities
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
    # Clean missing true data
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
    # Clean missing true data
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
    }
    
    mat_segments_true_response <- matrix(
      c(segments_response_true$response_count, segments_response_true$no_response_count),
      ncol = 2
    )
  } else {
    mat_segments_true_response <- c(0,0)
  }
  
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
  
  # print(segments_response)
  
  return(list(mat_segments_response, mat_segments_true_response))
}

# 

# ----------------------------- Compute Segment Sizes shortcut real dataset -----------------------------
compute_segment_sizes_real_dataset <- function(dataset, segmentation) {
  if (segmentation == "gender") {
    segments_response_true <- dataset[[1]] %>%
      group_by(true_gender) %>%
      summarise(
        response_count = sum(response_count),
        no_response_count = sum(no_response_count),
        .groups = 'drop'
        )
    segments_response <- dataset[[2]] %>%
      group_by(estimated_gender) %>%
      summarise(
        response_count = sum(response_count),
        no_response_count = sum(no_response_count),
        .groups = 'drop'
        )
  } else if (segmentation == "age") {
    segments_response_true <- dataset[[1]] %>%
      group_by(true_age) %>%
      summarise(
        response_count = sum(response_count),
        no_response_count = sum(no_response_count),
        .groups = 'drop'
      )
    segments_response <- dataset[[2]] %>%
      group_by(estimated_age) %>%
      summarise(
        response_count = sum(response_count),
        no_response_count = sum(no_response_count),
        .groups = 'drop')
  } else if (segmentation == "demo") {
    segments_response_true <- dataset[[1]] %>%
      group_by(true_demo) %>%
      summarise(
        response_count = sum(response_count),
        no_response_count = sum(no_response_count),
        .groups = 'drop'
      )
    segments_response <- dataset[[2]] %>%
      group_by(estimated_demo) %>%
      summarise(
        response_count = sum(response_count),
        no_response_count = sum(no_response_count),
        .groups = 'drop')
  }
  
  mat_segments_true_response <- matrix(
    c(segments_response_true$response_count, segments_response_true$no_response_count),
    ncol = 2
  )
  
  mat_segments_response <- matrix(
    c(segments_response$response_count, segments_response$no_response_count),
    ncol = 2
  )
  
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
optimize_loglikelihood <- function(dataset, segmentation, with_prior, use_true_seperate = FALSE, print_result = TRUE, simulation = FALSE) {
    # Compute P(Z | S)
  if (simulation) {
    if (with_prior) {
      p_z_given_s <- compute_p_z_given_s_including_prior(dataset, segmentation)
    } else {
      p_z_given_s <- compute_p_z_given_s(dataset, segmentation)
    }
  } else {
    if (with_prior) {
      p_z_given_s <- get(paste0("p_z_given_s_", segmentation, "_with_prior"))
    } else {
      p_z_given_s <- get(paste0("p_z_given_s_", segmentation, "_without_prior"))
    }
  }
  
  print(p_z_given_s)
  
  if (simulation) {
    response_segments <- compute_segment_sizes(dataset, segmentation)
  } else {
    response_segments <- compute_segment_sizes_real_dataset(dataset, segmentation)
  }
  
  
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


