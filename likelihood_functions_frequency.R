# ------------- Compute number of observations per variable combination ----------

compute_segment_sizes_frequency <- function(dataset, segmentation, use_true_seperate = TRUE, simulation = FALSE) {
  true_col <- paste0("true_", segmentation)
  estimated_col <- paste0("estimated_", segmentation)
  
  if (simulation) {
    segments_exposures_true <- tibble()
    
    segments_exposures <- dataset %>%
      group_by(!!sym(estimated_col), total_exposures) %>%
      summarise(
        num_individuals = n(),
        .groups = 'drop'
      )
  }
  else {
    if (use_true_seperate) {
      dataset_true_only <- dataset[[1]]
      
      segments_exposures_true <- dataset_true_only %>%
        group_by(!!sym(true_col), total_exposures) %>%
        summarise(
          num_individuals = sum(individual_count),  # Count number of people with this exposure level
          .groups = 'drop'
        ) %>% filter(total_exposures > 0)
    } else {
      segments_exposures_true <- tibble()  # Empty tibble if not using true segmentation
    }
    
    # Compute estimated segment exposures
    segments_exposures <- dataset[[2]] %>%
      group_by(!!sym(estimated_col), total_exposures) %>%
      summarise(
        num_individuals = sum(individual_count),
        .groups = 'drop'
      ) %>% filter(total_exposures > 0)
  }
  
  return(list(segments_exposures, segments_exposures_true))
}


# ----------------------------- Log-Likelihood Function Poisson -----------------------------

loglikelihood_segments_based_frequency <- function(beta, p_z_given_s, segment_exposures, segmentation, use_binomial = FALSE) {
  log_value <- 0  # Initialize log-likelihood
  segment_levels <- unique(segment_exposures[[1]][[paste0("estimated_", segmentation)]])  # Get unique segment values
  estimated_exposures <- segment_exposures[[1]]
  true_exposures <- segment_exposures[[2]]
  
  
  # Compute log-likelihood for estimated segments weighted by P(Z | S)
  for (i in segment_levels) {
    lambda_i <- exp(beta[i])  # Poisson rate for estimated segment i
    
    # Get exposures for the current estimated segment
    exposure_subset <- estimated_exposures %>% filter(!!sym(paste0("estimated_", segmentation)) == i)
    
    for (j in 1:nrow(exposure_subset)) {
      num_exposed <- exposure_subset$num_individuals[j]
      observed_exposure <- log(exposure_subset$total_exposures[j])
      
      # Sum over all possible true segments Z, weighted by P(Z | S)
      weighted_likelihood <- 0
      if (use_binomial) {
        theta_mat <- matrix(beta, ncol = 2)
        for (z in unique(estimated_exposures[[paste0("estimated_", segmentation)]])) {
          p_z_given_s_val <- p_z_given_s[i, z]  # Transition probability P(Z | S)
          r_z <- theta_mat[z,1]
          p_z <- theta_mat[z,2]
          
          weighted_likelihood <- weighted_likelihood + p_z_given_s_val * (gamma(observed_exposure + r_z)/
                                                          (factorial(observed_exposure)*gamma(r_z))*
                                                          p_z**r_z*(1-p_z)**observed_exposure)
        }
      } else {
        for (z in unique(estimated_exposures[[paste0("estimated_", segmentation)]])) {
          p_z_given_s_val <- p_z_given_s[i, z]  # Transition probability P(Z | S)
          
          weighted_likelihood <- weighted_likelihood + 
            p_z_given_s_val * (exp(-1 * beta[z])/factorial(observed_exposure)*beta[z]**observed_exposure)
        }
      }
      
      if (observed_exposure < 7) {
        log_value <- log_value + num_exposed * log(weighted_likelihood)
      }
    }
    if (sum(true_exposures) > 0) {
      # Get exposures for true_segments
      exposure_true_subset <- true_exposures %>% filter(!!sym(paste0("true_", segmentation)) == i)
      
      for (j in 1:nrow(exposure_true_subset)) {
        num_exposed <- exposure_true_subset$num_individuals[j]
        observed_exposure <- log(exposure_true_subset$total_exposures[j])
        
        # Sum over all possible true segments Z, weighted by P(Z | S)
        if (use_binomial) {
          theta_mat <- matrix(beta, ncol = 2)
          r_z <- theta_mat[i,1]
          p_z <- theta_mat[i,2]
          weighted_likelihood <- (gamma(observed_exposure + r_z)/ 
                                    (factorial(observed_exposure)*gamma(r_z))* 
                                    p_z**r_z*(1-p_z)**observed_exposure)
        } else {
            weighted_likelihood <- 
              (exp(-1 * beta[z])/factorial(observed_exposure)*beta[z]**observed_exposure)
          }
        
        if (observed_exposure < 7) {
          log_value <- log_value + num_exposed * log(weighted_likelihood)
        }
      }
    }
  }
  
  return(-1 * log_value)  # Return negative log-likelihood for optimization
}


# ----------------------------- Optimize Poisson Log-Likelihood -----------------------------

optimize_loglikelihood_frequency <- function(dataset, segmentation, with_prior = TRUE, 
                                             print_result = TRUE, use_binomial = FALSE, simulation = FALSE) {
  # Collect P(Z | S)
  if (with_prior) {
    p_z_given_s <- get(paste0("p_z_given_s_", segmentation, "_with_prior"))
  } else {
    p_z_given_s <- get(paste0("p_z_given_s_", segmentation, "_without_prior"))
  }
  
  segment_exposures <- compute_segment_sizes_frequency(dataset, segmentation, simulation = simulation)
  
  segmentation_col <- paste0("estimated_", segmentation)
  
  if (!simulation) {
    # Compute mean exposures for better initial values
    exposure_means <- segment_exposures[[1]] %>%
      group_by(!!sym(segmentation_col)) %>%
      summarise(full_total_exposure = sum(total_exposures*num_individuals), 
                total_count = sum(num_individuals),
                .groups = "drop") %>%
      mutate(mean_exposure = full_total_exposure/total_count)
  } else {
    exposure_means <- segment_exposures[[1]] %>%
      group_by(!!sym(segmentation_col)) %>%
      summarise(full_total_exposure = sum(total_exposures*num_individuals), 
                total_count = sum(num_individuals),
                .groups = "drop") %>%
      mutate(mean_exposure = full_total_exposure/total_count)
  }
  
  print(exposure_means$mean_exposure)
  
  segment_count <- length(exposure_means$mean_exposure)
  # Prevent log(0) errors
  if (use_binomial) {
    initial_par <- c(rep(0.5,2*segment_count))
    
    best_result <- optim(
      par = initial_par,
      fn = loglikelihood_segments_based_frequency,
      p_z_given_s = p_z_given_s,
      segment_exposures = segment_exposures,
      segmentation = segmentation,
      use_binomial = TRUE,
      method = "L-BFGS-B",
      lower = c(rep(0.1, segment_count), rep(0.01, segment_count)),     
      upper =  c(rep(50, segment_count), rep(0.90, segment_count))
    )
    beta_final <- exp(best_result$par[1:segment_count]*(1-best_result$par[(segment_count+1):(2*segment_count)])/ best_result$par[(segment_count+1):(2*segment_count)])
    
  } else {
    initial_beta <- pmax(exposure_means$mean_exposure, 0.01)
    initial_par <- rep(0.5, segment_count)
    print(initial_par)
    best_result <- optim(
      par = initial_par,
      fn = loglikelihood_segments_based_frequency,
      p_z_given_s = p_z_given_s,
      segment_exposures = segment_exposures,
      segmentation = segmentation,
      method = "L-BFGS-B",
      lower = c(0.1),     
      upper = c(10)
    )
    beta_final <- exp(best_result$par)  
  }
  
  print(initial_par)
  
  if (print_result) {
    cat("\n✅ Optimization converged successfully!\n")
    cat("\nOptimized Model Parameters:")
    cat("\nBeta:", best_result$par, "\n")
    cat("\nPredicted Exposures:", beta_final, "\n")  
  }
  
  return(list(beta = beta_final, predicted_exposures = beta_final, optimization = best_result))
}

# ----------------------------- Compute Expected Frequency Given Reach -----------------------------

evaluation_poisson <- function(dataset, beta, segmentation) {
  segmentation_col_est <- paste0("estimated_", segmentation)
  segmentation_col_true <- paste0("true_", segmentation)
  
  pred_expected_exposures <- beta
  
  dataset <- dataset %>%
    filter(total_exposures > 0)

  
  # Compute mean exposures for estimated segments
  exposure_for_estimated_seg <- dataset %>%
    group_by(!!sym(segmentation_col_est)) %>%
    summarise(
      mean_exposures = mean(total_exposures, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(predicted_exposures = pred_expected_exposures)
  
  # Compute mean exposures for true segments
  exposure_for_true_seg <- dataset %>%
    filter(!is.na(!!sym(segmentation_col_true))) %>%
    group_by(!!sym(segmentation_col_true)) %>%
    summarise(
      mean_exposures = mean(total_exposures, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("\n🔍 Estimated Segments - Mean & Predicted Exposures:\n")
  print(exposure_for_estimated_seg)
  
  cat("\n🔍 True Segments - Mean Exposures:\n")
  print(exposure_for_true_seg)
  
  return(list(estimated_exposures = exposure_for_estimated_seg, true_exposures = exposure_for_true_seg))
}

# # ----------------------------- Run the Poisson Model -----------------------------
# 
# segmentation <- "age"  # Choose "gender", "age", or "demo"
# 
# # ----------------------- For simulated data ----------------------------
# 
# n <- 10000000
# dataset_frequency <- simulate_exposure_dataset(n, c(0.516, 0.484), c(1.8, 1.2), matrix(c(0.7, 0.3, 0.2, 0.8), ncol = 2),
#                                                segmentation = "gender", dispersion = 1, distribution = "poisson",
#                                                seed_number = 0, print_simulation = TRUE)
# dataset_frequency <- simulate_exposure_dataset(
#   n_test = 1e7,
#   fraction_per_segment = c(0.25, 0.25, 0.25, 0.25),
#   exposure_per_segment = c(3, 1.4, 1.2, 1.8),
#   estimation_correctness = matrix(
#     c(0.7, 0.1, 0.1, 0.1,
#       0.05, 0.75, 0.1, 0.1,
#       0.1, 0.1, 0.7, 0.1,
#       0.1, 0.1, 0.1, 0.7),
#     ncol = 4, byrow = TRUE),
#   segmentation = "age",
#   dispersion = 1,
#   distribution = "poisson",
#   seed_number = 0,
#   print_simulation = TRUE
# )
# 
# dataset_frequency_filtered <- dataset_frequency %>%
#   filter(total_exposures > 0)
# 
# 
# pois_results <- optimize_loglikelihood_frequency(dataset_frequency_filtered, segmentation, with_prior = FALSE, simulation = TRUE)
# eval_results <- evaluation_poisson(dataset_frequency, pois_results$beta, segmentation)
# print(eval_results)
# 
# # ----------------------- For real data ----------------------------
# 
# real_dataset <- read_exposures()
# real_dataset_filtered <- real_dataset
# 
# pois_results <- optimize_loglikelihood_poisson(real_dataset, segmentation, with_prior = TRUE)
# eval_results <- evaluation_poisson(real_dataset_filtered, pois_results$beta, segmentation)
# print(eval_results)