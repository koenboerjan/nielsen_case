# ------------- Compute number of observations per variable combination ----------

compute_segment_sizes_poisson <- function(dataset, segmentation, use_true_seperate = TRUE) {
  true_col <- paste0("true_", segmentation)
  estimated_col <- paste0("estimated_", segmentation)
  
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
  
  print(segments_exposures)
  return(list(segments_exposures, segments_exposures_true))
}

# ----------------------------- Log-Likelihood Function Poisson -----------------------------

loglikelihood_segments_based_poisson <- function(beta, p_z_given_s, segment_exposures, segmentation) {
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
      observed_exposure <- exposure_subset$total_exposures[j]
      
      # Sum over all possible true segments Z, weighted by P(Z | S)
      weighted_likelihood <- 0
      for (z in unique(estimated_exposures[[paste0("estimated_", segmentation)]])) {
        p_z_given_s_val <- p_z_given_s[i, z]  # Transition probability P(Z | S)
        
        weighted_likelihood <- weighted_likelihood + 
          p_z_given_s_val * (exp(-1 * beta[z])/factorial(observed_exposure)*beta[z]**observed_exposure)
      }
      log_value <- log_value + num_exposed * log(weighted_likelihood)

      if (observed_exposure < 500) {
        log_value <- log_value + num_exposed * log(weighted_likelihood)
      }
    }
  }
  
  return(-1 * log_value)  # Return negative log-likelihood for optimization
}


# ----------------------------- Optimize Poisson Log-Likelihood -----------------------------

optimize_loglikelihood_poisson <- function(dataset, segmentation, with_prior = TRUE, print_result = TRUE) {
  segmentation_col <- paste0("estimated_", segmentation)
  
  # if (!segmentation_col %in% names(dataset)) {
  #   stop(paste("Column", segmentation_col, "not found in dataset. Check segmentation input!"))
  # }
  
  # Compute P(Z | S)
  if (with_prior) {
    p_z_given_s <- get(paste0("p_z_given_s_", segmentation, "_with_prior"))
  } else {
    p_z_given_s <- get(paste0("p_z_given_s_", segmentation, "_without_prior"))
  }
  
  # dataset <- dataset %>%
  #   filter(total_exposures > 0)
  # 
  # segment_levels <- unique(dataset[[segmentation_col]])  # Get unique segment levels
  # 
  # Compute segment exposures (grouped total exposures)
  # segment_exposures <- compute_segment_sizes_poisson(dataset, segmentation)
  segment_exposures <- segments_poisson
  
  
  # Compute mean exposures for better initial values
  exposure_means <- dataset[[2]] %>%
    group_by(!!sym(segmentation_col)) %>%
    summarise(mean_exposure = mean(total_exposures, na.rm = TRUE), .groups = "drop") %>%
    pull(mean_exposure)
  
  # Prevent log(0) errors
  initial_beta <- log(pmax(exposure_means, 0.01))  
  initial_par <- initial_beta
  
  print(initial_beta)
  

  # Optimize Poisson log-likelihood function
  best_result <- optim(
    par = initial_beta,
    fn = loglikelihood_segments_based_poisson,
    p_z_given_s = p_z_given_s,
    segment_exposures = segment_exposures,
    segmentation = segmentation,
    method = "L-BFGS-B",
    lower = c(-3),     
    upper = c(10)
  )
  
  beta_final <- best_result$par  
  
  if (print_result) {
    cat("\n✅ Optimization converged successfully!\n")
    cat("\nOptimized Poisson Model Parameters:")
    cat("\nBeta:", beta_final, "\n")
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

# ----------------------------- Run the Poisson Model -----------------------------

segmentation <- "age"  # Choose "gender", "age", or "demo"

# ----------------------- For simulated data ----------------------------

n <- 10000000
dataset_frequency <- simulate_exposure_dataset(n, c(0.516, 0.484), c(3.2, 3.0), matrix(c(0.7, 0.3, 0.2, 0.8), ncol = 2),
                                               segmentation = "gender", dispersion = 1, distribution = "poisson",
                                               seed_number = 0, print_simulation = TRUE)
dataset_frequency <- simulate_exposure_dataset(
  n_test = 1e7, 
  fraction_per_segment = c(0.25, 0.25, 0.25, 0.25), 
  exposure_per_segment = c(6, 3.0, 3.5, 4.0), 
  estimation_correctness = matrix(
    c(0.6, 0.2, 0.1, 0.1,
      0.2, 0.6, 0.1, 0.1,
      0.1, 0.1, 0.6, 0.2,
      0.1, 0.1, 0.2, 0.6), 
    ncol = 4, byrow = TRUE),
  segmentation = "age",
  dispersion = 1,
  distribution = "poisson",
  seed_number = 0,
  print_simulation = TRUE
)

dataset_frequency_filtered <- dataset_frequency %>%
  filter(total_exposures > 0)


pois_results <- optimize_loglikelihood_poisson(dataset_frequency_filtered, segmentation, with_prior = FALSE)
eval_results <- evaluation_poisson(dataset_frequency, pois_results$beta, segmentation)
print(eval_results)

# ----------------------- For real data ----------------------------

real_dataset <- read_exposures()
real_dataset_filtered <- real_dataset

pois_results <- optimize_loglikelihood_poisson(real_dataset_filtered, segmentation, with_prior = TRUE)
eval_results <- evaluation_poisson(real_dataset_filtered, pois_results$beta, segmentation)
print(eval_results)