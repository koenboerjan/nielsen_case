loglikelihood_segments_based_poisson <- function(beta, p_z_given_s, segment_responses) {
  log_value <- 0
  segment_count <- dim(p_z_given_s)[1]
  
  for (i in 1:segment_count) {
    for (j in 1:length(segment_responses)) {
      exposure_j <- segment_responses[j]  # Observed exposures
      
      value_segment <- 0
      for (s in 1:segment_count) {
        lambda_s <- exp(beta[s])  # Expected exposures for segment `s`
        p_pois <- dpois(exposure_j, lambda_s, log = TRUE)  # Poisson probability
        value_segment <- value_segment + p_z_given_s[i, s] * exp(p_pois)  # Weighted likelihood
      }
      
      log_value <- log_value + log(value_segment + 1e-10)  # Avoid log(0)
    }
  }
  
  # 🔹 Regularization Term (Prevents Extreme Beta Values)
  reg_penalty <- sum(beta^2) * 0.01  # Small penalty (Ridge effect)
  
  return(-1 * (log_value - reg_penalty))  # Minimize negative log-likelihood
}

# ----------------------------- Optimize Poisson Log-Likelihood -----------------------------
optimize_loglikelihood_poisson <- function(dataset, segmentation, with_prior = TRUE, print_result = TRUE) {
  segmentation_col <- paste0("estimated_", segmentation)  
  
  if (!segmentation_col %in% names(dataset)) {
    stop(paste("Column", segmentation_col, "not found in dataset. Check segmentation input!"))
  }
  
  dataset <- dataset %>%
    filter(response == 1, total_exposures <= 10)  
  
  # Compute P(Z | S)
  if (with_prior) {
    p_z_given_s <- compute_p_z_given_s_including_prior(dataset, segmentation)
  } else {
    p_z_given_s <- compute_p_z_given_s(dataset, segmentation)
  }
  
  segment_count <- dim(p_z_given_s)[1]
  segment_responses <- dataset$total_exposures
  
  # 🔹 Improved Initial Beta using Log Mean Exposures
  exposure_means <- dataset %>%
    group_by(!!sym(segmentation_col)) %>%
    summarise(mean_exposure = mean(total_exposures, na.rm = TRUE), .groups = "drop") %>%
    pull(mean_exposure)
  
  initial_beta <- log(exposure_means + 1e-3)  # Prevent log(0)
  initial_par <- initial_beta  
  
  # 🔹 Prevent Optimization From Going Wild (Bounds)
  best_result <- optim(
    par = initial_par,
    fn = function(params) {
      loglikelihood_segments_based_poisson(params, p_z_given_s, segment_responses)
    },
    method = "L-BFGS-B",
    lower = rep(log(0.1), segment_count),  # Min lambda = 0.1
    upper = rep(log(10), segment_count),   # Max lambda = 10
    control = list(fnscale = -1, maxit = 1000, factr = 1e7)
  )
  
  beta_final <- best_result$par  
  
  if (print_result) {
    cat("\n✅ Optimization converged successfully!\n")
    cat("\nOptimized Poisson Model Parameters:")
    cat("\nBeta:", beta_final, "\n")
  }
  
  return(list(beta = beta_final, optimization = best_result))
}

# ----------------------------- Compute Expected Frequency Given Reach -----------------------------
evaluation_poisson <- function(dataset, beta, segmentation) {
  segmentation_col <- paste0("estimated_", segmentation)
  
  pred_expected_exposures <- exp(beta)
  
  exposure_for_estimated_seg <- dataset %>%
    filter(response == 1) %>%
    group_by(!!sym(segmentation_col)) %>%
    summarise(
      mean_exposures = mean(total_exposures, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(predicted_exposures = pred_expected_exposures)
  
  return(exposure_for_estimated_seg)
}

# ----------------------------- Run the Poisson Model -----------------------------
segmentation <- "gender"  # Choose "gender", "age", or "demo"
pois_results <- optimize_loglikelihood_poisson(real_dataset, segmentation, with_prior = TRUE)
eval_results <- evaluation_poisson(real_dataset, pois_results$beta, segmentation)

# Print results
print(eval_results)
