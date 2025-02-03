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
  
  return(-log_value_true - 0*log_value_est )
}