# ----------------------------- Evaluation ------------------------------
evaluate_model <- function(panel_demos, aggregated_exposure, beta) {
  # Compute predicted probabilities
  pred_probabilities <- exp(beta) / (1 + exp(beta))
  
  # Create unique age segments and map predictions
  unique_age_segments <- panel_demos %>%
    group_by(panel_true_age) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(
      true_age_binary = panel_true_age # Binary encoding for age
    )
  
  unique_age_segments <- unique_age_segments %>%
    rowwise() %>%
    mutate(
      prediction = pred_probabilities[true_age_binary]
    ) %>%
    ungroup()
  
  # Compute true fraction response
  true_fraction_response <- panel_demos %>%
    left_join(aggregated_exposure, by = "person_id") %>%
    mutate(
      total_exposures = coalesce(total_exposures, 0),
      response = ifelse(total_exposures > 0, 1, 0) # Binary response
    ) %>%
    group_by(panel_true_age) %>%
    summarise(
      true_fraction_response = mean(response, na.rm = TRUE), # Fraction of response per age
      .groups = 'drop'
    ) %>%
    mutate(true_age_binary = panel_true_age)
  
  # Merge predictions and true response data
  evaluation_data <- unique_age_segments %>%
    left_join(true_fraction_response, by = "true_age_binary")
  
  # Calculate evaluation metrics
  evaluation <- evaluation_data %>%
    summarise(
      mean_absolute_error = mean(abs(prediction - true_fraction_response), na.rm = TRUE),
      mean_squared_error = mean((prediction - true_fraction_response)^2, na.rm = TRUE)
    )
  
  # Print evaluation data and metrics
  print(evaluation_data)
  print(evaluation)
  
  # ----------------------------- True Fraction of Response for Estimated Age ------------------------------
  true_fraction_response_estimated_age <- merged_exposure_data %>%
    group_by(graph_age) %>%
    summarise(
      true_fraction_response = mean(response, na.rm = TRUE), # Fraction of response for estimated age
      .groups = 'drop'
    )
  
  # Print results
  print("True Fraction of Response for Estimated Age:")
  print(true_fraction_response_estimated_age)
  
  print("Evaluation Metrics:")
  print(evaluation)
  
  # Return evaluation metrics
  return(list(evaluation_data = evaluation_data, evaluation = evaluation))
}