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

create_exposure_plot <- function(dataset, segmentation) {
  segment_sizes <- compute_segment_sizes_frequency(dataset, segmentation, simulation = TRUE)
  segment_levels <- unique(segment_sizes[[1]][[paste0("estimated_", segmentation)]])
  
  for (i in segment_levels) {
    print(i)
    segment_values <- segment_sizes[[1]] %>% filter(!!sym(paste0("estimated_", segmentation)) == i)
    print(segment_values)
    plot1 <- ggplot(segment_values, aes(x=total_exposures, y = num_individuals)) +
      geom_point()+
      theme_minimal() +
      xlab("Frequency") +
      ylab("Individual count") +
      labs(title = (paste("Amount of individuals for exposure frequency in ", segmentation, " segment ", i)))
    print(plot1)
  }
}

# Function to create bar plots
plot_bar_chart <- function(df, title, x_label) {
  colnames(df) <- c("Category", "Reach")
  ggplot(df, aes(x = Category, y = Reach)) +
    geom_bar(stat = "identity", position = "dodge", fill = 'indianred2') +
    labs(title = title, x = x_label, y = "Frequency") +
    theme_minimal()
}

# Gender tables
gender_simulation_poiss <- matrix(c(2.156, 1.717, 2.016, 1.876), nrow = 2, byrow = TRUE,
                            dimnames = list(c("Simulated", "Poisson estimates"), c("Male", "Female")))

gender_simulation_neg_binom <- matrix(c(2.156, 1.717, 1.991, 1.962), nrow = 2, byrow = TRUE,
                                  dimnames = list(c("Simulated", "Poisson estimates"), c("Male", "Female")))

age_simulation_poiss <- matrix(c(3.157, 1.858, 1.716, 2.156, 2.772, 1.979, 1.877, 2.115), nrow = 2, byrow = TRUE,
                                      dimnames = list(c("Simulated", "Poisson estimates"), c("lt35", "gt35-lt50", "gt50-lt65", "gt65")))

age_simulation_neg_binom <- matrix(c(3.157, 1.858, 1.716, 2.156, 2.694, 2.049, 2.005, 2.134), nrow = 2, byrow = TRUE,
                               dimnames = list(c("Simulated", "Poisson estimates"), c("lt35", "gt35-lt50", "gt50-lt65", "gt65")))

plot_bar_chart_sim <- function(data, title, x_label) {
  df <- melt(data)
  colnames(df) <- c("Group", "Category", "Percentage")
  ggplot(df, aes(x = Category, y = Percentage, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = title, x = x_label, y = "Frequency") +
    theme_minimal() +
    scale_fill_manual(values = c("Simulated" = "dodgerblue1", "Poisson estimates" = "indianred2"))
}


# Plot for nlsn389882
p1_gender_pois <- plot_bar_chart_sim(gender_simulation_poiss, "Gender exposure frequency (poisson)", "Gender")
p1_gender_neg_binom <- plot_bar_chart_sim(gender_simulation_neg_binom, "Gender exposure frequency (negative binomial)", "Gender")

p2_age_pois <- plot_bar_chart_sim(age_simulation_poiss, "Age exposure frequency (poisson)", "Age")
p2_age_neg_binom <- plot_bar_chart_sim(age_simulation_neg_binom, "Age exposure frequency (negative binomial)", "Age")

plot(p1_gender_pois)
plot(p1_gender_neg_binom)
plot(p2_age_pois)
plot(p2_age_neg_binom)

