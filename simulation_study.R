simulate_evaluate_method <- function() {
  n <- 10**7
  dataset <- simulate_dataset(n, c(0.516, 0.484), c(0.02,0.06), matrix(c(0.7,0.3,0.2,0.8), ncol = 2),
                              "gender", seed_number = 0, print_simulation = TRUE) 
  beta <- optimize_loglikelihood(real_dataset, segmentation = 'gender', with_prior = FALSE)
  
  
  
}
 
simulate_evaluate_fraction_2_segments <- function() {
  fraction_gender_1_log <- c(seq(-10, -2, 0.5), seq(-1.5, -0.1, 0.1))
  fraction_gender_1 <- exp(fraction_gender_1_log)
  exposure <- c(0.03, 0.05)
  estimation_correctness <- t(matrix(c(0.75, 0.25, 0.3, 0.7), ncol = 2))
  n_test <- 10**6

  final_results <- data.frame(segment1 = numeric(), segment2 = numeric())
  for (r in 1:10) {
    print(r)
    results <- data.frame(segment1 = numeric(), segment2 = numeric())
    for (i in 1:length(fraction_gender_1)) {
      dataset <- simulate_dataset(n_test = n_test, fraction_per_segment = c(fraction_gender_1[i], 1 - fraction_gender_1[i]), 
                                  exposure_per_segment = exposure, estimation_correctness = estimation_correctness,
                                  segmentation = "gender", seed_number = (i*10 + 15*r), print_simulation = FALSE)
      optimal_beta <- optimize_loglikelihood(dataset, segmentation = 'gender', with_prior = FALSE, print_result = FALSE)
      results <- rbind(results, exp(optimal_beta$par) / (1 + exp(optimal_beta$par)))
      
    }
    if (r == 1) {
      final_results <- results
    } else {
      final_results <- final_results + results
    }
  }
  
  final_results <- final_results/10
  final_results$x_axis <- fraction_gender_1_log
  
  colnames(final_results) <- c('segment_1', 'segment_2', 'x_axis')
  
  plot1 <- ggplot(final_results, aes(x=x_axis)) + 
    geom_line(aes(y = segment_1, colour = 'Segment 1')) + 
    geom_line(aes(y = segment_2, colour = 'Segment 2'), linetype = 5) +
    geom_hline(aes(yintercept = exposure[1], color = c("True segment 1")), linetype = 6) +
    geom_hline(aes(yintercept = exposure[2], color = c("True segment 2")), linetype = 6) + 
    labs(
      title = "Prediction of segment reach given logarithm of fraction",
      x = "log of fraction of segment 1",
      y = "Segment reach"
    )
  print(plot1)
}

simulate_evaluate_fraction_3_segments <- function() {
  fraction_age_1_log <- c(seq(-6, -2, 1), seq(-1.5, -0.1, 0.3))
  fraction_age_2_log <- c(seq(-6, -2, 1), seq(-1.5, -0.1, 0.3))
  fraction_age_1 <- round(exp(fraction_age_1_log),6)
  fraction_age_2 <- round(exp(fraction_age_2_log),6)
  fractions <- data.frame(segment1 = numeric(), segment2 = numeric(), segment3 = numeric())
  for (i in 1:length(fraction_age_1)) {
    for (j in 1:length(fraction_age_1)) {
      fraction_age_3 <- 1 - (fraction_age_1[i] + fraction_age_2[j])
      if (fraction_age_3 > 0) {
        fractions <- rbind(fractions, c(fraction_age_1[i], fraction_age_2[j], fraction_age_3))
      }
    }
  }
  exposure <- c(0.03, 0.05, 0.07)
  estimation_correctness <- t(matrix(c(0.75, 0.20, 0.05, 0.1, 0.7, 0.2, 0.1, 0.25, 0.65), ncol = 3))
  
  n_test <- 10**6
  fractions <- as.matrix(fractions)
  
  final_results <- data.frame(segment1 = numeric(), segment2 = numeric(), segment3 = numeric())
  for (r in 1:10) {
    print(r)
    results <- data.frame(segment1 = numeric(), segment2 = numeric(), segment3 = numeric())
    for (i in 1:dim(fractions)[1]) {
      print(i)
      dataset <- simulate_dataset(n_test = n_test, fraction_per_segment = fractions[i,], 
                                  exposure_per_segment = exposure, estimation_correctness = estimation_correctness,
                                  segmentation = "age", seed_number = (i*10+15*r), print_simulation = FALSE)
      optimal_beta <- optimize_loglikelihood(dataset, segmentation = 'age', with_prior = FALSE, print_result = FALSE, simulation = TRUE)
      results <- rbind(results, exp(optimal_beta$par) / (1 + exp(optimal_beta$par)))
      
    }
    if (r == 1) {
      final_results <- results
    } else {
      final_results <- final_results + results
    }
  }
  
  final_results <- final_results/10
  final_results$age_1_axis <- fractions[,1]
  final_results$age_2_axis <- fractions[,2]
  
  colnames(final_results) <- c('segment_1', 'segment_2', 'segment_3', 'segment_1_fraction', 'segment_2_fraction')
  
  x_plane <- c(0, 0, 1, 1)
  y_plane <- c(0, 1, 0, 1)
  z_plane_1 <- rep(0.03, 4)
  z_plane_2 <- rep(0.05, 4)
  z_plane_3 <- rep(0.07, 4)
  
  # Create plot1
  plot1 <- plot_ly() %>%
    # Add 3D Scatter Points
    add_trace(
      data = final_results,
      x = ~segment_1_fraction, y = ~segment_2_fraction, z = ~segment_1,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 5),
      name = "Segment 1"
    ) %>%
    # Add a Transparent Plane
    add_trace(
      x = x_plane, y = y_plane, z = z_plane_1,
      type = "mesh3d",
      opacity = 0.5,
      color = "red",
      showscale = FALSE
    ) %>%
    layout(title = "Reach of segment 1 given segment fractions",
           scene = list(xaxis = list(title = "Fraction 1"),
                        yaxis = list(title = "Fraction 2"),
                        zaxis = list(title = "Reach segment 1", range = c(0, 0.1))))
  print(plot1)
  
  # Create plot2
  plot2 <- plot_ly() %>%
    # Add 3D Scatter Points
    add_trace(
      data = final_results,
      x = ~segment_1_fraction, y = ~segment_2_fraction, z = ~segment_2,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 5),
      name = "Segment 2"
    ) %>%
    # Add a Transparent Plane
    add_trace(
      x = x_plane, y = y_plane, z = z_plane_2,
      type = "mesh3d",
      opacity = 0.5,
      color = "red",
      showscale = FALSE
    ) %>%
    layout(title = "Reach of segment 2 given segment fractions",
           scene = list(xaxis = list(title = "Fraction 1"),
                        yaxis = list(title = "Fraction 2"),
                        zaxis = list(title = "Reach segment 2", range = c(0, 0.1))))
  print(plot2)
  
  # Create plot3
  plot3 <- plot_ly() %>%
    # Add 3D Scatter Points
    add_trace(
      data = final_results,
      x = ~segment_1_fraction, y = ~segment_2_fraction, z = ~segment_3,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 5),
      name = "Segment 3"
    ) %>%
    # Add a Transparent Plane
    add_trace(
      x = x_plane, y = y_plane, z = z_plane_3,
      type = "mesh3d",
      opacity = 0.5,
      color = "red",
      showscale = FALSE
    ) %>%
    layout(title = "Reach of segment 3 given segment fractions",
           scene = list(xaxis = list(title = "Fraction 1"),
                        yaxis = list(title = "Fraction 2"),
                        zaxis = list(title = "Reach segment 3", range = c(0, 0.1))))
  print(plot3)
}

### Simulation study exposure
simulate_evaluate_exposure_2_segments <- function() {
  fraction_gender <- c(0.516, 0.484)
  exposure_1_log <- c(seq(-10, -3, 1),seq(-2.5, -1, 0.5))
  exposure_2_log <- c(seq(-10, -3, 1),seq(-2.5, -1, 0.5))
  exposure_1 <- exp(exposure_1_log)
  exposure_2 <- exp(exposure_2_log)
  exposure <- data.frame(segment_1 = numeric(), segment_2 = numeric())

  for (i in 1:length(exposure_1)) {
    for (j in 1:length(exposure_2)) {
        exposure <- rbind(exposure, c(exposure_1[i], exposure_2[j]))
      }
    }
  exposure <- as.matrix(exposure)

  estimation_correctness <- t(matrix(c(0.75, 0.25, 0.3, 0.7), ncol = 2))
  n_test <- 10**6
  
  final_results <- data.frame(segment1 = numeric(), segment2 = numeric())
  for (r in 1:10) {
    print(r)
    results <- data.frame(segment1 = numeric(), segment2 = numeric())
    for (i in 1:dim(exposure)[1]) {
      print(i)
      dataset <- simulate_dataset(n_test = n_test, fraction_per_segment = fraction_gender, 
                                  exposure_per_segment = exposure[i,], estimation_correctness = estimation_correctness,
                                  segmentation = "gender", seed_number = (i*10 + 15*r), print_simulation = FALSE)
      optimal_beta <- optimize_loglikelihood(dataset, segmentation = 'gender', with_prior = FALSE, print_result = FALSE, simulation = TRUE)
      results <- rbind(results, exp(optimal_beta$par) / (1 + exp(optimal_beta$par)))
      
    }
    if (r == 1) {
      final_results <- results
    } else {
      final_results <- final_results + results
    }
  }
  
  final_results <- final_results/10
  final_results$segment_1_axis <- exposure[,1]
  final_results$segment_2_axis <- exposure[,2]
  
  colnames(final_results) <- c('segment_1', 'segment_2', 'segment_1_exposure', 'segment_2_exposure')
  
  x_plane <- exposure_1[c(1,1,12,12)]
  y_plane <- exposure_2[c(1, 12, 1, 12)]
  z_plane_1 <- exposure_1[c(1,1,12,12)]
  z_plane_2 <- exposure_2[c(1,12,1,12)]
  
  # Create plot1
  plot1 <- plot_ly() %>%
    # Add 3D Scatter Points
    add_trace(
      data = final_results,
      x = ~segment_1_exposure, y = ~segment_2_exposure, z = ~segment_1,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 5),
      name = "Segment 1"
    ) %>%
    add_trace(
      x = x_plane, y = y_plane, z = z_plane_1,
      type = "mesh3d",
      opacity = 0.5,
      color = "red",
      showscale = FALSE
    ) %>%
    layout(title = "Reach of segment 1 given segment exposures",
           scene = list(xaxis = list(title = "Exposure 1"),
                        yaxis = list(title = "Exposure 2"),
                        zaxis = list(title = "Reach segment 1", range = c(0, 1))))
  print(plot1)
  
  # Create plot2
  plot2 <- plot_ly() %>%
    # Add 3D Scatter Points
    add_trace(
      data = final_results,
      x = ~segment_1_exposure, y = ~segment_2_exposure, z = ~segment_2,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 5),
      name = "Segment 2"
    ) %>%
    add_trace(
      x = x_plane, y = y_plane, z = z_plane_2,
      type = "mesh3d",
      opacity = 0.5,
      color = "red",
      showscale = FALSE
    ) %>%
    layout(title = "Reach of segment 2 given segment exposures",
           scene = list(xaxis = list(title = "Exposure 1"),
                        yaxis = list(title = "Exposure 2"),
                        zaxis = list(title = "Reach segment 2", range = c(0, 1))))
  print(plot2)
  
}


### Simulation estimation accuracy exposure
simulate_evaluate_accuracy_2_segments <- function() {
  fraction_gender <- c(0.516, 0.484)
  exposure <- c(0.03, 0.05)
  accuracy_estimate_1 <- c(seq(0.5, 0.65, 0.025), seq(0.7, 0.95, 0.05))
  accuracy_estimate_2 <- c(seq(0.5, 0.65, 0.025), seq(0.7, 0.95, 0.05))
  accuracy_estimation <- data.frame(segment_1_correct = numeric(), segment_1_wrong = numeric(), 
                                    segment_2_correct = numeric(),segment_2_wrong = numeric())
  
  for (i in 1:length(accuracy_estimate_1)) {
    for (j in 1:length(accuracy_estimate_2)) {
      accuracy_estimation <- rbind(accuracy_estimation, c(accuracy_estimate_1[i], 1 - accuracy_estimate_1[i], 
                                                          1- accuracy_estimate_2[j], accuracy_estimate_2[j]))
    }
  }
  accuracy_estimation <- as.matrix(accuracy_estimation)
   
  n_test <- 10**6
  
  final_results <- data.frame(segment1 = numeric(), segment2 = numeric())
  for (r in 1:10) {
    print(r)
    results <- data.frame(segment1 = numeric(), segment2 = numeric())
    for (i in 1:dim(accuracy_estimation)[1]) {
      print(i)
      estimation_correctness_mat <- as.matrix(t(matrix(accuracy_estimation[i,], ncol = 2)))
      dataset <- simulate_dataset(n_test = n_test, fraction_per_segment = fraction_gender, 
                                  exposure_per_segment = exposure, estimation_correctness = estimation_correctness_mat,
                                  segmentation = "gender", seed_number = (i*10 + 15*r), print_simulation = FALSE, simulation = TRUE)
      optimal_beta <- optimize_loglikelihood(dataset, segmentation = 'gender', with_prior = FALSE, print_result = FALSE)
      results <- rbind(results, exp(optimal_beta$par) / (1 + exp(optimal_beta$par)))
      
    }
    if (r == 1) {
      final_results <- results
    } else {
      final_results <- final_results + results
    }
  }
  
  final_results <- final_results/10
  final_results$correctness_1 <- accuracy_estimation[,1]
  final_results$correctness_2 <- accuracy_estimation[,4]
  
  colnames(final_results) <- c('segment_1', 'segment_2', 'segment_1_accuracy', 'segment_2_accuracy')
  
  x_plane <- c(0.5, 0.5, 1, 1)
  y_plane <- c(0.5, 1, 0.5, 1)
  z_plane_1 <- rep(0.03, 4)
  z_plane_2 <- rep(0.05, 4)
  
  # Create plot1
  plot1 <- plot_ly() %>%
    # Add 3D Scatter Points
    add_trace(
      data = final_results,
      x = ~segment_1_accuracy, y = ~segment_2_accuracy, z = ~segment_1,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 5),
      name = "Segment 1"
    ) %>%
    # Add a Transparent Plane
    add_trace(
      x = x_plane, y = y_plane, z = z_plane_1,
      type = "mesh3d",
      opacity = 0.5,
      color = "red",
      showscale = FALSE
    )  %>%
    layout(title = "Reach of segment 1 given segment estimation accuracy",
           scene = list(xaxis = list(title = "Accuracy 1"),
                        yaxis = list(title = "Accuracy 2"),
                        zaxis = list(title = "Reach segment 1", range = c(0, 0.1))))
  print(plot1)
  
  # Create plot2
  plot2 <- plot_ly() %>%
    # Add 3D Scatter Points
    add_trace(
      data = final_results,
      x = ~segment_1_accuracy, y = ~segment_2_accuracy, z = ~segment_2,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 5),
      name = "Segment 2"
    ) %>%
    # Add a Transparent Plane
    add_trace(
      x = x_plane, y = y_plane, z = z_plane_2,
      type = "mesh3d",
      opacity = 0.5,
      color = "red",
      showscale = FALSE
    ) %>%
    layout(title = "Reach of segment 2 given segment estimation accuracy",
           scene = list(xaxis = list(title = "Accuracy 1"),
                        yaxis = list(title = "Accuracy 2"),
                        zaxis = list(title = "Reach segment 2", range = c(0, 0.1))))
  print(plot2)
  
}



