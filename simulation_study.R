simulate_evaluate_fraction_2_segments <- function() {
  fraction_gender_1_log <- c(seq(-10, -2, 0.5), seq(-1.5, -0.1, 0.1))
  fraction_gender_1 <- exp(fraction_gender_1_log)
  exposure <- c(0.03, 0.05)
  estimation_correctness <- t(matrix(c(0.75, 0.25, 0.3, 0.7), ncol = 2))
  n_test <- 10**6

  final_results <- data.frame(segment1 = numeric(), segment2 = numeric())
  for (r in 1:10) {
    results <- data.frame(segment1 = numeric(), segment2 = numeric())
    for (i in 1:length(fraction_gender_1)) {
      dataset <- simulate_dataset(n_test = n_test, fraction_per_segment = c(fraction_gender_1[i], 1 - fraction_gender_1[i]), 
                                  exposure_per_segment = exposure, estimation_correctness = estimation_correctness,
                                  segmentation = "gender", seed_number = i, print_simulation = FALSE)
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
    geom_line(aes(y = segment_2, colour = 'Segment 2'), linetype="twodash") +
    geom_hline(yintercept = exposure, color = "blue", linetype = "dotted") +
    labs(
      title = "Prediction of segment reach given logarithm of fraction",
      x = "log of fraction of segment 1",
      y = "Segment reach"
    )
  print(plot1)
}
