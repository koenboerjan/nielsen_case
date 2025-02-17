simulate_exposure_dataset <- function(n_test, fraction_per_segment, exposure_per_segment, estimation_correctness,
                                      segmentation, dispersion = 1, distribution = "poisson",
                                      seed_number = 0, print_simulation = TRUE) {
  set.seed(seed_number)
  segments_count <- length(fraction_per_segment)
  if (sum(fraction_per_segment) != 1 || sum(rowSums(estimation_correctness)) != segments_count) {
    stop("Fractions of segments do not add up to 1 or estimation matrix is incorrect.")
  }
  segments <- seq(segments_count)
  total_exposures <- c()
  estimate_segments <- c()
  true_segments <- c()
  true_sampled_mean_exposure <- c()
  for (i in 1:segments_count) {
    segment_size <- round(n_test * fraction_per_segment[i])
    # 🔹 Simulate exposures using Poisson or Negative Binomial
    if (distribution == "poisson") {
      exposure_segment <- rpois(segment_size, lambda = exposure_per_segment[i])
    } else if (distribution == "negative_binomial") {
      exposure_segment <- rnbinom(segment_size, size = dispersion, mu = exposure_per_segment[i])
    } else {
      stop("Invalid distribution. Choose 'poisson' or 'negative_binomial'.")
    }
    total_exposures <- c(total_exposures, exposure_segment)
    estimate_segments <- c(estimate_segments, sample(seq(segments), segment_size, 
                                                     replace = TRUE, prob = estimation_correctness[i, ]))
    true_segments <- c(true_segments, rep(i, segment_size))
    true_sampled_mean_exposure <- c(true_sampled_mean_exposure, mean(exposure_segment))
  }
  sample_df <- data.frame(estimate_segments, true_segments, total_exposures)
  # 🔹 Rename columns based on segmentation type
  if (segmentation == 'gender') {
    colnames(sample_df) <- c("estimated_gender", "true_gender", "total_exposures")
  } else if (segmentation == 'age') {
    colnames(sample_df) <- c("estimated_age", "true_age", "total_exposures")
  } else if (segmentation == 'demo') {
    colnames(sample_df) <- c("estimated_demo", "true_demo", "total_exposures")
  } else {
    stop('Segmentation does not exist')
  }
  # 🔹 Print simulation details
  if (print_simulation) {
    cat(paste0('\nDataset for ', segmentation, ' segmentation has been generated. Sample size: ', n_test, "\n"))
    for (i in 1:segments_count) {
      cat('\nSegment:', i,
          '\n Fraction of segment:', fraction_per_segment[i], 
          '\n Mean exposure of segment:', exposure_per_segment[i], 
          '\n Correct estimation segments:', paste(estimation_correctness[i,]), 
          '\n True sampled mean exposure:', true_sampled_mean_exposure[i])
    }
  }
  return(sample_df)
}
