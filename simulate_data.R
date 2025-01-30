simulate_dataset <- function(n_test, fraction_per_segment, exposure_per_segment, estimation_correctness,
                             segmentation, seed_number = 0) {

  set.seed(seed_number)
  segments_count <- length(fraction_per_segment)
  
  if (sum(fraction_per_segment) != 1 || sum(estimation_correctness) != segments_count) {
    print(sum(estimation_correctness))
    stop("Fractions of segments do not add up to 1")
  }
  
  
  segments <- seq(segments_count)
  responses <- c()
  estimate_segments <- c()
  true_segments <- c()
  true_sampled_reach <- c()
  
  for (i in 1:segments_count) {
    response_segment <- rbinom(n_test*fraction_per_segment[i], 1, exposure_per_segment[i])
    responses <- c(responses, response_segment)
    estimate_segments <- c(estimate_segments, sample(seq(segments), n_test*fraction_per_segment[i], 
                                                     replace = TRUE, prob = estimation_correctness[i,]))
    true_segments <- c(true_segments, rep(i,n_test*fraction_per_segment[i]))
    true_sampled_reach <- c(true_sampled_reach, sum(response_segment)/(n_test*fraction_per_segment[i]))
    print(true_sampled_reach)
  }
  
  sample_df <- data.frame(estimate_segments, true_segments, responses)
  if (segmentation == 'gender') {
    colnames(sample_df) <- c("estimated_gender", "true_gender", "response")
  } else if (segmentation == 'age') {
    colnames(sample_df) <- c("estimated_age", "true_age", "response")
  } else if (segmentation == 'demo') {
    colnames(sample_df) <- c("estimated_demo", "true_demo", "response")
  } else {
    close('Segmentation does not exist')
  }
  
  cat(paste('Dataset for male and female devision has been generated, input parameters are: \n sample size: ', n_test))
  for (i in 1:segments_count) {
    cat('\nSegment: ', i,'\n fractions of segment:', fraction_per_segment[i], '\n exposure of segment:', exposure_per_segment[i], 
              '\n correct estimation segments: ', paste(estimation_correctness[i,]), '\n true sampled reach is: ', true_sampled_reach[i])
  }

  
  return(sample_df)
}

