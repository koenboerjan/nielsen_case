
#--------------------------------Compute True Segment Sizes-------------------------------------
# compute_true_segment_sizes<-function(dataset,segmentation){
#   #Delete the observations with NA values
#   dataset <- dataset %>%
#     filter(!is.na(true_gender), !is.na(true_age), !is.na(true_demo))
#   print(head(dataset))
#   # Summarize response and non-response counts
#   if (segmentation == "gender") {
#     segments_response <- dataset %>%
#       group_by(true_gender) %>%
#       summarise(
#         response_count = sum(response),
#         no_response_count = sum(1 - response),
#         .groups = 'drop'
#       )
#   } else if (segmentation == "age") {
#     segments_response <- dataset %>%
#       group_by(true_age) %>%
#       summarise(
#         response_count = sum(response),
#         no_response_count = sum(1 - response),
#         .groups = 'drop'
#       )
#   } else if (segmentation == "demo") {
#     segments_response <- dataset %>%
#       group_by(true_demo) %>%
#       summarise(
#         response_count = sum(response),
#         no_response_count = sum(1 - response),
#         .groups = 'drop'
#       )
#   }
# 
#   mat_segments_response <- matrix(
#     c(segments_response$response_count, segments_response$no_response_count),
#     ncol = 2
#   )
# 
#   print(segments_response)
# 
#   return(mat_segments_response)
# }

# ----------------------------- Compute Log Likelihood for Segments -----------------------------
loglikelihood_segments_based <- function(beta, p_z_given_s, segment_responses, true_segments, true_weight, est_weight) {
  log_value_est <- 0
  log_value_true <- 0
  segment_count <- dim(p_z_given_s)[1]
  
  # Compute likelihood for estimated demographics
  for (i in 1:segment_count) {
    for (j in 1:2) {
      response_j <- ifelse(j == 1, 1, 0)
      value_segment <- 0
      
      for (s in 1:segment_count) {
        e_z_beta <- plogis(beta[s])
        value_segment <- value_segment + p_z_given_s[i, s] * (e_z_beta^response_j) * ((1 - e_z_beta)^(1 - response_j))
        
      }
      
      log_value_est <- log_value_est + segment_responses[i, j] * log(max(value_segment, 1e-10))
    }
  }
  
  # Compute likelihood for true demographics
  for (z in 1:segment_count) {
    for (j in 1:2) {
      response_j <- ifelse(j == 1, 1, 0)
      e_z_beta <- plogis(beta[z])  # Convert beta to probability
      
      # Compute p(y_j | z)
      p_y_given_z <- (e_z_beta^response_j) * ((1 - e_z_beta)^(1 - response_j))
      
      # Log likelihood contribution for this (y_j, z) pair
      log_value_true <- log_value_true + true_segments[z, j] * log(max(p_y_given_z, 1e-10))
      
      # Debugging output
      # print(paste("DEBUG: log_value_true at step (z =", z, ", j =", j, "):", 
      #             "true_segments =", true_segments[z, j], 
      #             "p_y_given_z =", p_y_given_z, 
      #             "log_value_true =", log_value_true))
    }
  }
  # Debugging print
  # print(paste("true_weight:", true_weight, "est_weight:", est_weight, 
  #             "log_value_true:", round(log_value_true, 4), 
  #             "log_value_est:", round(log_value_est, 4)))
  
  return(-true_weight * log_value_true - est_weight * log_value_est)
  
  
}
#----------------------------------------------------------------------------------------------------------------------------------------#
#Load libraries and functions
prepare_setup <- function() {
  source("requirements.R")
  install_requirements()
}


prepare_setup()


# Load data
real_dataset<-read_exposures()

# Compute conditional probabilities for gender segmentation
conditional_probs <- compute_p_z_given_s_including_prior(dataset = real_dataset, segmentation = 'gender')
print(conditional_probs)
results <- compute_segment_sizes(dataset = real_dataset, segmentation = 'gender',use_true_seperate = TRUE)
print(results)
# Extract segment responses and true segments from the results
segment_responses <- results[[1]]
print(segment_responses)
true_segments <- results[[2]]
print(true_segments)
segment_count <- dim(conditional_probs)[1]

# Define empty dataframe to store results
results_df <- data.frame(true_weight = numeric(),
                         est_weight = numeric(),
                         pr_1_given_female = numeric(),
                         pr_1_given_male = numeric())

# Compute true response fractions
true_fractions <- real_dataset %>% 
  filter(!is.na(true_gender), !is.na(true_age), !is.na(true_demo)) %>%
  group_by(true_gender, response) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(true_gender) %>%
  mutate(fraction = count / sum(count)) %>%
  filter(response == 1) %>%  # Keep only response=1
  select(true_gender, fraction)

true_female_response <- true_fractions$fraction[true_fractions$true_gender == 0]
true_male_response <- true_fractions$fraction[true_fractions$true_gender == 1]
print(true_fractions)
print(true_male_response)
print(true_female_response)

# Grid search over different weight combinations
for (true_weight in seq(1, 1000, by = 5)) {
  for (est_weight in seq(1, 20, by = 5)) {
    
    # ⚠️ Randomize initial betas for each iteration to avoid local minimal
    initial_betas <- rnorm(segment_count, mean = 0, sd = 1)
    
    
    maximization <- optim(
      par = initial_betas,
      fn = loglikelihood_segments_based,
      p_z_given_s = conditional_probs,
      segment_responses = segment_responses,
      true_segments = true_segments,
      true_weight = true_weight,
      est_weight = est_weight,
      method = "BFGS",
      control = list(maxit = 1000, reltol = 1e-6)
    )
    
    optimal_betas <- maximization$par
    print(paste("true_weight:", true_weight, "est_weight:", est_weight, 
                "Beta:", paste(round(optimal_betas, 4), collapse = ", ")))
    
    pr_1_given_female <- plogis(optimal_betas[2])
    pr_1_given_male <- plogis(optimal_betas[1])
    
    # Clip probabilities to avoid extreme values
    pr_1_given_female <- pmin(pr_1_given_female, 0.15)
    pr_1_given_male <- pmin(pr_1_given_male, 0.15)
    
    # Store results correctly in results_df
    results_df <- rbind(results_df, data.frame(true_weight, est_weight, pr_1_given_female, pr_1_given_male))
  }
}

# Calculate Euclidean distance for each point
results_df <- results_df %>%
  mutate(distance = sqrt((pr_1_given_male - true_male_response)^2 + 
                           (pr_1_given_female - true_female_response)^2))

# Find the row with the minimum distance
closest_point <- results_df[which.min(results_df$distance), ]
print(closest_point)
# Plot with the closest point highlighted
ggplot(results_df, aes(x = pr_1_given_male, y = pr_1_given_female)) +
  geom_point(aes(color = distance), alpha = 0.5) +  # Color by distance (closer points are darker)
  scale_color_gradient(low = "red", high = "blue") +  # Red for closer points
  # True response values (red triangle)
  annotate("point", x = true_male_response, y = true_female_response, 
           color = "green", size = 4, shape = 17) +  
  # Closest point to the true fraction (highlighted in blue)
  annotate("point", x = closest_point$pr_1_given_male, y = closest_point$pr_1_given_female, 
           color = "red", size = 4, shape = 18) +
  labs(title = "Effect of Different Weights on Estimated Probabilities",
       x = "Estimated Pr(Response | Male)",
       y = "Estimated Pr(Response | Female)",
       color = "Distance to True Fraction") +
  theme_minimal()


#Do the same for age
conditional_probs_age <- compute_p_z_given_s_including_prior(dataset = real_dataset, segmentation = 'age')
segment_responses_age <- compute_segment_sizes(dataset = real_dataset, segmentation = 'age')
true_segments_age <- compute_true_segment_sizes(real_dataset, 'age')
segment_count <- dim(conditional_probs_age)[1]

# Define empty dataframe to store results
results_age_df <- data.frame(true_weight = numeric(),
                         est_weight = numeric(),
                         pr_1_given_lt35 = numeric(),
                         pr_1_given_gt35_lt50 = numeric(),
                         pr_1_given_gt50_lt65=numeric(),
                         pr_1_given_gt65=numeric())

# Grid search over different weight combinations
for (true_weight in seq(1, 1000, by = 5)) {
  for (est_weight in seq(1, 20, by = 5)) {
    
    # ⚠️ Randomize initial betas for each iteration to avoid local minimal
    initial_betas <- rnorm(segment_count, mean = 0, sd = 1)
    
    maximization <- optim(
      par = initial_betas,
      fn = loglikelihood_segments_based,
      p_z_given_s = conditional_probs_age,
      segment_responses = segment_responses_age,
      true_segments = true_segments_age,
      true_weight = true_weight,
      est_weight = est_weight,
      control = list(maxit = 1000, reltol = 1e-6),
      method = 'BFGS'
    )
    
    optimal_betas <- maximization$par
    print(paste("true_weight:", true_weight, "est_weight:", est_weight, 
                "Beta:", paste(round(optimal_betas, 4), collapse = ", ")))
    
    pr_1_given_lt35 <- plogis(optimal_betas[1])
    pr_1_given_gt35_lt50 <- plogis(optimal_betas[2])
    pr_1_given_gt50_lt65 <- plogis(optimal_betas[3])
    pr_1_given_gt65<-plogis(optimal_betas[4])
    
    
    # Store results correctly in results_df
    results_age_df <- rbind(results_age_df, data.frame(true_weight, est_weight, pr_1_given_lt35, pr_1_given_gt35_lt50,pr_1_given_gt50_lt65,
                                               pr_1_given_gt65))
  }
}
results_age_df <- results_age_df %>%
  rowwise() %>%
  mutate(distance = sqrt((pr_1_given_lt35 - true_fractions_age$fraction[1])^2 +
                           (pr_1_given_gt35_lt50 - true_fractions_age$fraction[2])^2 +
                           (pr_1_given_gt50_lt65 - true_fractions_age$fraction[3])^2 +
                           (pr_1_given_gt65 - true_fractions_age$fraction[4])^2)) %>%
  ungroup()

# Find the row with the minimum distance
closest_point <- results_age_df[which.min(results_age_df$distance), ]
print(closest_point)

# Compute true response fractions
true_fractions_age <- real_dataset %>% 
  filter(!is.na(true_gender), !is.na(true_age), !is.na(true_demo)) %>%
  group_by(true_age, response) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(true_age) %>%
  mutate(fraction = count / sum(count)) %>%
  filter(response == 1) %>%  # Keep only response=1
  select(true_age, fraction)

print(true_fractions_age)

# Compute conditional probabilities for 'demo'
conditional_probs_demo <- compute_p_z_given_s_including_prior(dataset = real_dataset, segmentation = 'demo')
segment_responses_demo <- compute_segment_sizes(dataset = real_dataset, segmentation = 'demo')
true_segments_demo <- compute_true_segment_sizes(real_dataset, 'demo')
segment_count_demo <- dim(conditional_probs_demo)[1]  # Should be 5 now

# Define empty dataframe to store results
results_demo_df <- data.frame(
  true_weight = numeric(),
  est_weight = numeric(),
  pr_1_given_demo1 = numeric(),
  pr_1_given_demo2 = numeric(),
  pr_1_given_demo3 = numeric(),
  pr_1_given_demo4 = numeric(),
  pr_1_given_demo5 = numeric()
)

# Grid search over different weight combinations
for (true_weight in seq(1, 1000, by = 5)) {
  for (est_weight in seq(1, 20, by = 5)) {
    
    # ⚠️ Randomize initial betas for each iteration to avoid local minima
    initial_betas_demo <- rnorm(segment_count_demo, mean = 0, sd = 1)
    
    maximization_demo <- optim(
      par = initial_betas_demo,
      fn = loglikelihood_segments_based,
      p_z_given_s = conditional_probs_demo,
      segment_responses = segment_responses_demo,
      true_segments = true_segments_demo,
      true_weight = true_weight,
      est_weight = est_weight,
      control = list(maxit = 1000, reltol = 1e-6),
      method = 'BFGS'
    )
    
    optimal_betas_demo <- maximization_demo$par
    print(paste("true_weight:", true_weight, "est_weight:", est_weight, 
                "Beta:", paste(round(optimal_betas_demo, 4), collapse = ", ")))
    
    pr_1_given_demo1 <- plogis(optimal_betas_demo[1])
    pr_1_given_demo2 <- plogis(optimal_betas_demo[2])
    pr_1_given_demo3 <- plogis(optimal_betas_demo[3])
    pr_1_given_demo4 <- plogis(optimal_betas_demo[4])
    pr_1_given_demo5 <- plogis(optimal_betas_demo[5])
    
    # Store results correctly in results_demo_df
    results_demo_df <- rbind(results_demo_df, data.frame(true_weight, est_weight, pr_1_given_demo1, pr_1_given_demo2, pr_1_given_demo3, pr_1_given_demo4, pr_1_given_demo5))
  }
}


# Compute true response fractions for demo
true_fractions_demo <- real_dataset %>% 
  filter(!is.na(true_gender), !is.na(true_age), !is.na(true_demo)) %>%
  group_by(true_demo, response) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(true_demo) %>%
  mutate(fraction = count / sum(count)) %>%
  filter(response == 1) %>%  # Keep only response=1
  select(true_demo, fraction)

print(true_fractions_demo)

results_demo_df <- results_demo_df %>%
  rowwise() %>%
  mutate(distance = sqrt(
    (pr_1_given_demo1 - true_fractions_demo$fraction[true_fractions_demo$true_demo == 1])^2 +
      (pr_1_given_demo2 - true_fractions_demo$fraction[true_fractions_demo$true_demo == 2])^2 +
      (pr_1_given_demo3 - true_fractions_demo$fraction[true_fractions_demo$true_demo == 3])^2 +
      (pr_1_given_demo4 - true_fractions_demo$fraction[true_fractions_demo$true_demo == 4])^2 +
      (pr_1_given_demo5 - true_fractions_demo$fraction[true_fractions_demo$true_demo == 5])^2
  )) %>%
  ungroup()

# Find the row with the minimum distance
closest_point_demo <- results_demo_df[which.min(results_demo_df$distance),]
print(closest_point_demo)

conditional_probs_full <- compute_p_z_given_s_full(dataset = real_dataset)
print(conditional_probs_full)
universe_estimates<-read_universe_estimates()