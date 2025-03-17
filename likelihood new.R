
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
      
    }
  }
  
  return(-true_weight * log_value_true - est_weight * log_value_est)
  
  
}
#----------------------------------------------------------------------------------------------------------------------------------------#
#Load libraries and functions
prepare_setup <- function() {
  source("requirements.R")
  install_requirements()
}


prepare_setup()


#Read and merge the datasets
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
segment_sizes(dataset = real_dataset, segmentation = 'demo')
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

<<<<<<< HEAD
conditional_probs_full <- compute_p_z_given_s_full(dataset = real_dataset)
=======


##########################################################################################
#COMBINED SEGMENTS ESTIMATIONS
#Now estimate the conditional probabilities for the full set of demographic combinations
<<<<<<< HEAD
<<<<<<< HEAD
conditional_probs_full <- compute_p_z_given_s_including_prior(dataset = real_dataset, segmentation = 'full')
>>>>>>> 2a58319eb8ab22955090655882a3b80ce6ccbbe9
print(conditional_probs_full)
=======
# conditional_probs_full <- compute_p_z_given_s_including_prior(dataset = real_dataset, segmentation = 'full')
# print(conditional_probs_full)
>>>>>>> 26cc42abda0b9418c5adf2a80f5f64c0b1340f34
=======
>>>>>>> 4d8661d09744ada2af47323259f20d4733f380d6

#Store the observations with known true demographics separately
segment_responses_full<-real_dataset[[2]]
true_segments_full<-real_dataset[[1]]
#Calculate the 40 x 40 conditional matrix for the combined segments
conditional_probs_full <-p_z_given_s_full()
print(true_segments_full)

#Make sure the demographic segments appear at the same order in both datasets
segment_responses_full <- segment_responses_full %>%
  arrange(estimated_demo, estimated_age,estimated_gender)
true_segments_full<-true_segments_full%>%
  arrange(true_demo, true_age,true_gender)
 
#In the dataset of users with known true demographics, calculate the fraction of exposed individuals
segment_data <- true_segments_full %>%
  as.data.frame() %>%
  mutate(fraction = ifelse((response_count + no_response_count) > 0,
                      response_count / (response_count + no_response_count),
                      0) 
  )
print(segment_data)

#Define the zero vector as the initial beta for the optimization
beta_length <- nrow(conditional_probs_full)  
initial_beta <- rep(0, beta_length) 

true_segments_full<-as.matrix(true_segments_full[, c("response_count", "no_response_count")])
segment_responses_full <- as.matrix(segment_responses_full[, c("response_count", "no_response_count")])
str(segment_responses_full)
str(true_segments_full)

# Define the different values of true_weight
true_weights <- c(1,10,100,300)  

#Create a dataframe to store results
probability_df <- segment_data[, c("true_gender", "true_age", "true_demo", "fraction")]

#Optimize the loglikelihood for every different 
for (w in true_weights) {
  
  optim_result <- optim(
    par = initial_beta,
    fn = loglikelihood_segments_based,
    p_z_given_s = conditional_probs_full,
    segment_responses = segment_responses_full,
    true_segments = true_segments_full,      
    true_weight = w,                               
    est_weight = 1,                            
    control = list(maxit = 1000),
    method = "L-BFGS-B",
    lower=-10,
    upper=0
  )
  
  estimated_probs <- plogis(optim_result$par)  
  
  probability_df[[paste0("estimated_prob_w", w)]] <- estimated_probs
}

#For every weight, define the differences between the estimated probability and the fraction of exposed individuals
differences <- data.frame(
  difference_w1   = probability_df$fraction - probability_df$estimated_prob_w1,
  difference_w10  = probability_df$fraction - probability_df$estimated_prob_w10,
  difference_w100 = probability_df$fraction - probability_df$estimated_prob_w100,
  difference_w300 = probability_df$fraction - probability_df$estimated_prob_w300
)

differences_long <- pivot_longer(differences, cols = everything(), 
                                 names_to = "Weight", values_to = "Difference")

# Create separate histograms for each difference
ggplot(differences_long, aes(x = Difference)) +
  geom_histogram(binwidth=0.02, fill = "skyblue", color = "black", alpha = 0.7) +
  facet_wrap(~ Weight, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms of Differences (Fraction - Estimated Probabilities)",
       x = "Difference",
       y = "Frequency")

# Print the final probability data frame
print(probability_df)

write.xlsx(probability_df, "probabilities_425736.xlsx", rowNames = FALSE)


