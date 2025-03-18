
# ----------------------------- Compute Log Likelihood for Segments -----------------------------
loglikelihood_segments_based <- function(beta, p_z_given_s, segment_responses, true_segments, true_weight, est_weight) {
  log_est <- 0
  log_true <- 0
  segment_count <- dim(p_z_given_s)[1]
  
  # Compute likelihood for estimated demographics
  for (i in 1:segment_count) {
    for (j in 1:2) {
      response_j <- ifelse(j == 1, 1, 0)
      value <- 0
      
      for (s in 1:segment_count) {
        #Calculate the probability of response based on the logit function
        e_z_beta <- plogis(beta[s])
        value <- value + p_z_given_s[i, s] * (e_z_beta^response_j) * ((1 - e_z_beta)^(1 - response_j))
        
      }
      #For every possible combination of segment z and binary variable y_j we count how many
      #Observations fall into that combination and multiply it with the log-likelihood of this 
      #combination. This effectively reduces the computational expense
      log_est <- log_est + segment_responses[i, j] * log(max(value_segment, 1e-10))
    }
  }
  
  # Compute likelihood for true demographics
  for (z in 1:segment_count) {
    for (j in 1:2) {
      response_j <- ifelse(j == 1, 1, 0)
      #Calculate the probability of response based on the logit function
      e_z_beta <- plogis(beta[z])  
      
      # Compute the conditional probability p(y_j | z)
      p_y_given_z <- (e_z_beta^response_j) * ((1 - e_z_beta)^(1 - response_j))
      
      # Log likelihood contribution for this (y_j, z) pair
      log_true <- log_true + true_segments[z, j] * log(max(p_y_given_z, 1e-10))
      
    }
  }
  #return the weighted log likelihood
  return(-true_weight * log_true - est_weight * log_est)
  
  
}
#----------------------------------------------------------------------------------------------------------------------------------------#
#Load libraries and functions
prepare_setup <- function() {
  source("requirements.R")
  install_requirements()
}


prepare_setup()


#Read and merge the datasets
#BEAR IN MIND: Since we investigate every exposure dataset separately, we only include the dataset 
#of one campaign in the exposures_all file every time
real_dataset<-read_exposures()


##########################################################################################
#COMBINED SEGMENTS ESTIMATIONS
#Now estimate the conditional probabilities for the full set of demographic combinations

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
  
  estimated_probabilities <- plogis(optim_result$par)  
  
  probability_df[[paste0("estimated_p_w", w)]] <- estimated_probabilities
}

#For every weight, define the differences between the estimated probability and the fraction of exposed individuals
differences <- data.frame(
  difference_w1   = probability_df$fraction - probability_df$estimated_p_w1,
  difference_w10  = probability_df$fraction - probability_df$estimated_p_w10,
  difference_w100 = probability_df$fraction - probability_df$estimated_p_w100,
  difference_w300 = probability_df$fraction - probability_df$estimated_p_w300
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