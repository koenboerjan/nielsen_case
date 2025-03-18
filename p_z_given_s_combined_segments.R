p_z_given_s_full<-function(){
  #Read the universe dataset, as well as the estimated and the true demographics dataset
  universe_estimates <- read_universe_estimates()
  id_graph_demos<-read_id_graph_demos()
  panel_demos<-read_panel_demos()
  print("Reading Done")
  
  #Remove the incomplete observations in the estimated demos dataset
  id_graph_demos_clean <- id_graph_demos %>%
    filter(complete.cases(.))
  
  #Merge the true and estimated demographics datasets
  colnames(id_graph_demos_clean)=c("person_id","estimated_age","estimated_gender","estimated_demo")
  colnames(panel_demos)=c("person_id","true_age","true_gender","true_demo")
  merged_demos <- merge(id_graph_demos_clean, panel_demos, by = 'person_id', all.x = TRUE)
  
  #Remove observations with unknown true demographics
  merged_demos <- merged_demos %>% filter(!is.na(true_demo) & !is.na(true_age) & !is.na(true_gender))
  
  #Change the from categorical demographic variables to numeric variables
  merged_demos <- merged_demos %>%
    mutate(
      estimated_gender = ifelse(estimated_gender == "male", 1, 0),
      estimated_demo = as.numeric(estimated_demo),
      estimated_age = case_when(
        estimated_age == "lt35" ~ 1,
        estimated_age == "gt35_lt50" ~ 2,
        estimated_age == "gt50_lt65" ~ 3,
        estimated_age == "gt65" ~ 4
      ),
      true_gender = ifelse(true_gender == "male", 1, 0),
      true_demo = as.numeric(true_demo),
      true_age = case_when(
        true_age == "lt35" ~ 1,
        true_age == "gt35_lt50" ~ 2,
        true_age == "gt50_lt65" ~ 3,
        true_age == "gt65" ~ 4
      )
    )
  print("merging done")
  
  #Calculate the prior probability of the true segment based on the universe estimates
  p_z <- universe_estimates %>%
    mutate(probability = num_persons / tot_persons)
  
  # Create all possible combinations for true and estimated demographic
  true_combinations <- expand.grid(true_gender = 0:1, true_age = 1:4, true_demo = 1:5)
  estimated_combinations <- expand.grid(estimated_gender = 0:1, estimated_age = 1:4, estimated_demo = 1:5)
  
  #Create all the different combinations of true and estimated demographics
  all_combinations <- merge(true_combinations, estimated_combinations, by = NULL)
  
  # Compute the conditional probability of estimated segment given the true one P(S | Z) based on the merged demos
  p_s_given_z <- merged_demos %>%
    group_by(estimated_gender, estimated_age, estimated_demo, true_gender, true_age, true_demo) %>%
    summarise(count = n(), .groups = 'drop') %>%
    replace_na(list(count = 0)) %>%
    group_by(true_gender, true_age, true_demo) %>%
    mutate(probability = count / sum(count)) 
  
  # Make sure that all combinations exist in our dataset, so that our final matrix will be 40 x 40
  p_s_given_z <- left_join(all_combinations, p_s_given_z,                           
                           by = c("true_gender", "true_age", "true_demo", "estimated_gender", "estimated_age", "estimated_demo")) %>%
    replace_na(list(count = 0, probability = 0))
  print(p_s_given_z)
  
  
  #Place the calculated probabilities into a square matrix
  dim_matrix <- sqrt(length(p_s_given_z$probability))
  p_s_given_z_matrix <- matrix(p_s_given_z$probability, nrow = dim_matrix, ncol = dim_matrix)
  
  # Compute the matrix with the probabilities P(S,Z)=P(S|Z)*P(Z) 
  p_s_and_z <- p_s_given_z_matrix * p_z$probability
  #P(S) is calculating by summing across the columns of P(S,Z)
  p_s <- colSums(p_s_and_z)
  #Store the inverses of P(S) into a matrix in order to carry out the division in the next step
  p_s_matrix <- t(matrix(rep((p_s**-1),dim_matrix), ncol = dim_matrix))
  
  #Compute the matrix with probabilities P(Z|S)=P(S,Z)/P(S)
  p_z_given_s <- t(p_s_and_z * p_s_matrix)
  
  #If a row has NA values, then for that row we assign equal probability (0.025) for every true value Z
  p_z_given_s[is.na(p_z_given_s)] <- 0.025
  
  return(p_z_given_s)
  
}


p_z_given_s_full()