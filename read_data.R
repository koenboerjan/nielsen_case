read_id_graph_demos <- function() {
  data <- fread("other_data/id_graph_demos.csv")
  return (data)
}

read_panel_demos <- function() {
  fread("other_data/panel_demos.csv")
}

read_universe_estimates <- function() {
  universe_estimates <- fread("other_data/universe_estimates.csv")
  universe_estimates <- universe_estimates %>%
    mutate(
      gender_bucket = ifelse(gender_bucket == "M", 2, 1),
      age_bucket = case_when(
        age_bucket == "18-20" ~ 1,
        age_bucket == "21-24" ~ 1,
        age_bucket == "25-29" ~ 1,
        age_bucket == "30-34" ~ 1,
        age_bucket == "35-39" ~ 2,
        age_bucket == "40-44" ~ 2,
        age_bucket == "45-49" ~ 2,
        age_bucket == "50-54" ~ 3,
        age_bucket == "55-59" ~ 3,
        age_bucket == "60-64" ~ 3,
        age_bucket == "65+" ~ 4
      ) 
    ) %>% 
    group_by(gender_bucket, age_bucket, demo3_bucket) %>%
    summarise(num_persons = sum(num_persons), tot_persons = mean(tot_persons), .groups = "keep")
  
  return (universe_estimates)
}

read_campaign_details <- function() {
  fread("other_data/campaign_details.csv")
}

read_site_details <- function() {
  fread("other_data/site_details.csv")
}

combine_panel_id_graph <- function() {
  id_graph_demos <- read_id_graph_demos()
  panel_demos <- read_panel_demos()
  
  # Clean missing observations from id_graph_demos
  id_graph_demos_clean <- id_graph_demos %>%
    filter(complete.cases(.))
  
  merged_demos <- merge(id_graph_demos_clean, panel_demos, by = 'person_id', all.x = TRUE) 
  
  colnames(merged_demos) <- c('person_id', 'estimated_age', 'estimated_gender', 'estimated_demo', 'true_age', 
                                             'true_gender', 'true_demo')
  
  # Standardize data types in `demographic_groups`
  merged_demos <- merged_demos %>%
    mutate(
      true_gender = ifelse(true_gender == "male", 2, 1),
      true_demo = as.numeric(true_demo),
      true_age = case_when(
        true_age == "lt35" ~ 1,
        true_age == "gt35_lt50" ~ 2,
        true_age == "gt50_lt65" ~ 3,
        true_age == "gt65" ~ 4),
      estimated_gender = ifelse(estimated_gender == "male", 2, 1),
      estimated_demo = as.numeric(estimated_demo),
      estimated_age = case_when(
        estimated_age == "lt35" ~ 1,
        estimated_age == "gt35_lt50" ~ 2,
        estimated_age == "gt50_lt65" ~ 3,
        estimated_age == "gt65" ~ 4)
      )
  
  return(merged_demos)
}

read_only_exposures <- function() {
  # Load and combine all exposure files
  exposure_files <- list.files(path = "exposures_all/", pattern = "exposures_nlsn.*\\.csv", full.names = TRUE)
  
  # Define a function to standardize and read files
  read_and_standardize <- function(file) {
    data <- fread(file)
    if ("site_id" %in% names(data)) {
      data <- data %>%
        mutate(site_id = as.character(site_id)) # Ensure consistent data type
    }
    return(data)
  }
  
  # Combine all exposures
  all_exposures <- exposure_files %>%
    map_dfr(read_and_standardize)
  
  aggregated_exposure <- all_exposures %>%
    group_by(person_id, site_id) %>%
    summarise(
      total_exposures = sum(num_exposures, na.rm = TRUE),
      .groups = 'drop'
    )
  
  summary_exposed <- aggregated_exposure %>%
    group_by(site_id) %>%
    summarise(
      total_exposed = n(),
      .groups = 'drop'
    )
  
  return(summary_exposed)
}

read_exposures <- function(site_id_input = NA) {
  id_graph_demos <- read_id_graph_demos()
  panel_demos <- read_panel_demos()
  
  # Clean missing observations from id_graph_demos
  id_graph_demos_clean <- id_graph_demos %>%
    filter(complete.cases(.))
  
  # Load and combine all exposure files
  exposure_files <- list.files(path = "exposures_all/", pattern = "exposures_nlsn.*\\.csv", full.names = TRUE)
  
  # Define a function to standardize and read files
  read_and_standardize <- function(file) {
    data <- fread(file)
    if ("site_id" %in% names(data)) {
      data <- data %>%
        mutate(site_id = as.character(site_id)) # Ensure consistent data type
    }
    return(data)
  }
  
  # Combine all exposures
  all_exposures <- exposure_files %>%
    map_dfr(read_and_standardize)
  
  
  if (!is.na(site_id_input)) {
    site_id_vector <- strsplit(as.character(site_id_input), ",\\s*")[[1]]  # Convert to vector
    all_exposures <- all_exposures %>%
      filter(site_id %in% site_id_vector)  # Use %in% for multiple values
    # Aggregate exposure data
    aggregated_exposure <- all_exposures %>%
      group_by(person_id, site_id) %>%
      summarise(
        total_exposures = sum(num_exposures, na.rm = TRUE),
        .groups = 'drop'
      )
  } else {
    # Aggregate exposure data
  aggregated_exposure <- all_exposures %>%
    group_by(person_id) %>%
    summarise(
      total_exposures = sum(num_exposures, na.rm = TRUE),
      .groups = 'drop'
    )
  }
  
  total_count <- dim(id_graph_demos_clean)[1]
  subset_size <- 250000
  
  loopings <- floor(total_count/subset_size) + 1
  for (r in 1:loopings) {
    print(r)
    if (r == loopings) {
      print("finito")
      merged_demos_subset <- merge(id_graph_demos_clean[(subset_size*(r-1)+1): total_count,], panel_demos, by = 'person_id', all.x = TRUE) 
    } else {
      merged_demos_subset <- merge(id_graph_demos_clean[(subset_size*(r-1)+1): (subset_size*r),], panel_demos, by = 'person_id', all.x = TRUE) 
    }
    
    merged_exposure_data_subset <- merge(
      merged_demos_subset,
      aggregated_exposure,
      by = "person_id",
      all.x = TRUE
    ) %>%
      mutate(total_exposures = coalesce(total_exposures, 0)) %>%  # Replace NA with 0
      mutate(response = ifelse(total_exposures > 0, 1, 0))       # Define binary response
    
    if (!is.na(site_id_input)) {
      colnames(merged_exposure_data_subset) <- c('person_id', 'estimated_age', 'estimated_gender', 'estimated_demo', 
                                                 'true_age', 'true_gender', 'true_demo', 'site_id', 'total_exposures', 
                                                 'response')
    } else {
      colnames(merged_exposure_data_subset) <- c('person_id', 'estimated_age', 'estimated_gender', 'estimated_demo', 
                                                 'true_age', 'true_gender', 'true_demo', 'total_exposures', 
                                                 'response')
    }
    dataset_true_only <- merged_exposure_data_subset %>% filter(!is.na(true_gender))
    
    segments_response_true <- dataset_true_only %>%
      group_by(true_gender, true_age, true_demo) %>%
      summarise(
        response_count = sum(response),
        no_response_count = sum(1 - response),
        .groups = 'drop'
      )
    
    dataset <- merged_exposure_data_subset %>% filter(is.na(true_gender))
    
    segments_response <- dataset %>%
      group_by(estimated_gender, estimated_age, estimated_demo) %>%
      summarise(
        response_count = sum(response),
        no_response_count = sum(1 - response),
        .groups = 'drop'
      )
    
    if (r == 1) {
      segments_response_total_true <- segments_response_true
      segments_response_total <- segments_response
    } else {
      segments_response_total_true <- merge(segments_response_total_true, segments_response_true, by = c('true_gender', 'true_age', 'true_demo'), all.x = TRUE) %>%
        mutate(response_count = ifelse(!is.na(response_count.x), response_count.x, 0) + ifelse(!is.na(response_count.y), response_count.y, 0)) %>%
        mutate(no_response_count = ifelse(!is.na(no_response_count.x), no_response_count.x, 0) + ifelse(!is.na(no_response_count.y), no_response_count.y, 0))
      segments_response_total_true <- segments_response_total_true[, c('true_gender', 'true_age', 'true_demo', 'response_count', 'no_response_count')]
      segments_response_total <- merge(segments_response_total, segments_response, by = c('estimated_gender', 'estimated_age', 'estimated_demo'), all.x = TRUE) %>%
        mutate(response_count = ifelse(!is.na(response_count.x), response_count.x, 0) + ifelse(!is.na(response_count.y), response_count.y, 0)) %>%
        mutate(no_response_count = ifelse(!is.na(no_response_count.x), no_response_count.x, 0) + ifelse(!is.na(no_response_count.y), no_response_count.y, 0))
      segments_response_total <- segments_response_total[, c('estimated_gender', 'estimated_age', 'estimated_demo', 'response_count', 'no_response_count')]
    }
  }
  
  # Standardize data types in `demographic_groups`
  segments_response_total_true <- segments_response_total_true %>%
    mutate(
      true_gender = ifelse(true_gender == "male", 2, 1),
      true_demo = as.numeric(true_demo),
      true_age = case_when(
        true_age == "lt35" ~ 1,
        true_age == "gt35_lt50" ~ 2,
        true_age == "gt50_lt65" ~ 3,
        true_age == "gt65" ~ 4
      )
    )
  # Standardize data types in `demographic_groups`
  segments_response_total <- segments_response_total %>%
    mutate(
      estimated_gender = ifelse(estimated_gender == "male", 2, 1),
      estimated_demo = as.numeric(estimated_demo),
      estimated_age = case_when(
        estimated_age == "lt35" ~ 1,
        estimated_age == "gt35_lt50" ~ 2,
        estimated_age == "gt50_lt65" ~ 3,
        estimated_age == "gt65" ~ 4
      )
    )
  
  print(paste("Total amount of exposed indivduals is: ", sum(segments_response_total$response_count, segments_response_total_true$response_count)))
  return(list(segments_response_total_true, segments_response_total))
}
