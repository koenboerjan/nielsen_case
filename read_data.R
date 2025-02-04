

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
      gender_bucket = ifelse(gender_bucket == "M", 1, 0),
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

read_exposures <- function(site_id_input = NA, platform_input = NA) {
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
  }
  
  if (!is.na(platform_input)) {
    platform_vector <- strsplit(as.character(platform_input), ",\\s*")[[1]]  # Convert to vector
    all_exposures <- all_exposures %>%
      filter(platform %in% platform_vector)
  }
  
  # Aggregate exposure data
  aggregated_exposure <- all_exposures %>%
    group_by(person_id, site_id, platform) %>%
    summarise(
      total_exposures = sum(num_exposures, na.rm = TRUE),
      .groups = 'drop'
    )
  
  merged_demos <- merge(id_graph_demos_clean, panel_demos, by = 'person_id', all.x = TRUE) 
  
  merged_exposure_data <- merge(
    merged_demos,
    aggregated_exposure,
    by = "person_id",
    all.x = TRUE  # Retain all rows from panel_demos
  ) %>%
    mutate(total_exposures = coalesce(total_exposures, 0)) %>%  # Replace NA with 0
    mutate(response = ifelse(total_exposures > 0, 1, 0))       # Define binary response


  
  colnames(merged_exposure_data) <- c('person_id', 'estimated_age', 'estimated_gender', 'estimated_demo', 'true_age', 
                                      'true_gender', 'true_demo', 'site_id', 'platform', 'total_exposures', 'response')
  
  # Standardize data types in `demographic_groups`
  merged_exposure_data <- merged_exposure_data %>%
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
  return(merged_exposure_data)
}

