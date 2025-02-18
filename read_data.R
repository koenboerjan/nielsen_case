
# ------------------------ Read Data Functions ------------------------

read_id_graph_demos <- function() {
  fread("other_data/id_graph_demos.csv", nrow = 5000000)
}

read_panel_demos <- function() {
  fread("other_data/panel_demos.csv")
}

read_universe_estimates <- function() {
  fread("other_data/universe_estimates.csv")
}

read_site_details <- function() {
  fread("other_data/site_details.csv")
}

# ------------------------ Read and Process Exposure Data ------------------------

read_exposures <- function(site_id_input = NA, platform_input = NA) {
  
  # Load demographic data
  id_graph_demos <- fread("other_data/id_graph_demos.csv") %>%
    filter(complete.cases(.))  # Remove incomplete records
  
  panel_demos <- fread("other_data/panel_demos.csv")
  
  # Load and combine all exposure files
  exposure_files <- list.files(path = "exposures_all/", pattern = "exposures_nlsn.*\\.csv", full.names = TRUE)
  
  read_and_standardize <- function(file) {
    fread(file, nrow = 1000000) %>%
      mutate(site_id = as.character(site_id))  # Ensure consistent site_id type
  }
  
  all_exposures <- map_dfr(exposure_files, read_and_standardize)
  
  # Apply site and platform filters if specified
  if (!is.na(site_id_input)) {
    site_id_vector <- strsplit(as.character(site_id_input), ",\\s*")[[1]]
    all_exposures <- all_exposures %>% filter(site_id %in% site_id_vector)
  }
  
  if (!is.na(platform_input)) {
    platform_vector <- strsplit(as.character(platform_input), ",\\s*")[[1]]
    all_exposures <- all_exposures %>% filter(platform %in% platform_vector)
  }
  
  # Aggregate exposure data by person_id, site_id, and platform
  aggregated_exposure <- all_exposures %>%
    group_by(person_id, site_id, platform) %>%
    summarise(
      total_exposures = sum(num_exposures, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Merge demographic data
  merged_demos <- merge(id_graph_demos, panel_demos, by = 'person_id', all.x = TRUE)
  
  # Merge with aggregated exposure data
  merged_exposure_data <- merge(
    merged_demos,
    aggregated_exposure,
    by = "person_id",
    all.x = TRUE
  ) %>%
    mutate(
      total_exposures = coalesce(total_exposures, 0),  # Replace NA with 0
      response = ifelse(total_exposures > 0, 1, 0)     # Define binary response
    )
  
  # Standardize column names
  colnames(merged_exposure_data) <- c(
    'person_id', 'estimated_age', 'estimated_gender', 'estimated_demo', 'true_age', 
    'true_gender', 'true_demo', 'site_id', 'platform', 'total_exposures', 'response'
  )
  
  # Convert demographic values to numeric categories
  merged_exposure_data <- merged_exposure_data %>%
    mutate(
      estimated_gender = ifelse(estimated_gender == "male", 1, 0),
      true_gender = ifelse(true_gender == "male", 1, 0),
      estimated_age = case_when(
        estimated_age == "18-24" ~ 0,
        estimated_age == "25-34" ~ 1,
        estimated_age == "35-44" ~ 2,
        estimated_age == "45-64" ~ 3,
        estimated_age == "gt65" ~ 4
      ),
      true_age = case_when(
        true_age == "18-24" ~ 0,
        true_age == "25-34" ~ 1,
        true_age == "35-44" ~ 2,
        true_age == "45-64" ~ 3,
        true_age == "gt65" ~ 4
      )
    )
  
  return(merged_exposure_data)
}
