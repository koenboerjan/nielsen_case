

read_id_graph_demos <- function() {
  read.csv("other_data/id_graph_demos_short.csv")
}

read_panel_demos <- function() {
  read.csv("other_data/panel_demos.csv")
}

read_universe_estimates <- function() {
  read.csv("other_data/universe_estimates.csv")
}

read_campaign_details <- function() {
  read.csv("other_data/campaign_details.csv")
}

read_site_details <- function() {
  read.csv("other_data/site_details.csv")
}

read_exposures <- function() {
  id_graph_demos <- read_id_graph_demos()
  panel_demos <- read_panel_demos()
  
  # Clean missing observations from id_graph_demos
  id_graph_demos_clean <- id_graph_demos %>%
    filter(complete.cases(.))
  
  # Load and combine all exposure files
  exposure_files <- list.files(path = "exposures_all/", pattern = "exposures_nlsn.*\\.csv", full.names = TRUE)
  
  # Define a function to standardize and read files
  read_and_standardize <- function(file) {
    data <- read.csv(file)
    if ("site_id" %in% names(data)) {
      data <- data %>%
        mutate(site_id = as.character(site_id)) # Ensure consistent data type
    }
    return(data)
  }
  
  # Combine all exposures
  all_exposures <- exposure_files %>%
    map_dfr(read_and_standardize)
  
  # Aggregate exposure data
  aggregated_exposure <- all_exposures %>%
    group_by(person_id) %>%
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
                                      'true_gender', 'true_demo', 'total_exposures', 'response')
  return(merged_exposure_data)
}

