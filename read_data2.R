read_id_graph_demos <- function(){
  return(read_parquet("other_data/id_graph_demos_v3.parquet"))
}
read_panel_demos <- function(){
  return(read_parquet("other_data/panel_demos_v3.parquet"))
}
read_universe_estimates <- function(){
  return(read_parquet("other_data/universe_estimates.parquet"))
}
read_campaign_details <- function(){
  return(read_parquet("other_data/campaign_details.parquet"))
}
read_site_details <- function(){
  return(read_parquet("other_data/site_details.parquet"))
}


# Function to process exposure files in chunks
process_exposure_file <- function(file) {
  data <- read_parquet(file, as_data_frame = FALSE)  # Read as Arrow Table
  df <- as.data.table(data)  # Convert to data.table
  
  # Ensure consistent data types
  if ("site_id" %in% names(df)) {
    df[, site_id := as.character(site_id)]
  }
  
  # Aggregate per file
  df <- df[, .(total_exposures = sum(num_exposures, na.rm = TRUE)), by = person_id]
  
  return(df)
}

# Function 1: Load and clean demographic data
load_demographics <- function() {
  id_graph_demos <- read_id_graph_demos
  panel_demos <- read_panel_demos
  
  setDT(id_graph_demos)
  id_graph_demos_clean <- id_graph_demos[complete.cases(id_graph_demos)]
  setDT(panel_demos)
  
  return(list(id_graph_demos_clean = id_graph_demos_clean, panel_demos = panel_demos))
}

# Function 2: Process a single parquet file
process_exposure_file <- function(file) {
  data <- read_parquet(file)
  setDT(data)
  
  if (!"person_id" %in% names(data) || !"num_exposures" %in% names(data)) {
    return(NULL)  # Skip invalid files
  }
  
  # Aggregate exposures per person
  chunk_agg <- data[, .(total_exposures = sum(num_exposures, na.rm = TRUE)), by = person_id]
  return(chunk_agg)
}

# Function 3: Read and aggregate all exposure data
read_and_aggregate_exposures <- function() {
  exposure_files <- list.files(path = "exposures_all/", pattern = "exposures_nlsn.*\\.parquet", full.names = TRUE)
  
  pb <- progress_bar$new(
    format = "Processing Files [:bar] :percent ETA: :eta",
    total = length(exposure_files),
    clear = FALSE,
    width = 60
  )
  
  aggregated_exposure <- data.table(person_id = integer(), total_exposures = numeric())
  
  for (file in exposure_files) {
    chunk_data <- process_exposure_file(file)
    
    if (!is.null(chunk_data)) {
      aggregated_exposure <- rbindlist(list(aggregated_exposure, chunk_data), use.names = TRUE, fill = TRUE)
      
      # Re-aggregate to avoid duplicate entries
      aggregated_exposure <- aggregated_exposure[, .(total_exposures = sum(total_exposures)), by = person_id]
    }
    
    pb$tick()
  }
  
  return(aggregated_exposure)
}

# Function 4: Merge all data together
merge_all_data <- function(id_graph_demos_clean, panel_demos, aggregated_exposure) {
  pb_merge <- progress_bar$new(
    format = "Merging Data [:bar] :percent ETA: :eta",
    total = 3,
    clear = FALSE,
    width = 60
  )
  
  setkey(id_graph_demos_clean, person_id)
  setkey(panel_demos, person_id)
  setkey(aggregated_exposure, person_id)
  
  pb_merge$tick()
  merged_demos <- panel_demos[id_graph_demos_clean, nomatch = 0]
  
  pb_merge$tick()
  merged_exposure_data <- aggregated_exposure[merged_demos, nomatch = 0]
  
  pb_merge$tick()
  merged_exposure_data[, total_exposures := fifelse(is.na(total_exposures), 0, total_exposures)]
  merged_exposure_data[, response := fifelse(total_exposures > 0, 1, 0)]
  
  colnames(merged_exposure_data) <- c('person_id', 'estimated_age', 'estimated_gender', 'estimated_demo', 
                                      'true_age', 'true_gender', 'true_demo', 'total_exposures', 'response')
  
  return(merged_exposure_data)
}

# Main function: Calls all the smaller functions
read_exposures <- function() {
  message("Loading demographic data...")
  demographics <- load_demographics()
  
  message("Reading and aggregating exposures...")
  aggregated_exposure <- read_and_aggregate_exposures()
  
  message("Merging all data together...")
  real_dataset <- merge_all_data(demographics$id_graph_demos_clean, demographics$panel_demos, aggregated_exposure)
  
  return(real_dataset)
}

real_dataset <- read_exposures()