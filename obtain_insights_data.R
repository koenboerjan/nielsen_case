create_summary_estimates <- function(dataset, demographics) {
  for (i in 1: length(demographics)) {
    if (demographics[i] == 'gender') {
      print(dataset %>%
              group_by(estimated_gender) %>%
              summarise(count = n(), .groups = 'drop'))
      print(ggplot(dataset, aes(x = estimated_gender)) +
        geom_bar(fill = "steelblue") +
        theme_minimal())
    }
    if (demographics[i] == 'age') {
      print(dataset %>%
              group_by(estimated_age) %>%
              summarise(count = n(), .groups = 'drop'))
      print(ggplot(dataset, aes(x = estimated_age)) +
              geom_bar(fill = "steelblue") +
              theme_minimal())
    }
    if (demographics[i] == 'demo') {
      print(dataset %>%
              group_by(estimated_demo) %>%
              summarise(count = n(), .groups = 'drop'))
      print(ggplot(dataset, aes(x = estimated_demo)) +
              geom_bar(fill = "steelblue") +
              theme_minimal())
    }
  }
}

create_summary_true <- function(dataset, demographics) {
  # Clean missing true values from dataset
  dataset_true_values_known <- dataset %>%
    filter(complete.cases(.))
  
  for (i in 1: length(demographics)) {
    if (demographics[i] == 'gender') {
      print(dataset_true_values_known %>%
              group_by(true_gender) %>%
              summarise(count = n(), .groups = 'drop'))
      print(ggplot(dataset_true_values_known, aes(x = true_gender)) +
              geom_bar(fill = "maroon") +
              theme_minimal())
    }
    if (demographics[i] == 'age') {
      print(dataset_true_values_known %>%
              group_by(true_age) %>%
              summarise(count = n(), .groups = 'drop'))
      print(ggplot(dataset_true_values_known, aes(x = true_age)) +
              geom_bar(fill = "maroon") +
              theme_minimal())
    }
    if (demographics[i] == 'demo') {
      print(dataset_true_values_known %>%
              group_by(true_demo) %>%
              summarise(count = n(), .groups = 'drop'))
      print(ggplot(dataset_true_values_known, aes(x = true_demo)) +
              geom_bar(fill = "maroon") +
              theme_minimal())
    }
  }
}

create_bubble_chart <- function(dataset) {
  segment_sizes <- dataset %>%
    group_by(estimated_gender, estimated_age) %>%
    summarise(count = n(), .groups = 'drop')
  
  print(ggplot(segment_sizes, aes(x = estimated_gender, y = estimated_age, size = count)) +
    geom_point(alpha = 0.6) + 
    scale_size(range = c(3, 15)) +
    theme_minimal() +
    labs(title = "Age versus gender",
         x = "Gender",
         y = "Age",
         size = "Individual count"))
  
  segment_sizes <- dataset %>%
    group_by(estimated_gender, estimated_demo) %>%
    summarise(count = n(), .groups = 'drop')
  
  print(ggplot(segment_sizes, aes(x = estimated_gender, y = estimated_demo, size = count)) +
    geom_point(alpha = 0.6) + 
    scale_size(range = c(3, 15)) +
    theme_minimal() + 
    labs(title = "Demo versus gender",
         x = "Gender",
         y = "Demo",
         size = "Individual count"))
  
  segment_sizes <- dataset %>%
    group_by(estimated_demo, estimated_age) %>%
    summarise(count = n(), .groups = 'drop')
  
  print(ggplot(segment_sizes, aes(x = estimated_demo, y = estimated_age, size = count)) +
    geom_point(alpha = 0.6) + 
    scale_size(range = c(3, 15)) +
    theme_minimal() +
    labs(title = "Age versus demo",
         x = "Demo",
         y = "Age",
         size = "Individual count"))
  
}

insight_exposed_individuals <- function(dataset, demographics) {
  only_exposed_individuals <- dataset %>%
    filter(response == 1)
  
  create_summary_estimates(only_exposed_individuals, demographics)
  create_summary_true(only_exposed_individuals, demographics)
}
