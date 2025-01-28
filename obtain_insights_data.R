create_bar_chart <- function(dataset, demographics) {
  for (i in 1: length(demographics)) {
    if (demographics[i] == 'gender') {
      print(ggplot(dataset, aes(x = estimated_gender)) +
        geom_bar(fill = "steelblue") +
        theme_minimal())
    }
    if (demographics[i] == 'age') {
      print(ggplot(dataset, aes(x = estimated_age)) +
              geom_bar(fill = "steelblue") +
              theme_minimal())
    }
    if (demographics[i] == 'demo') {
      print(ggplot(dataset, aes(x = estimated_demo)) +
              geom_bar(fill = "steelblue") +
              theme_minimal())
    }
  }
}

insight_true_demos <- function(dataset, demographics) {
  # Clean missing true values from dataset
  dataset_true_values_known <- dataset %>%
    filter(complete.cases(.))
  for (i in 1: length(demographics)) {
    if (demographics[i] == 'gender') {
      print(ggplot(dataset_true_values_known, aes(x = true_gender)) +
              geom_bar(fill = "maroon") +
              theme_minimal())
    }
    if (demographics[i] == 'age') {
      print(ggplot(dataset_true_values_known, aes(x = true_age)) +
              geom_bar(fill = "maroon") +
              theme_minimal())
    }
    if (demographics[i] == 'demo') {
      print(ggplot(dataset_true_values_known, aes(x = true_demo)) +
              geom_bar(fill = "maroon") +
              theme_minimal())
    }
  }
}
