library(readxl)
probabilities_389882<-read_excel("probabilities_389882_new.xlsx")
probabilities_458088<-read_excel("probabilities_458088.xlsx")

differences1 <- data.frame(
  difference_w1   = probabilities_389882$fraction -  probabilities_389882$estimated_prob_w1,
  difference_w10  =  probabilities_389882$fraction -  probabilities_389882$estimated_prob_w10,
  difference_w100 =  probabilities_389882$fraction -  probabilities_389882$estimated_prob_w100,
  difference_w300 =  probabilities_389882$fraction -  probabilities_389882$estimated_prob_w300
)

differences_long1 <- pivot_longer(differences1, cols = everything(), 
                                 names_to = "Weight", values_to = "Difference")

# Create separate histograms for each difference
ggplot(differences_long1, aes(x = Difference)) +
  geom_histogram(binwidth = 0.02, fill = "skyblue", color = "black", alpha = 0.7) +
  facet_wrap(~ Weight, scales = "fixed") +
  theme_minimal() +
  labs(title = "Histograms of Differences -Campaign 389882",
       x = "Difference",
       y = "Count")

differences2 <- data.frame(
  difference_w1   = probabilities_458088$fraction -  probabilities_458088$estimated_prob_w1,
  difference_w10  =  probabilities_458088$fraction -  probabilities_458088$estimated_prob_w10,
  difference_w100 =  probabilities_458088$fraction -  probabilities_458088$estimated_prob_w100,
  difference_w300 =  probabilities_458088$fraction -  probabilities_458088$estimated_prob_w300
)

differences_long2 <- pivot_longer(differences2, cols = everything(), 
                                 names_to = "Weight", values_to = "Difference")

# Create separate histograms for each difference
ggplot(differences_long2, aes(x = Difference)) +
  geom_histogram(binwidth = 0.02, fill = "skyblue", color = "black", alpha = 0.7) +
  facet_wrap(~ Weight, scales = "fixed") +
  theme_minimal() +
  labs(title = "Histograms of Differences- Campaign 458088",
       x = "Difference",
       y = "Count")