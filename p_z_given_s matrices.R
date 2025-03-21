p_z_given_s_gender_without_prior <- compute_p_z_given_s(merged_demos, "gender")
p_z_given_s_gender_with_prior <- compute_p_z_given_s_including_prior(merged_demos, "gender")
p_z_given_s_age_without_prior <- compute_p_z_given_s(merged_demos, "age")
p_z_given_s_age_with_prior <- compute_p_z_given_s_including_prior(merged_demos, "age")
p_z_given_s_demo_without_prior <- compute_p_z_given_s(merged_demos, "demo")
p_z_given_s_demo_with_prior <- compute_p_z_given_s_including_prior(merged_demos, "demo")

# Define the probability matrices

p_z_given_s_gender_without_prior <- matrix(c(
  0.7463356, 0.2536644,
  0.4714094, 0.5285906
), nrow = 2, byrow = TRUE)

p_z_given_s_gender_with_prior <- matrix(c(
  0.6492846, 0.3507154,
  0.3594495, 0.6405505
), nrow = 2, byrow = TRUE)

p_z_given_s_age_without_prior <- matrix(c(
  0.6161240, 0.1595249, 0.07746539, 0.1468857,
  0.2523317, 0.4872976, 0.09237073, 0.1680000,
  0.2914582, 0.2112857, 0.34592698, 0.1513291,
  0.2804664, 0.1563858, 0.05859244, 0.5045553
), nrow = 4, byrow = TRUE)

p_z_given_s_age_with_prior <- matrix(c(
  0.5233834, 0.1595561, 0.12749693, 0.1895636,
  0.2002178, 0.4552587, 0.14200559, 0.2025180,
  0.2023501, 0.1727152, 0.46532005, 0.1596147,
  0.2085784, 0.1369364, 0.08442486, 0.5700603
), nrow = 4, byrow = TRUE)

p_z_given_s_age_with_prior <- matrix(c(
  0.85, 0.05, 0.05 ,0.05,
  0.05, 0.85, 0.05, 0.05,
  0.05, 0.05, 0.85, 0.05,
  0.05, 0.05, 0.05, 0.85
), nrow = 4, byrow = TRUE)


p_z_given_s_demo_without_prior <- matrix(c(
  0.7627596, 0.07846589, 0.10935700, 0.01623943, 0.03317812,
  0.3637739, 0.47312031, 0.09831726, 0.02911704, 0.03567152,
  0.3904692, 0.07005661, 0.48369282, 0.01160085, 0.04418047,
  0.4159415, 0.11573268, 0.10616081, 0.30804038, 0.05412461,
  0.6665010, 0.09567594, 0.11108350, 0.04125249, 0.08548708
), nrow = 5, byrow = TRUE)

p_z_given_s_demo_with_prior <- matrix(c(
  0.7273641, 0.1308560, 0.08402524, 0.04067852, 0.01707617,
  0.2662787, 0.6056546, 0.05798741, 0.05598637, 0.01409292,
  0.4079969, 0.1280170, 0.40722899, 0.03184124, 0.02491581,
  0.2696967, 0.1312345, 0.05546341, 0.52466400, 0.01894143,
  0.6183730, 0.1552390, 0.08304208, 0.10053796, 0.04280798
), nrow = 5, byrow = TRUE)



