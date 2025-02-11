data_preprocessing <- function(file_name, filter){
  
  data <- read_excel(file_name)
  filter_data <- subset(data, data["Name"] == filter)
  # Replace -999.0 with NA
  data_cleaned <- mutate(filter_data, across(I:XII, ~ ifelse(. == -999.0, NA, .)))
  if(filter == "Tashkent") data_cleaned <- data_cleaned %>% filter(Years > 1893)

  sorted_data <- arrange(data_cleaned, Years)
  
  sorted_data[, c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII")] <- lapply(
    sorted_data[, c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII")], 
    function(x) as.numeric(as.character(x)))
  
  sorted_data$Average <- rowMeans(sorted_data[, c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII")], na.rm = TRUE)
  print(sorted_data)
  df <- select(sorted_data, Years, Average)
  df$YearInterval <- df$Years - min(df$Years)+1
  
  return(df)
}


get_exceedances  <- function(df, threshold) {
  
  # Identify years where the average exceeds the mean
  years_above_mean <- df$YearInterval[!is.na(df$Average) & df$Average > threshold]
  
  cat("Years where the average exceeds the mean:", length(years_above_mean), "\n")
  
  return(years_above_mean)
}

get_exceedance_cumsum  <- function(df, threshold) {
  
  exceedances_df <- df %>%
    mutate(Exceedance = ifelse(is.na(Average), 0, ifelse(Average > threshold, 1, 0)),
           CumulativeCases = cumsum(Exceedance))  # Cumulative sum of exceedance cases
  
  
  exceedances_df <- exceedances_df %>%
    mutate(Increment = c(1, diff(CumulativeCases))) %>%
    filter(Increment > 0)
  
  return(exceedances_df)
}

get_observation_period <- function(df){
  
  return(length(df$Years))
}

get_threshold <- function(df){
  
  # mean of the average values - threshold
  mean_average <- mean(df$Average, na.rm = TRUE)
  cat("Mean of the average values:", round(mean_average), "\n")
  
  return(round(mean_average))
}


display_statistics <- function(samples){
  hist(samples[, 1], main = "Posterior of Alpha", xlab = "alpha", col = "lightblue", border = "darkblue")
  hist(samples[, 2], main = "Posterior of Sigma", xlab = "sigma", col = "lightgreen", border = "darkgreen")
  
  mean_alpha <- mean(samples[, 1])
  mean_sigma <- mean(samples[, 2])
  
  cat("Post. mean alpha:", mean_alpha, "\n")
  cat("Post. mean sigma:", mean_sigma, "\n")
  
  # standard deviations
  sd_alpha <- sd(samples[, "alpha"])
  sd_sigma <- sd(samples[, "sigma"])
  
  # 95% credible intervals
  ci_alpha <- quantile(samples[, "alpha"], probs = c(0.025, 0.975))
  ci_sigma <- quantile(samples[, "sigma"], probs = c(0.025, 0.975))
  
  # Print the results
  cat("Standard deviation for alpha:", sd_alpha, "\n")
  cat("95% credible interval for alpha:", ci_alpha, "\n")
  cat("Standard deviation for sigma:", sd_sigma, "\n")
  cat("95% credible interval for sigma:", ci_sigma, "\n")
  
  return(c(mean_alpha, mean_sigma))
  
}


display_statistics2 <- function(samples){
  # Histograms for the posterior distributions
  hist(samples[, 1], main = "Posterior of Alpha1", xlab = "alpha1", col = "lightblue", border = "darkblue")
  hist(samples[, 2], main = "Posterior of Sigma1", xlab = "sigma1", col = "lightgreen", border = "darkgreen")
  hist(samples[, 3], main = "Posterior of Alpha2", xlab = "alpha2", col = "lightcoral", border = "darkred")
  hist(samples[, 4], main = "Posterior of Sigma2", xlab = "sigma2", col = "lightyellow", border = "darkorange")
  hist(samples[, 5], main = "Posterior of Tau", xlab = "tau", col = "lightgray", border = "black")
  
  # Mean of the posterior distributions
  mean_alpha1 <- mean(samples[, 1])
  mean_sigma1 <- mean(samples[, 2])
  mean_alpha2 <- mean(samples[, 3])
  mean_sigma2 <- mean(samples[, 4])
  mean_tau <- mean(samples[, 5])
  median_tau <- median(samples[, 5])

  cat("Post. mean alpha1:", mean_alpha1, "\n")
  cat("Post. mean sigma1:", mean_sigma1, "\n")
  cat("Post. mean alpha2:", mean_alpha2, "\n")
  cat("Post. mean sigma2:", mean_sigma2, "\n")
  cat("Post. mean tau:", mean_tau, "\n")
  cat("Post. median tau:", median_tau, "\n")
  
  sd_alpha1 <- sd(samples[, "alpha1"])
  sd_sigma1 <- sd(samples[, "sigma1"])
  sd_alpha2 <- sd(samples[, "alpha2"])
  sd_sigma2 <- sd(samples[, "sigma2"])
  sd_tau <- sd(samples[, "tau"])

  # 95% credible intervals
  ci_alpha1 <- quantile(samples[, "alpha1"], probs = c(0.025, 0.975))
  ci_sigma1 <- quantile(samples[, "sigma1"], probs = c(0.025, 0.975))
  ci_alpha2 <- quantile(samples[, "alpha2"], probs = c(0.025, 0.975))
  ci_sigma2 <- quantile(samples[, "sigma2"], probs = c(0.025, 0.975))
  ci_tau <- quantile(samples[, "tau"], probs = c(0.025, 0.975))
  
  # Print the results
  cat("Standard deviation for alpha1:", sd_alpha1, "\n")
  cat("95% credible interval for alpha1:", ci_alpha1, "\n")
  cat("Standard deviation for sigma1:", sd_sigma1, "\n")
  cat("95% credible interval for sigma1:", ci_sigma1, "\n")
  
  cat("Standard deviation for alpha2:", sd_alpha2, "\n")
  cat("95% credible interval for alpha2:", ci_alpha2, "\n")
  cat("Standard deviation for sigma2:", sd_sigma2, "\n")
  cat("95% credible interval for sigma2:", ci_sigma2, "\n")
  
  cat("Standard deviation for tau:", sd_tau, "\n")
  cat("95% credible interval for tau:", ci_tau, "\n")
  
  return(c(mean_alpha1, mean_sigma1, mean_alpha2, mean_sigma2, mean_tau))
  
}

approximate_missing <- function(vec) {
  for (i in 1:length(vec)) {
    if (is.na(vec[i])) {
      left_index <- max(1, i - 1)
      right_index <- min(length(vec), i + 1)
      left_value <- vec[left_index]
      right_value <- vec[right_index]
      if (!is.na(left_value) && !is.na(right_value)) {
        vec[i] <- (left_value + right_value) / 2

      }
      else{
        left_index <- max(1, left_index - 1)
        right_index <- min(length(vec), right_index + 1)
        left_value <- vec[left_index]
        right_value <- vec[right_index]
        if (!is.na(left_value) && !is.na(right_value)) {
          vec[i] <- (left_value + right_value) / 2
        }
        else{print(i)}
        
      }
    }
  }
  return(vec)
}