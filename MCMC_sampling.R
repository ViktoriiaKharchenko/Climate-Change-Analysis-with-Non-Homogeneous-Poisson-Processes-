# Log-likelihood function
likelihood_function <- function(alpha, sigma, t_i, T) {
  log_lambda_t_i <- log((alpha / sigma) * (t_i / sigma)^(alpha - 1))
  m_T <- (T / sigma)^alpha
  log_likelihood <- sum(log_lambda_t_i) - m_T
  return(log_likelihood)
}



mcmc_sampling <- function(t_i, T, n_iter, burn_in, thinning = 100) {
  
  alpha_prior <- function() runif(1, 0, 5)
  sigma_prior <- function() runif(1, 0, 10000)
  
  alpha <- alpha_prior()
  sigma <- sigma_prior()
  
  # Calculate the number of samples to keep after burn-in and thinning
  n_samples <- ceiling((n_iter - burn_in) / thinning)
  
  # Storage for thinned samples after burn-in
  samples <- matrix(nrow = n_samples, ncol = 2)
  colnames(samples) <- c("alpha", "sigma")
  
  sample_count <- 0
  stored_samples <- 0
  
  for (i in 1:n_iter) {
    # Metropolis step for alpha
    alpha_proposed <- alpha_prior()
    likelihood_current <- likelihood_function(alpha, sigma, t_i, T)
    likelihood_proposed <- likelihood_function(alpha_proposed, sigma, t_i, T)
    
    if (log(runif(1)) < (likelihood_proposed - likelihood_current)) {
      alpha <- alpha_proposed
    }
    
    
    # Metropolis step for sigma
    sigma_proposed <- sigma_prior()
    likelihood_current <- likelihood_function(alpha, sigma, t_i, T)
    likelihood_proposed <- likelihood_function(alpha, sigma_proposed, t_i, T)
    
    if (log(runif(1)) < (likelihood_proposed - likelihood_current)) {
      sigma <- sigma_proposed
    }
    
    
    sample_count <- sample_count + 1
    
    # Store thinned samples after burn-in
    if (i > burn_in && sample_count %% thinning == 0) {
      stored_samples <- stored_samples + 1
      samples[stored_samples, ] <- c(alpha, sigma)
    }
  }
  
  return(samples)
}


likelihood_function_with_changepoint <- function(alpha1, sigma1, alpha2, sigma2, t_i, T, tau) {
  # Points before the change-point
  t_i_before <- t_i[t_i < tau]
  log_lambda_t_i_before <- log((alpha1 / sigma1) * (t_i_before / sigma1)^(alpha1 - 1))
  m_tau <- (tau / sigma1)^alpha1
  
  # Points after the change-point
  t_i_after <- t_i[t_i >= tau]
  log_lambda_t_i_after <- log((alpha2 / sigma2) * (t_i_after / sigma2)^(alpha2 - 1))
  m_T_minus_m_tau <- (T / sigma2)^alpha2 - (tau / sigma2)^alpha2
  
  # Compute log-likelihood
  log_likelihood <- sum(log_lambda_t_i_before) - m_tau + sum(log_lambda_t_i_after) - m_T_minus_m_tau
  return(log_likelihood)
}

tau_prior <- function(T) runif(1, 1, 124)
alpha_prior <- function() rgamma(1, 1.3, 1)
sigma_prior <- function() runif(1, 1, 10)

mcmc_sampling_with_changepoint <- function(t_i, T, n_iter = 200000, burn_in = 11000, thinning = 100) {

  
  alpha1 <- alpha_prior()
  sigma1 <- sigma_prior()
  alpha2 <- alpha_prior()
  sigma2 <- sigma_prior()
  tau <- tau_prior(T)
  
  # Calculate the number of samples to keep after burn-in and thinning
  n_samples <- ceiling((n_iter - burn_in) / thinning)
  
  # Storage for samples
  samples <- matrix(nrow = n_samples, ncol = 5)
  colnames(samples) <- c("alpha1", "sigma1", "alpha2", "sigma2", "tau")
  
  sample_count <- 0
  stored_samples <- 0
  
  for (i in 1:n_iter) {
    
    # Metropolis step for alpha1
    alpha1_proposed <- alpha_prior()
    likelihood_current <- likelihood_function_with_changepoint(alpha1, sigma1, alpha2, sigma2, t_i, T, tau)
    likelihood_proposed <- likelihood_function_with_changepoint(alpha1_proposed, sigma1, alpha2, sigma2, t_i, T, tau)
    
    if (log(runif(1)) < (likelihood_proposed - likelihood_current)) {
      alpha1 <- alpha1_proposed
    }
    
    # Metropolis step for sigma1
    sigma1_proposed <- sigma_prior()
    likelihood_current <- likelihood_function_with_changepoint(alpha1, sigma1, alpha2, sigma2, t_i, T, tau)
    likelihood_proposed <- likelihood_function_with_changepoint(alpha1, sigma1_proposed, alpha2, sigma2, t_i, T, tau)
    
    if (log(runif(1)) < (likelihood_proposed - likelihood_current)) {
      
      sigma1 <- sigma1_proposed
    }
    # Metropolis step for alpha2
    alpha2_proposed <- alpha_prior()
    likelihood_current <- likelihood_function_with_changepoint(alpha1, sigma1, alpha2, sigma2, t_i, T, tau)
    likelihood_proposed <- likelihood_function_with_changepoint(alpha1, sigma1, alpha2_proposed, sigma2, t_i, T, tau)
    
    if (log(runif(1)) < (likelihood_proposed - likelihood_current)) {
      
      alpha2 <- alpha2_proposed
    }
    
    # Metropolis step for sigma2
    sigma2_proposed <- sigma_prior()
    likelihood_current <- likelihood_function_with_changepoint(alpha1, sigma1, alpha2, sigma2, t_i, T, tau)
    likelihood_proposed <- likelihood_function_with_changepoint(alpha1, sigma1, alpha2, sigma2_proposed, t_i, T, tau)
    
    if (log(runif(1)) < (likelihood_proposed - likelihood_current)) {
      
      sigma2 <- sigma2_proposed
    } 
    
    # Metropolis step for tau
    tau_proposed <- tau_prior(T)
    likelihood_current <- likelihood_function_with_changepoint(alpha1, sigma1, alpha2, sigma2, t_i, T, tau)
    likelihood_proposed <- likelihood_function_with_changepoint(alpha1, sigma1, alpha2, sigma2, t_i, T, tau_proposed)
    
    if (log(runif(1)) < (likelihood_proposed - likelihood_current)) {
      tau <- tau_proposed
    }
    
    sample_count <- sample_count + 1
    
    # Store thinned samples after burn-in
    if (i > burn_in && sample_count %% thinning == 0) {
      stored_samples <- stored_samples + 1
      samples[stored_samples, ] <- c(alpha1, sigma1, alpha2, sigma2, tau)
    }
    
  }
  
  return(samples)
}
