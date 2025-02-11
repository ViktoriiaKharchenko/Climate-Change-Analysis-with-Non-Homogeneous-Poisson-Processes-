lambda <- function(t, alpha, sigma) {
  (alpha / sigma) * (t / sigma)^(alpha - 1)
}

m_t <- function(t, alpha, sigma) {
  (t / sigma)^alpha
}

# Define the intensity function 
lambda2 <- function(t, theta1, theta2, tau) {
  ifelse(t < tau, 
         (theta1[1] / theta1[2]) * (t / theta1[2])^(theta1[1] - 1), 
         (theta2[1] / theta2[2]) * (t / theta2[2])^(theta2[1] - 1))
}

# Define the mean value function
m_t2 <- function(t, theta1, theta2, tau) {
  ifelse(t < tau,
         (t / theta1[2])^theta1[1],
         (tau / theta1[2])^theta1[1] + ((t / theta2[2])^theta2[1] - (tau / theta2[2])^theta2[1]))
}