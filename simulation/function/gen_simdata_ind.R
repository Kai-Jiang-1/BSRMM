gen_simdata_ind <- function(n, p, gamma, b, theta, sigmaX, sigma) {
  # Calculate beta as element-wise multiplication of gamma and b
  beta <- gamma * b
  # Generate matrix X of dimensions n x p, with N(0, sigmaX^2)
  X <- matrix(theta + sigmaX * rnorm(n*p),n, p)
  # Generate noise epsilon with N(0, sigma^2)
  epsilon <- sigma * rnorm(n)
  
  return(list(X = X, beta = beta, epsilon = epsilon))
}