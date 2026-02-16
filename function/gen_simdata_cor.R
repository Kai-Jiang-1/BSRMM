gen_simdata_cor <- function(n, gamma, b, theta, sigmaX, Xcor, sigma){
  beta<- gamma*b
  
  X <- MASS:: mvrnorm(n = n, mu = theta, Sigma = diag(sigmaX) %*% Xcor %*% diag(sigmaX))
  
  epsilon <- sigma * rnorm(n)
  
  return(list(X = X, beta = beta, epsilon = epsilon))
}