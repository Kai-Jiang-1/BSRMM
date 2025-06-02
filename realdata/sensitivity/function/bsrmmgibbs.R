bsrmmgibbs <- function(nburnin, niter, p, nop, Y, X, N, a, Q, n, tau, nu, omega, seed, predict_result, bayestmetab, display) {
  # Setting the seed
  set.seed(seed)
  
  # Intitalize section 
  # Initialize index, nop is the number of initial 1's
  index <- sample(1:p, nop, replace = FALSE)
  index <- sort(index)
  nop <- length(index)
  index
  
  # Initialize gamma - variable selection
  gamma <- matrix(0, nrow = nburnin + niter + 1, ncol = p) 
  gamma[1, index] <- 1
  
  # Initialized sigma2 - variable selection 
  sigma2 <- numeric(nburnin + niter + 1)

  # Initialize alpha  - missing value imputation 
  alpha <- numeric(nburnin + niter + 1)
  alpha[1] <- runif(1, min = 0, max = 1)
  
  
  # Update section for beta, sigma2, beta_0 
  # Ari = X^T (I - J/(n+1)) X + sigma^(-2) tau^(-2) T^T T
  # Let Lambda = (I - J/(n+1))
  Xri <- X[, index] 
  Tri <- N[, index] 
  Lambda <- diag(1, nrow = n, ncol = n) - matrix(1, nrow = n, ncol = n)/(n + 1)
  Ari <- t(Xri) %*% Lambda %*% Xri + tau^(-2) * (t(Tri) %*% Tri)
  
  # Ari^-1 
  invAri<-chol2inv(chol(Ari))
  
  # by using cholesky decomposition, we can obtain the determinant of the matrix 
  # the determinant of the matrix is the product of the diagonal elements of triangular matrix
  Lri <- t(chol(invAri))  
  
  # C_gamma <- Y^T (I - J/(n+1)) Y - Y^T (I - J/(n+1)) X Ari^-1 X^T (I - J/(n+1)) Y
  # resAi is the C_gamma in the supplemental document
  keep <- list()
  keep$resAi <- t(Y) %*% Lambda %*% Y - t(Y) %*% Lambda %*% X[, index] %*% invAri %*% t(X[, index]) %*% Lambda %*% Y 
  keep$sqrtdetinvAi <- sum(log(diag(Lri))) 
  
  # beta follows multivariate t-distribution with mean as Ari^-1 X^T (I - J/(n+1)) Y and 
  # scale matrix as (keep$resAi + nu * omega)/(n + nu) * Ari^-1
  tembeta <- as.matrix(t(rmvt(n = 1, mu = invAri %*% t(X[, index]) %*% Lambda %*% Y, 
                              sigma = (as.numeric(keep$resAi) + nu * omega)/(n + nu) * invAri, df = n + nu)))
  
  # save the result for tembeta
  tembeta_save <- matrix(NA, nrow = ncol(X), ncol = nburnin + niter + 1)
  tembeta_save[index, 1] <- tembeta
  
  # sigma2 is a inverse gamma distribution with scale parameter as 0.5 * C_gamma + nu * omega 
  # and shape parameter as 0.5 * (n + nu)
  # sigma2 = (0.5*C_gamma + nu*omega)/(0.5*(n+nu)-1)
  sigma2[1] <- rinvgamma(n = 1, shape = 0.5*(n+nu), scale = 0.5*(as.numeric(keep$resAi)+nu*omega))

  # beta_0 ~ N(sum(y-xb)/(n + 1), sigma^2/(n+1))
  beta_0 <- numeric(nburnin + niter + 1)
  beta_0[1] <- rnorm(1, mean = sum(Y - X[, index] %*% tembeta)/(n + 1), sd = sqrt(sigma2[1]/(n + 1)))
  
  # Update section for missing value imputation 
  # Missingness indicator: 0=Miss 1=Observed
  Ri <- rep(0,nrow(Y))
  Ri[Y != 0 ] <- 1
  
  # Limit of detection (LOD)
  xi <- sort(unique(Y))[1]
  
  if (bayestmetab == TRUE) {
    # probability of lower than lod 
    phi <- matrix(NA,nrow = n, ncol = nburnin + niter + 1)
    phi[which(Ri == 0), 1] <-
      pnorm(
        xi,
        mean = beta_0[1] + X[which(Ri == 0), index] %*% tembeta,
        sd = rep(sqrt(sigma2[1]), sum(Ri == 0)),
        lower.tail = T
      )
    
    # whether latent value is greater than lod or not, 
    # if z = 1, the latent true value is lower than lod 
    # if z = 0, the latent true value is greater than lod
    # we use rbinom to sample from the probability of lower than lod 
    phi_values <- phi[which(Ri == 0), 1]
    prob_values <- 1 / (1 + (alpha[1] * (1 - phi_values) / phi_values ))
    z <- matrix(data = NA,nrow = n, ncol = nburnin + niter + 1)
    z[which(Ri == 0), 1] <- rbinom(length(prob_values), size = 1, prob_values)
                            
    # missing value index 
    missing_index <- matrix(NA, nrow = n, ncol = nburnin + niter + 1) 
    missing_index[, 1] <- sapply(1:n, function(j) {
      if (Ri[j] == 1) {
        return(2) # observed
      } else if (Ri[j] == 0 & z[j, 1] == 1) {
        return(1) # lower than lod
      } else {
        return(0) # greater than lod
      }
    })
    
    # missing value imputation 
    y_mvi <- matrix(NA, nrow = n, ncol = nburnin + niter + 1)
    y_mvi[, 1] <- sapply(1:n, function(s) {
      if (missing_index[s, 1] == 1) {
        return(
          rtruncnorm(
            1,
            a = -Inf,
            b = xi,
            mean = beta_0[1] + X[s, index] %*% tembeta,
            sd = sqrt(sigma2[1])
          ) 
        )
      } else if (missing_index[s, 1] == 0) {
        return(
          rtruncnorm(
            1,
            a = xi,
            b = Inf,
            mean = beta_0[1] + X[s, index] %*% tembeta,
            sd = sqrt(sigma2[1])
          ) 
        )
      } else {
        return(Y[s, 1])
      }
    })
    Y <- as.matrix(y_mvi[,1])
    alpha_a <- sum(missing_index[, 1] == 0) + 1
    alpha_b <- sum(missing_index[, 1] == 2) + 1
    alpha[1] <-  rbeta(1, alpha_a, alpha_b)
  }
  
  nselect <- numeric(nburnin + niter + 1)
  nselect[1] <- nop
  
  if (display) {
    k <- 1
    cat("Gibbs Sampling is starting and will print every 5000 iterations.\n")
  }
  
  for(i in 1:(nburnin+niter)){
    proposeindx <- sample(1:p, 1, replace = TRUE)
    gamma[i+1,] <- gamma[i,]
    
    flag <- any(index == proposeindx) & (length(index) > 1)

    
    if (flag) {
      # proposed variable is inside existing index and there are more than 1 variable in the model
      # because the propose variable is inside the model
      # We consider if the variable should be removed or not
      indxtemp <- index[index != proposeindx] # indxtemp does not include proposeindx
      Xri <- X[, indxtemp]
      Xi <- X[, proposeindx]
      XIi <- cbind(Xri, Xi)
      Tri <- N[, indxtemp]
      Ti <- N[, proposeindx]
      TIi <- cbind(Tri, Ti)
      invAritemp <- invAri
      idx <- list()
      tn <- length(index)
      seq <- 1:tn
      if (proposeindx == max(index)) {
        idx[[1]] <- seq
        idx[[2]] <- seq
      } else if (proposeindx == min(index)) {
        idx[[1]] <- c(2:tn, 1)
        idx[[2]] <- c(2:tn, 1)
      } else {
        ti <- seq[index == proposeindx]
        idx[[1]] <- c(1:(ti - 1), (ti+1):tn, ti)
        idx[[2]] <- c(1:(ti - 1), (ti+1):tn, ti)
      }
      
      invAitemp <- invAri[idx[[1]], idx[[2]]]
      # swap the ti variable to the last one of the list
      invAri1 <- invAitemp[1:(tn-1), 1:(tn-1)]
      
      # print(dim(invAri1)[1])
      invAri2 <- invAitemp[1:(tn-1), tn, drop = FALSE]
      invAri3 <- invAitemp[tn, tn]
      invAri <- invAri1 - invAri2 %*% t(invAri2) / invAri3
      # invAri_1 <- invAritemp[1:(tn-1), 1:(tn-1)] - invAritemp[1:(tn-1), tn] %*% invAritemp[tn, 1:(tn-1)] / invAritemp[tn, tn]
      result1 <- bsrmmbf(Y, Xri, Xi, XIi, Tri, Ti, TIi, invAri, n, tau, nu, omega, flag, keep)
      BF1 <- result1$BF1
      keep <- result1$keep
      pgammai1 <- exp(a[proposeindx] + Q[proposeindx, indxtemp] %*% gamma[i, indxtemp])
      pcond <- 1 / (1 + 1 / (BF1 * pgammai1))
      newgamma <- rbinom(1, 1, pcond)
      gamma[i+1, proposeindx] <- newgamma
      if ((newgamma == 0) && (tn > 1)) {
        index <- indxtemp # propose variable is removed
        keep$resAi <- keep$resAri
        keep$sqrtdetinvAi <- keep$sqrtdetinvAri
        # beta = A^-1 X^T (I - J/(n+1)) Y
        tembeta <- as.matrix(t(rmvt(n = 1, mu = invAri %*% t(X[, index]) %*% Lambda %*% Y, , 
                                    sigma = (as.numeric(keep$resAi) + nu * omega)/(n + nu) * invAri, df = n + nu)))
        
        tembeta_save[index, i + 1] <- tembeta
        
        # sigma2 = (0.5*C_gamma + nu*omega)/(0.5*(n+nu)-1)
        sigma2[i + 1] <- rinvgamma(n = 1, shape = 0.5*(n+nu), scale = 0.5*(as.numeric(keep$resAi)+nu*omega))
        
        # beta_0 ~ N(sum(y-xb)/(n + 1), sigma^2/(n+1))
        beta_0[i + 1] <- rnorm(1, mean = sum(Y - X[, index] %*% tembeta)/(n + 1), sd = sqrt(sigma2[i+1]/(n + 1)))
        
        nselect[i + 1] <- tn - 1
      }else if ((newgamma != 0) || (tn == 1)) {
        # propose variable is included 
        invAri <- invAritemp
        nselect[i + 1] <- nselect[i]
        
        # sampling beta
        tembeta <- as.matrix(t(rmvt(n = 1, mu = invAri %*% t(X[, index]) %*% Lambda %*% Y, 
                                    sigma = (as.numeric(keep$resAi) + nu * omega)/(n + nu) * invAri, df = n + nu)))
        tembeta_save[index, i + 1] <- tembeta
        
        # sampling sigma2
        sigma2[i + 1] <- rinvgamma(n = 1, shape = 0.5*(n+nu), scale = 0.5*(as.numeric(keep$resAi)+nu*omega))
        
        # sampling beta_0 
        beta_0[i+1] <- rnorm(1, mean = sum(Y - X[, index] %*% tembeta)/(n + 1), sd = sqrt(sigma2[i+1]/(n + 1)))
      }
    }else if (any(index != proposeindx)) {
      indxtemp <- index # indxtemp not include proposeindx  
      Xri <- X[, indxtemp]
      Xi <- X[, proposeindx]
      XIi <- cbind(Xri, Xi)
      Tri <- N[, indxtemp]
      Ti <- N[, proposeindx]
      TIi <- cbind(Tri, Ti)
      BayesResult <- bsrmmbf(Y, Xri, Xi, XIi, Tri, Ti, TIi, invAri, n, tau, nu, omega, flag, keep)
      BF <- BayesResult$BF
      keep <- BayesResult$keep
      invAi <- BayesResult$invAi
      pgammai1 <- exp(a[proposeindx] + Q[proposeindx, indxtemp] %*% gamma[i, indxtemp])
      pcond <- 1 / (1 + 1/(BF * pgammai1))
      newgamma <- rbinom(1, 1, pcond)
      gamma[i + 1, proposeindx] <- newgamma
      if (newgamma == 1) {
        #include the proposed variable
        index <- sort(c(indxtemp, proposeindx)) # include proposeindx 
        invAri <- invAi
        idx <- list()
        tn <- length(index)
        seq <- 1:tn
        if (proposeindx > max(indxtemp)) {
          idx[[1]] <- seq
          idx[[2]] <- seq
        } else if (proposeindx < min(indxtemp)) {
          idx[[1]] <- c(tn, 1:(tn - 1))
          idx[[2]] <- c(tn, 1:(tn - 1))
        } else {
          ti <- seq[index == proposeindx]
          idx[[1]] <- c(1:(ti - 1), tn, ti:(tn - 1))
          idx[[2]] <- c(1:(ti - 1), tn, ti:(tn - 1))
        }
        invAri <- invAri[idx[[1]], idx[[2]]]
        if (predict_result) {
          # sample beta 
          tembeta <- as.matrix(t(rmvt(n = 1, mu = invAri %*% t(X[, index]) %*% Lambda %*% Y, 
                                      sigma = (as.numeric(keep$resAi) + nu * omega)/(n + nu) * invAri, df = n + nu)))
          tembeta_save[index, i + 1] <- tembeta
          
          # sample sigma2
          # sigma2 = (0.5*C_gamma + nu*omega)/(0.5*(n+nu)-1)
          sigma2[i + 1] <- rinvgamma(n = 1, shape = 0.5*(n+nu), scale = 0.5*(as.numeric(keep$resAi)+nu*omega))
          
          # sample beta_0
          beta_0[i+1] <- rnorm(1, mean = sum(Y - X[, index] %*% tembeta)/(n + 1), sd = sqrt(sigma2[i+1]/(n + 1)))
          
          # update nselect
          nselect[i + 1] <- tn
        }
      } else {
        keep$resAi <- keep$resAri
        keep$sqrtdetinvAi <- keep$sqrtdetinvAri
        if (predict_result) {
          # sample beta
          tembeta <- as.matrix(t(rmvt(n = 1, mu = invAri %*% t(X[, index]) %*% Lambda %*% Y, 
                                      sigma = (as.numeric(keep$resAi) + nu * omega)/(n + nu) * invAri, df = n + nu)))
          tembeta_save[index, i + 1] <- tembeta
          
          # sample sigma2
          sigma2[i + 1] <- rinvgamma(n = 1, shape = 0.5*(n+nu), scale = 0.5*(as.numeric(keep$resAi)+nu*omega))
          
          # sample beta_0
          beta_0[i+1] <- rnorm(1, mean = sum(Y - X[, index] %*% tembeta)/(n + 1), sd = sqrt(sigma2[i+1]/(n + 1)))
          
          # update nselect
          nselect[i + 1] <- nselect[i]
        }
      }
    } else {
      keep$resAi <- keep$resAri
      keep$sqrtdetinvAi <- keep$sqrtdetinvAri
      if (predict_result) {
        # sample beta
        tembeta <- as.matrix(t(rmvt(n = 1, mu = invAri %*% t(X[, index]) %*% Lambda %*% Y, 
                                    sigma = (as.numeric(keep$resAi) + nu * omega)/(n + nu) * invAri, df = n + nu)))
        tembeta_save[index, i + 1] <- tembeta
        
        # sample sigma2
        sigma2[i + 1] <- rinvgamma(n = 1, shape = 0.5*(n+nu), scale = 0.5*(as.numeric(keep$resAi)+nu*omega))
        
        # sample beta_0
        beta_0[i+1] <- rnorm(1, mean = sum(Y - X[, index] %*% tembeta)/(n + 1), sd = sqrt(sigma2[i+1]/(n + 1)))
        
        # update nselect
        nselect[i + 1] <- nselect[i]
      }
    }
    
    # Missing value imputation
    if (bayestmetab == TRUE) {
      # cdf for less than xi
      phi[which(Ri == 0), i + 1] <-
        pnorm(
          xi,
          mean = beta_0[i + 1] + X[which(Ri == 0), index] %*% tembeta,
          sd = rep(sqrt(sigma2[i + 1]), sum(Ri == 0)),
          lower.tail = T
        )
      
      # sample from rbinorm to decide if latent value is greater than lod or not.
      # if z = 1, the latent true value is lower than lod
      # if z = 0, the latent true value is greater than lod
      phi_values <- phi[which(Ri == 0), i + 1]
      prob_values <- 1 / (1 + (alpha[i] * (1 - phi_values) / phi_values))
      z[which(Ri == 0), i + 1] <- rbinom(length(prob_values), size = 1, prob_values)
    
      # given the z, we can decide the missing value index
      # 0=MAR 1=MNAR 2=Observed
      missing_index[, i + 1] <- sapply(1:n, function(j) {
        if (Ri[j] == 1) {
          return(2) # observed
        } else if (Ri[j] == 0 & z[j, i + 1] == 1) {
          return(1) # lower than lod
        } else {
          return(0) # greater than lod
        }
      })
      
      # Missing value imputation for iteration i + 1
      y_mvi[, i + 1] <- sapply(1:n, function(s) {
        if (missing_index[s, i + 1] == 1) {
          return(
            rtruncnorm(
              1,
              a = -Inf,
              b = xi,
              mean = beta_0[i + 1] + X[s, index] %*% tembeta,
              sd = sqrt(sigma2[i + 1])
            ) 
            # * stand$Sy[[i]] + stand$muy[[i]] # less than lod
          )
        } else if (missing_index[s, i + 1] == 0) {
          return(
            rtruncnorm(
              1,
              a = xi,
              b = Inf,
              mean = beta_0[i + 1] + X[s, index] %*% tembeta,
              sd = sqrt(sigma2[i + 1])
            ) 
            # * stand$Sy[[i]] + stand$muy[[i]] # greater than lod
          )
        } else {
          return(y_mvi[s, i])
        }
      })
      
      Y <- as.matrix(y_mvi[, i + 1])
      # Ri = 0 = missing
      # Ri = 1 = observed
      alpha_a <- sum(missing_index[, i + 1] == 0) + 1
      alpha_b <- sum(missing_index[, i + 1] == 2) + 1
      alpha[i + 1] <-  rbeta(1, alpha_a, alpha_b)
    }
    
    if (display) {
      if (k %% 5000 == 0) {
        print(paste("iteration is ", k))
      }
      k <- k + 1
      # print(k)
    }

  }
  if (display) {
    print("Gibbs Sampling is ending and will print the frequency of select variables.")
  }
  return(
    list(
      gamma = gamma, # For PPI 
      nselect = nselect, # For number of selected variables
      alpha = alpha, # posterior of alpha 
      phi = phi, # probability of lower than lod
      z = z, # whether latent value is greater than lod or not
      missing_index = missing_index, # missing value index
      Y_mvi = y_mvi, # missing value imputation
      Ri = Ri, # missingness indicator
      stand = stand, # standardization
      sigma2 = sigma2, # sigma2
      beta_0 = beta_0, # beta_0,
      tembeta_save = tembeta_save # save the result for tembeta
    )
  )
}
