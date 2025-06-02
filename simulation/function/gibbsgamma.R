gibbsgamma <- function(nburnin, niter, p, nop, Y, X, N, a, Q, n, tau, nu, omega, seed, predict, stand, display) {

  # Setting the seed
  set.seed(seed)

  # Initialize index, nop is the number of initial 1's
  index <- sample(1:p, nop, replace = FALSE)
  index <- sort(index)
  nop <- length(index)

  index
  
  # Initialize gamma
  gamma <- matrix(0, nrow = nburnin + niter + 1, ncol = p)


  gamma[1, index] <- 1


  if (is.null(stand)) {
    invSx <- diag(p)
    Sy <- 1
    Yobs <- Y
  } else {
    invSx <- diag(1 / stand$Sx)
    Sy <- stand$Sy
    Yobs <- Y * stand$Sy + stand$muy
    # Xobs <- X * diag(stand$Sx) + stand$mux
  }

  proposeindx <- sample(1:p, nburnin + niter + 1, replace = TRUE)

  # Initialize Ari and keep
  Xri <- X[, index]
  Tri <- N[, index]
  Ari <- t(Xri) %*% Xri + tau^(-2) * (t(Tri) %*% Tri)
  invAri<-chol2inv(chol(Ari))
  Lri <- t(chol(invAri))

  keep <- list()
  keep$resAi <- t(Y) %*% Y - t(Y) %*% Xri %*% invAri %*% t(Xri) %*% Y
  keep$sqrtdetinvAi <- sum(log(diag(Lri)))

  if (predict) {
    # initialize beta
    betahat <- matrix(0, nrow = p, ncol = nburnin + niter + 1)
    tembeta <- invAri %*% t(Xri) %*% Y
    betahat[index, 1] <- Sy * invSx[index, index] %*% tembeta
    Yhat <- matrix(0, nrow = n, ncol = nburnin + niter + 1)
    Yhat[, 1] <- stand$muy + Sy * X[, index] %*% tembeta
    MSE <- numeric(nburnin + niter + 1)
    MSE[1] <- 1 / n * t(Yobs - Yhat[, 1]) %*% (Yobs - Yhat[, 1])
    nselect <- numeric(nburnin + niter + 1)
    nselect[1] <- nop
  }


  if (display) {
    k <- 1
    cat("Gibbs Sampling is starting and will print every 5000 iterations.\n")
  }

  for(i in 1:(nburnin+niter)){
    gamma[i+1,] <- gamma[i,]
    betahat[,i+1] <- rep(0, p)
    flag <- any(index == proposeindx[i]) & (length(index) > 1)
    if (flag) {
      indxtemp <- index[index != proposeindx[i]]
      Xri <- X[, indxtemp]
      Xi <- X[, proposeindx[i]]
      XIi <- cbind(Xri, Xi)
      Tri <- N[, indxtemp]
      Ti <- N[, proposeindx[i]]
      TIi <- cbind(Tri, Ti)
      invAritemp <- invAri
      idx <- list()
      tn <- length(index)
      seq <- 1:tn
      if (proposeindx[i] == max(index)) {
        idx[[1]] <- seq
        idx[[2]] <- seq
      } else if (proposeindx[i] == min(index)) {
        idx[[1]] <- c(2:tn, 1)
        idx[[2]] <- c(2:tn, 1)
      } else {
        ti <- seq[index == proposeindx[i]]
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
      result1 <- BayesFactor(Y, Xri, Xi, XIi, Tri, Ti, TIi, invAri, n, tau, nu, omega, flag, keep)
      F1 <- result1$F1
      keep <- result1$keep
      pgammai1 <- exp(a[proposeindx[i]] + Q[proposeindx[i], indxtemp] %*% gamma[i, indxtemp])
      pcond <- 1 / (1 + 1 / (F1 * pgammai1))
      # print(pcond)
      newgamma <- rbinom(1, 1, pcond)
      gamma[i+1, proposeindx[i]] <- newgamma
      if ((newgamma == 0) && (tn > 1)) {
        index <- indxtemp
        keep$resAi <- keep$resAri
        keep$sqrtdetinvAi <- keep$sqrtdetinvAri
        if (predict) {
          tembeta <- invAri %*% t(Xri) %*% Y
          betahat[index, i + 1] <- Sy * invSx[index, index] %*% tembeta
          Yhat[, i + 1] <- stand$muy + Sy * X[, index] %*% tembeta
          MSE[i + 1] <- 1 / n * t(Yobs - Yhat[, i + 1]) %*% (Yobs - Yhat[, i + 1])
          nselect[i + 1] <- tn - 1
        }

      }else if ((newgamma != 0) || (tn == 1)) {
        invAri <- invAritemp
        if (predict) {
          betahat[, i + 1] <- betahat[, i]
          Yhat[, i + 1] <- Yhat[, i]
          MSE[i + 1] <- MSE[i]
          nselect[i + 1] <- nselect[i]
        }
      }


    }else if (any(index != proposeindx[i])) {
      indxtemp <- index
      Xri <- X[, indxtemp]
      Xi <- X[, proposeindx[i]]
      XIi <- cbind(Xri, Xi)
      Tri <- N[, indxtemp]
      Ti <- N[, proposeindx[i]]
      TIi <- cbind(Tri, Ti)
      BayesResult <- BayesFactor(Y, Xri, Xi, XIi, Tri, Ti, TIi, invAri, n, tau, nu, omega, flag, keep)
      F <- BayesResult$F
      keep <- BayesResult$keep
      invAi <- BayesResult$invAi
      pgammai1 <- exp(a[proposeindx[i]] + Q[proposeindx[i], indxtemp] %*% gamma[i, indxtemp])
      pcond <- 1 / (1 + 1/(F * pgammai1))
      newgamma <- rbinom(1, 1, pcond)
      gamma[i + 1, proposeindx[i]] <- newgamma
      if (newgamma == 1) {
        index <- sort(c(indxtemp, proposeindx[i]))
        invAri <- invAi
        idx <- list()
        tn <- length(index)
        seq <- 1:tn
        if (proposeindx[i] > max(indxtemp)) {
          idx[[1]] <- seq
          idx[[2]] <- seq
        } else if (proposeindx[i] < min(indxtemp)) {
          idx[[1]] <- c(tn, 1:(tn - 1))
          idx[[2]] <- c(tn, 1:(tn - 1))
        } else {
          ti <- seq[index == proposeindx[i]]
          idx[[1]] <- c(1:(ti - 1), tn, ti:(tn - 1))
          idx[[2]] <- c(1:(ti - 1), tn, ti:(tn - 1))
        }
        invAri <- invAri[idx[[1]], idx[[2]]]
        if (predict) {
          tembeta <- invAri %*% t(X[, index]) %*% Y
          betahat[index, i + 1] <- Sy * invSx[index, index] %*% tembeta
          Yhat[, i + 1] <- stand$muy + Sy * X[, index] %*% tembeta
          MSE[i + 1] <- 1 / n * t(Yobs - Yhat[, i + 1]) %*% (Yobs - Yhat[, i + 1])
          nselect[i + 1] <- tn
        }
      }else {
        keep$resAi <- keep$resAri
        keep$sqrtdetinvAi <- keep$sqrtdetinvAri

        if (predict) {
          betahat[, i+1] <- betahat[, i]
          Yhat[, i+1] <- Yhat[, i]
          MSE[i+1] <- MSE[i]
          nselect[i+1] <- nselect[i]
        }
      }
    } else {
      keep$resAi <- keep$resAri
      keep$sqrtdetinvAi <- keep$sqrtdetinvAri
      
      if (predict) {
        betahat[, i + 1] <- betahat[, i]
        Yhat[, i + 1] <- Yhat[, i]
        MSE[i + 1] <- MSE[i]
        nselect[i + 1] <- nselect[i]
      }
    }
    
    if (display) {
      if (k %% 500 == 0) {
        print(paste("iteration is ", k))
      }
      k <- k + 1
    }



  }
  if (display) {
    print("Gibbs Sampling is ending and will print the frequency of select variables.")
    gammaresult <- colSums(gamma[nburnin:(nburnin + niter), ])
  }
  return(list(gamma = gamma, betahat = betahat, gammaresult = gammaresult, MSE = MSE, nselect = nselect, Yhat=Yhat))
}


