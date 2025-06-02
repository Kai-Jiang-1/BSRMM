BayesFactor <- function(Y, Xri, Xi, XIi, Tri, Ti, TIi, invAri, n, tau, nu, omega, flag, keep) {
  keep <- keep
  # lri1 <- Rfast::cholesky(invAri, parallel = TRUE)
  # Assuming invAri is your matrix
  lri1 <- chol(invAri)
  # print(max(eigen(lri1)$values))
  Lri <- t(lri1)
  
  sqrtdetinvAri <- sum(log(diag(Lri)))
  resAri <- t(Y) %*% Y - t(Y) %*% Xri %*% invAri %*% t(Xri) %*% Y

  if (flag) {
    # copy invAi from previous step
    # keep1 <- keep
    sqrtdetinvAi <- keep$sqrtdetinvAi
    resAi <- keep$resAi

    keep$resAri <- resAri
    keep$sqrtdetinvAri <- sqrtdetinvAri

    logratiodetT <- sum(log(diag(chol(t(Tri) %*% Tri)))) - sum(log(diag(chol(t(TIi) %*% TIi))))
    
    # cat(paste("log(tau) is ", log(tau), "\n"))
    # cat(paste("(sqrtdetinvAri - sqrtdetinvAi) is ", (sqrtdetinvAri - sqrtdetinvAi), "\n"))
    # cat(paste("logratiodetT is ", logratiodetT, "\n"))
    # cat(paste("(n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi)) is ", 
    #       (n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi)), "\n"))
    
    # invAi <- invAi
    F1 <- -log(tau) + (sqrtdetinvAri - sqrtdetinvAi) + logratiodetT + (n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi))
    F1 <- exp(F1)
    return(list(F1 = F1, keep = keep))
  } else {
    # update invAi from invAri
    # keep <- keep
    Srii <- t(Xri) %*% Xi + tau^(-2) * (t(Tri) %*% Ti)
    sii <- t(Xi) %*% Xi + tau^(-2) * (t(Ti) %*% Ti)

    # calculate inverse Ai^(-1) given Ari^(-1)
    v1 <- sqrt(1 / (sii * (1 - t(Srii) %*% invAri %*% Srii / sii)))
    v <- v1[1,1] * invAri %*% Srii
    #v <- sqrt(1 / (sii * (1 - t(Srii) %*% invAri %*% Srii / sii))) %*% invAri %*% Srii

    A11 <- invAri + v %*% t(v)
    A12 <- -A11 %*% (Srii / sii[1,1])
    A21 <- -(1 / sii) %*% t(Srii) %*% A11
    A22 <- 1 / sii[1,1] + (1 / sii[1,1]) %*% t(Srii) %*% A11 %*% Srii / sii[1,1]
    invAi <- rbind(cbind(A11, A12), cbind(A21, A22))
    resAi <- t(Y) %*% Y - t(Y) %*% XIi %*% invAi %*% t(XIi) %*% Y

    # calculate determinant of Ai^(-1) and Ari^(-1)
    # Rank 1 update to Cholesky factorization
    tilLri <- t(chol(Lri %*% t(Lri) + v %*% t(v)))
    # print(tilLri)
    Lrii <- forwardsolve(tilLri, A12)
    # Lrii <- chol2inv(tilLri, A12)
    lii <- sqrt(A22 - t(Lrii) %*% Lrii)
    Li <- rbind(cbind(tilLri, matrix(0, nrow(tilLri), 1)), cbind(t(Lrii), lii))
    # print(Li)
    sqrtdetinvAi <- sum(log(diag(Li)))

    # keep quantities for future steps
    keep$resAri <- resAri
    keep$sqrtdetinvAri <- sqrtdetinvAri

    keep$resAi <- resAi
    keep$sqrtdetinvAi <- sqrtdetinvAi
    logratiodetT <- -sum(log(diag(chol(t(Tri) %*% Tri)))) + sum(log(diag(chol(t(TIi) %*% TIi))))
    
    # cat(paste("log(tau) is ", log(tau), "\n"))
    # cat(paste("(sqrtdetinvAri - sqrtdetinvAi) is ", (sqrtdetinvAri - sqrtdetinvAi), "\n"))
    # cat(paste("logratiodetT is ", logratiodetT, "\n"))
    # cat(paste("(n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi)) is ", 
    #       (n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi)), "\n"))


    F <- -log(tau) + (sqrtdetinvAri - sqrtdetinvAi) + logratiodetT + (n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi))
    F <- exp(F)
    return(list(F = F, keep = keep, invAi = invAi))
  }

  # calculate the Bayes factor: bf = P(Y|gammai=0,gamma_-i) / P(Y|gammai=1,gamma_-i)
  # F <- -log(tau) + (sqrtdetinvAri - sqrtdetinvAi) + logratiodetT + (n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi))
  # F <- exp(F)

  # return(list(F = F, keep = keep, invAi = invAi))
}
