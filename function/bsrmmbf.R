bsrmmbf <- function(Y, Xri, Xi, XIi, Tri, Ti, TIi, invAri, n, tau, nu, omega, flag, keep) {
  keep <- keep
  # lri1 <- Rfast::cholesky(invAri, parallel = TRUE)
  # Assuming invAri is your matrix
  lri1 <- chol(invAri)
  # print(max(eigen(lri1)$values))
  Lri <- t(lri1)
  Lambda <- diag(1, nrow = n, ncol = n) - matrix(1, nrow = n, ncol = n)/(n + 1)
  sqrtdetinvAri <- sum(log(diag(Lri)))
  resAri <- t(Y) %*% Lambda %*% Y - t(Y) %*% Lambda %*% Xri %*% invAri %*% t(Xri) %*% Lambda %*% Y
  # if (resAri < 0 ) {
  #     resAri <- 1e-8
  #   }
  if (flag) {
    # copy invAi from previous step
    # keep1 <- keep
    sqrtdetinvAi <- keep$sqrtdetinvAi
    resAi <- keep$resAi
    
    keep$resAri <- resAri
    keep$sqrtdetinvAri <- sqrtdetinvAri
    
    logratiodetT <- sum(log(diag(chol(t(Tri) %*% Tri)))) - sum(log(diag(chol(t(TIi) %*% TIi))))
    
    # BF1 <- exp(-log(tau) + (sqrtdetinvAri - sqrtdetinvAi) + logratiodetT) * ((nu * omega + resAri) / (nu * omega + resAi))^((n + nu) / 2)
    # invAi <- invAi
    BF1 <- -log(tau) + (sqrtdetinvAri - sqrtdetinvAi) + logratiodetT + (n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi))
    BF1 <- exp(BF1)
    return(list(BF1 = BF1, keep = keep, lri1 = lri1, Lri = Lri, Lambda = Lambda, sqrtdetinvAi = sqrtdetinvAi, resAri = resAri, 
    resAi = resAi, logratiodetT = logratiodetT))

    # return(list(BF1 = BF1, keep = keep))
  } else {
    # update invAi from invAri
    # keep <- keep
    Srii <- t(Xri) %*% Lambda %*% Xi + tau^(-2) * (t(Tri) %*% Ti)
    sii <- t(Xi) %*% Lambda %*% Xi + tau^(-2) * (t(Ti) %*% Ti)
    
    # calculate inverse Ai^(-1) given Ari^(-1)
    v1 <- sqrt(1 / (sii * (1 - t(Srii) %*% invAri %*% Srii / sii)))
    v <- v1[1,1] * invAri %*% Srii
    #v <- sqrt(1 / (sii * (1 - t(Srii) %*% invAri %*% Srii / sii))) %*% invAri %*% Srii
    
    A11 <- invAri + v %*% t(v)
    A12 <- -A11 %*% (Srii / sii[1,1])
    A21 <- -(1 / sii) %*% t(Srii) %*% A11
    A22 <- 1 / sii[1,1] + (1 / sii[1,1]) %*% t(Srii) %*% A11 %*% Srii / sii[1,1]
    invAi <- rbind(cbind(A11, A12), cbind(A21, A22))
    resAi <- t(Y) %*% Lambda %*% Y - t(Y) %*% Lambda %*% XIi %*% invAi %*% t(XIi) %*% Lambda %*% Y
    
    # if (resAi < 0 ) {
    #   resAi <- 1e-8
    # }
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
    
    # BF <- exp(-log(tau) + (sqrtdetinvAri - sqrtdetinvAi) + logratiodetT) * ((nu * omega + resAri) / (nu * omega + resAi))^((n + nu) / 2)
    BF <- -log(tau) + (sqrtdetinvAri - sqrtdetinvAi) + logratiodetT + (n + nu) / 2 * log((nu * omega + resAri) / (nu * omega + resAi))
    BF <- exp(BF)
    return(list(BF = BF, keep = keep, invAi = invAi, lri1 = lri1, Lri = Lri, Lambda = Lambda, 
    sqrtdetinvAri = sqrtdetinvAri, resAri = resAri, Srii = Srii, v1 = v1, v = v, A11 = A11, 
    A12 = A12, A21 = A21, A22 = A22, resAi = resAi))
    # return(list(BF = BF, keep = keep, invAi = invAi))
  }
}
