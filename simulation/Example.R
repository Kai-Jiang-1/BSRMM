rm(list =ls())
# Exmaple
# Simulation: 
setwd("/Users/kjiang/Desktop/BSRMM/simulation/")

# function
source("./function/gen_simdata_cor.R")
source("./function/gen_missing_value_cor.R")
source("./function/bsrmmgibbs.R")
source("./function/bsrmmbf.R")

# library 
library(truncnorm)
library(mvnfast)
library(MASS)
library(coda)
library(MCMCpack)

# parameter 
seed <- 22
missing_type <- "mnar0.33"
n <- 300
p <- 1000
snr <- 1
proportion <- 0.40
imputation_method <- "low" # our compared methods require low or mean imputation, so we generate the data with low or mean imputation. 

set.seed(seed)

mydata <-
  gen_missing_value_cor(
    proportion = proportion,
    n = n,
    p = p,
    snr = snr,
    type = missing_type
  )

# BSRMM 
# extract predictor and outcome 
predictor <- mydata$X
outcome <- mydata$Y_miss

# scale and center the predictor
stand <- list()
stand$mux <- colMeans(predictor)
stand$Sx <- apply(predictor, 2, sd)
predictor <- apply(predictor, 2, function(x) (x - mean(x))/sd(x))

# we already performed mean or low missing value imputation during that simulation generation
# we need to re-set the value into zero
index <- which(mydata$Y_miss != mydata$Y_true)
index <- sort(index)
outcome[index] <- 0 

# save dt for following prediction error analysis: 0 = miss, 1 = observed
dt <- data.frame(outcome)
dt$Ri <- 0
dt$Ri[which(dt$outcome != 0)] <- 1

# split the data into train and test
set.seed(123)
seed <- 123

train_ind <- sample(nrow(predictor),size = round(0.7*nrow(predictor)),replace = FALSE)
train_ind <- sort(train_ind)
train_ind

predictor_train <- predictor[train_ind,]
predictor_test <- predictor[-train_ind,]

outcome_train <- outcome[train_ind,]
outcome_test <- outcome[-train_ind,]

dt_train <- dt[train_ind,]
dt_test <- dt[-train_ind,]

X <- as.matrix(predictor_train)

Y <- as.matrix(outcome_train)

# preparation for gibbs sampling
n <- nrow(X)
p <- ncol(X)
c <- 100 # larger c to maintain the compositionality of the predictors 

# construct the T matrix, because we cannot assign a matrix to T due to R's limitation 
# we name it as N matrix 
N <- rbind(diag(x = 1, nrow = p, ncol = p), matrix(data = c, nrow = 1, ncol = p))

# number of sampling procedure 
nburnin <- 10000
niter <- 20000

# initial number of predictors 
nop <- floor(n/2)

# construct the Q matrix for both correlated and independent data 
correlated <- "cor" # "ind" for independent, "cor" for correlated
if (correlated == "ind") {
  Q <- diag(0,nrow = p) 
} else if (correlated == "cor") {
  Q <- 0.002 * (matrix(1, nrow = p, ncol = p) - diag(1, nrow = p))
  
  for(i in seq(180, 380, by = 20)) {
    for(j in seq(i + 20, 400, by = 20)) {
      Q[i, j] = 4
    }
  }
  
  for (i in seq(580, 780, by = 20)) {
    for (j in seq(i + 20, 800, by = 20)) {
      Q[i, j] = 4
    }
  }
  
  for (i in 445:459) {
    for (j in (i + 1):460) {
      Q[i, j] = 4
    }
  }
  
  for (i in 945:959) {
    for (j in (i + 1):960) {
      Q[i, j] = 4
    }
  }
  
  for (i in 45:59) {
    for (j in (i + 1):60) {
      Q[i, j] = 4
    }
  }
  
  Q = (Q + t(Q)) / 2
}

# initial hyperparameters 
tau <- 1 
nu <- 0
omega <- 0 

a0 <- -12
a <- rep(a0,p)

predict_result = TRUE
display = TRUE
bayestmetab = TRUE

bsrmm <- bsrmmgibbs(nburnin, niter, p, nop, Y, X, N, a, Q, n, tau, nu, omega, seed, predict_result, bayestmetab, display)

# # select features
PPI <- colMeans(bsrmm$gamma[(nburnin+2):(nburnin+niter+1),]) # posterior probabilities of inclusion
beta_select <- rep(0, ncol(X))
beta_select[which(PPI > 0.5)] <- 1
which(PPI > 0.5)

# save the result
if (sum(PPI > 0.5) == 0) {
  TPR = 0
  FNR = 1
  TNR = 1
  FPR = 0
  TP = 0
  FN = 1
  TN = 1
  FP = 0
  
  # obtain the missing value imputation 
  Y_mvi <- as.matrix(rowMeans(bsrmm$Y_mvi[, (nburnin+2):(nburnin+niter+1)]))
  
  Y_true <- as.matrix(mydata$Y_true)
  Y_true_train <- as.matrix(Y_true[train_ind, ])
  Y_true_test <- as.matrix(Y_true[-train_ind, ])
  
  # normalized root mean squared error 
  nrmse <- sqrt(mean((Y_mvi[which(bsrmm$Ri==0), ] - Y_true_train[which(bsrmm$Ri==0), ])^2)/var(Y_true_train[which(bsrmm$Ri==0), ]))
  nrmse
  
  # obtain the predict Y
  Y_pred <- mean(bsrmm$beta_0[(nburnin+2):(nburnin+niter+1)])
  
  # observed mse 
  observed_data_mse <- mean((Y_pred - dt_test$outcome[which(dt_test$Ri!=0)])^2)
  
  # obtain the F1 score
  F1_score <- TP/(TP + 0.5*(FP + FN))
  
  beta_0 <- mean(bsrmm$beta_0[(nburnin+2):(nburnin+niter+1)])
  # if no variable is selected, the beta will be zero, but we still have the intercept term 
  coef_l2 <- sqrt((0 - beta_0)^2 + sum((mydata$coeffs[, 2] - 0)^2))
  
  # store result
  result <- list(TPR = TPR, FNR = FNR, TNR = TNR, FPR = FPR, TP = TP, FN = FN, TN = TN, FP = FP, 
                 Y_mvi = Y_mvi, nrmse = nrmse, 
                 observed_data_mse = observed_data_mse, F1_score = F1_score, 
                 coef_l2 = coef_l2, alpha = bsrmm$alpha, 
                 beta_0 = mean(bsrmm$beta_0[(nburnin+2):(nburnin+niter+1)]), 
                 sigma2 = bsrmm$sigma2)
  
} else {
  # confusion matrix 
  true_beta <- mydata$coeffs[, "gammatrue"]
  t <- table(beta_select, true_beta)
  TPR <- t[2,2]/(sum(t[,2]))
  FNR <- 1- TPR
  TNR <- t[1,1]/(sum(t[,1]))
  FPR <- 1 - TNR
  
  TP <- t[2,2]
  
  FN <- t[1,2]
  
  TN <- t[1,1]
  
  FP <- t[2,1]
  
  # obtain the missing value imputation 
  Y_mvi <- as.matrix(rowMeans(bsrmm$Y_mvi[, (nburnin+2):(nburnin+niter+1)]))
  
  # recalculate beta
  Xri <- X[, which(PPI > 0.5)]
  Tri <- N[, which(PPI > 0.5)]
  Lambda <- diag(1, nrow = n, ncol = n) - matrix(1, nrow = n, ncol = n)/(n + 1)
  Ari <- t(Xri) %*% Lambda %*% Xri + tau^(-2) * t(Tri) %*% Tri
  invAri <- chol2inv(chol(Ari))
  # beta = Ari^(-1) X^T (I - J/(n+1)) Y
  beta <- invAri %*% t(Xri) %*% Lambda %*% Y_mvi
  
  Y_true <- as.matrix(mydata$Y_true)
  Y_true_train <- as.matrix(Y_true[train_ind, ])
  Y_true_test <- as.matrix(Y_true[-train_ind, ])
  
  # normalized root mean squared error
  nrmse <- sqrt(mean((Y_mvi[which(bsrmm$Ri==0), ] - Y_true_train[which(bsrmm$Ri==0), ])^2)/var(Y_true_train[which(bsrmm$Ri==0), ]))
  
  # obtain the predict Y
  X_test <- as.matrix(predictor_test)
  Y_pred <- X_test[, which(PPI > 0.5)] %*% beta + mean(bsrmm$beta_0[(nburnin+2):(nburnin+niter+1)])
  
  # observed mse
  observed_data_mse <- mean((Y_pred[which(dt_test$Ri!= 0), ] - as.matrix(dt_test$outcome)[which(dt_test$Ri!=0), ])^2)  
  
  # obtain the F1 score
  F1_score <- TP/(TP + 0.5*(FP + FN))
  
  beta_final <- rep(0, ncol(X))
  beta_final[which(PPI > 0.5)] <- beta/as.matrix(stand$Sx[which(PPI > 0.5)])
  
  # beta_0
  beta_0 <- mean(bsrmm$beta_0[(nburnin+2):(nburnin+niter+1)])
  coef_l2 <- sqrt((0 - beta_0)^2 + sum((mydata$coeffs[, 2] - beta_final)^2))
  
  # store result
  result <- list(TPR = TPR, FNR = FNR, TNR = TNR, FPR = FPR, TP = TP, FN = FN, TN = TN, FP = FP, 
                 Y_mvi = Y_mvi, nrmse = nrmse, 
                 observed_data_mse = observed_data_mse, F1_score = F1_score, 
                 coef_l2 = coef_l2, alpha = bsrmm$alpha, 
                 beta_0 = mean(bsrmm$beta_0[(nburnin+2):(nburnin+niter+1)]), 
                 sigma2 = bsrmm$sigma2)
  
}

