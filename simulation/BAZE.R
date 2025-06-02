# set to your own directory
setwd("/Users/kjiang/Desktop/BSRMM/simulation/")

# function
source("./function/gibbsgamma.R")
source("./function/BayesFactor.R")

# required library
library(truncnorm)

## ------------------------
## initial parameter
## ------------------------
jobid <- as.numeric(Sys.getenv("LSB_JOBINDEX"))
jobid 

# missing type
missing_type <- Sys.getenv("missing_type")
missing_type

# number of samples 
n <- as.numeric(Sys.getenv("n"))
n 

# number of features 
p <- as.numeric(Sys.getenv("p"))
p

# snr 
snr <- as.numeric(Sys.getenv("snr"))
snr 

# proportion of missing value
proportion <- as.numeric(Sys.getenv("proportion"))
proportion

# correlated data 
correlated <- Sys.getenv("correlated")
correlated

# imputation method
imputation_method <- Sys.getenv("imputation_method")
imputation_method

# # temporary parameter 
# imputation_method <- "low"
# missing_type <- "mnar0.33"
# correlated <- "cor"
# n <- 300
# p <- 1000
# snr <- 1
# proportion <- 0.40
# jobid <- 22

## ------------------------
## data/result path
## ------------------------
if (correlated == "cor") {
    mydata_path <- paste0("./raw_data/", imputation_method, "_", missing_type, "_cor_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
    result_path <- paste0("./output/baze_", imputation_method, "_imputation_", missing_type, "_cor_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
} else {
    mydata_path <- paste0("./raw_data/", imputation_method, "_", missing_type, "_ind_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
    result_path <- paste0("./output/baze_", imputation_method, "_imputation_", missing_type, "_ind_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
}

# import data
mydata <- readRDS(mydata_path)

# extract predictor and outcome 
predictor <- mydata$X
outcome <- mydata$Y_miss
dt <- data.frame(outcome)

# scale and center the predictor 
stand <- list()
stand$mux <- colMeans(predictor)
stand$Sx <- apply(predictor, 2, sd)
predictor <- apply(predictor, 2, function(x) (x - mean(x))/sd(x))

# scale and center the outcome 
stand$muy <- mean(outcome)
stand$Sy <- sd(outcome)
outcome <- (outcome - stand$muy)/stand$Sy

# create Ri 
dt$Ri <- 0
dt$Ri[which(mydata$Y_true == mydata$Y_miss)] <- 1

# split the data into train and test
set.seed(1996)
seed <- 1996

train_ind <- sample(nrow(predictor),size = round(0.7*nrow(predictor)),replace = FALSE)
train_ind <- sort(train_ind)
train_ind

predictor_train <- predictor[train_ind,]
predictor_test <- predictor[-train_ind,]

outcome_train <- outcome[train_ind,]
outcome_test <- outcome[-train_ind,]

dt_train <- dt[train_ind,]
dt_test <- dt[-train_ind,]

# calcaulte normalized root mean squared error for the outcome_train 
nrmse_dt <- data.frame(mydata$Y_true, mydata$Y_miss)
nrmse_dt <- nrmse_dt[train_ind,]
idx <- which(nrmse_dt$mydata.Y_true != nrmse_dt$mydata.Y_miss)
nrmse <- sqrt(mean((nrmse_dt$mydata.Y_true[idx] - nrmse_dt$mydata.Y_miss[idx])^2)/var(nrmse_dt$mydata.Y_true[idx]))
nrmse 

# create the data matrix for BAZE
X <- as.matrix(predictor_train)

Y <- as.matrix(outcome_train)

n <- nrow(X)
p <- ncol(X)

c <- 100
N <- rbind(diag(x = 1, nrow = p, ncol = p), matrix(data = c, nrow = 1, ncol = p))

nburnin <- 10000
niter <- 20000
nop <- floor(n/2)

# construct the Q matrix for both correlated and independent data 
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

# set hyperparameters 
tau <- 1
nu <- 0
omega <- 0

# set shrinkage parameter for Ising Prior 
a0=-12
a <- rep(a0,p)

# run the BAZE model
baze <- gibbsgamma(nburnin, niter, p, nop, Y, X, N, a, Q, n, tau, nu, omega, seed, TRUE, stand, TRUE)

# selected features 
PPI <- colMeans(baze$gamma[(nburnin+2):(nburnin+niter+1),]) 
beta_select <- rep(0,ncol(X))
beta_select[which(PPI > 0.5)] <- 1

# TPR, FPR, TNR, FNR and MSE
if (sum(PPI > 0.5) == 0) {
  TPR = 0
  FNR = 1
  TNR = 1
  FPR = 0
  TP = 0
  FN = 1
  TN = 1
  FP = 0

  # obtain y_true
  y_true <- mydata$Y_true
  y_true_test <- y_true[-train_ind]
  y_true_train <- y_true[train_ind]
  
  # obtain the predict y
  y_pred <- stand$muy
  
  # mse calculated by 
  observed_data_mse <- mean((dt_test$outcome[which(dt_test$Ri!=0)] - y_pred)^2)

  # obtain the F1 score
  F1_score <- TP/(TP + 0.5*(FP + FN))

  # l2 loss
  coef_l2 = sqrt(sum((mydata$coeffs[,2] - 0)^2))
  
  # collect result 
  result <-
    list(
      TPR = TPR,
      FPR = FPR,
      TNR = TNR,
      FNR = FNR,
      TP = TP,
      FN = FN,
      TN = TN,
      FP = FP,
      F1_score = F1_score,
      PPI = PPI,
      observed_data_mse = observed_data_mse,
      coef_l2 = coef_l2
    )
  
} else if (sum(PPI > 0.5) >= 1) {
  true_beta <- mydata$coeffs[, "gammatrue"]
  t <- table(beta_select, true_beta)
  t
  TPR <- t[2,2]/(sum(t[,2]))
  FNR <- 1- TPR
  TNR <- t[1,1]/(sum(t[,1]))
  FPR <- 1 - TNR

  TP <- t[2,2]

  FN <- t[1,2]

  TN <- t[1,1]

  FP <- t[2,1]
  
  # calculate the beta_hat 
  Xri <- X[, which(PPI > 0.5)]
  Tri <- N[, which(PPI > 0.5)]
  Ari <- t(Xri) %*% Xri + tau^(-2) * (t(Tri) %*% Tri)
  invAri<-chol2inv(chol(Ari))

  # beta_hat
  beta_hat <- invAri %*% t(Xri) %*% Y
  
  # scale predictor test 
  X_test <- predictor_test

  # obtain the predict y
  Y_pred <- stand$Sy * as.matrix(X_test[,which(PPI > 0.5)]) %*% beta_hat + stand$muy
  
  # mse calculated by 
  observed_data_mse <- mean((Y_pred[which(dt_test$Ri!= 0), ] - as.matrix(dt_test$outcome)[which(dt_test$Ri!=0), ])^2)

  # obtain the F1 score
  F1_score <- TP/(TP + 0.5*(FP + FN))
  
  # l2 loss
  beta_final <- rep(0, p)
  beta_final[which(PPI > 0.5)] <- stand$Sy * beta_hat * (1/stand$Sx[which(PPI > 0.5)])
  coef_l2 = sqrt(sum((mydata$coeffs[,2] - beta_final)^2))

  # collect result 
  result <-
    list(
      TPR = TPR,
      FPR = FPR,
      TNR = TNR,
      FNR = FNR,
      TP = TP,
      FN = FN,
      TN = TN,
      FP = FP,
      F1_score = F1_score,
      PPI = PPI,
      observed_data_mse = observed_data_mse,
      coef_l2 = coef_l2
    )
} 

result$nrmse <- nrmse

saveRDS(result, file = result_path)
