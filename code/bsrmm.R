# BRSMM (Bayesian sparse regression for microbiome-metabolite data integration)
# Kai Jiang

# set directory 
rm(list = ls())
setwd("/Users/kjiang/Desktop/BSRMM/")

# functions 
source("./function/bsrmmgibbs.R")
source("./function/bsrmmbf.R")

# required library
library(truncnorm)
library(mvnfast)
library(MASS)
library(coda)
library(MCMCpack)
library(reshape2)

# parameter
missing_type <- "mnar0.33"
correlated <- "ind"
n <- 300
p <- 1000
snr <- 1
proportion <- 0.30
imputation_method <- "low"
jobid <- 1


## ------------------------
## data/result path
## ------------------------
if (correlated == "cor") {
  mydata_path <- paste0("./data/", imputation_method, "_", missing_type, "_cor_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
  result_path <- paste0("./result/bsrmm_", imputation_method, "_imputation_", missing_type, "_cor_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
} else {
  mydata_path <- paste0("./data/", imputation_method, "_", missing_type, "_ind_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
  result_path <- paste0("./result/bsrmm_", imputation_method, "_imputation_", missing_type, "_ind_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
}

mydata_path
result_path

# import data
mydata <- readRDS(mydata_path)

# extract predictor and outcome 
predictor <- mydata$X
outcome <- mydata$Y_miss

# consider scale and center the predictor before split
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
set.seed(33)
seed <- 33

train_ind <- sample(nrow(predictor),size = round(0.7*nrow(predictor)),replace = FALSE)
train_ind <- sort(train_ind)
# train_ind

predictor_train <- predictor[train_ind,]
predictor_test <- predictor[-train_ind,]

outcome_train <- outcome[train_ind,]
outcome_test <- outcome[-train_ind,]

dt_train <- dt[train_ind,]
dt_test <- dt[-train_ind,]

X <- as.matrix(predictor_train)

# check the outcome 
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

quantile(Q)
sum(Q)

# initial hyperparameters 
tau <- 1 
nu <- 0
omega <- 0 

a0 <- -12
a <- rep(a0,p)

predict_result = TRUE
display = FALSE
bayestmetab = TRUE

# save the running time: 
bsrmm <- bsrmmgibbs(nburnin, niter, p, nop, Y, X, N, a, Q, n, tau, nu, omega, 
                    alpha_a_shape = 1, alpha_b_shape = 1, seed, predict_result, bayestmetab, display, lod_ratio = 1)

# Result
PPI <- colMeans(bsrmm$gamma[(nburnin+2):(nburnin+niter+1),]) # posterior probabilities of inclusion
beta_select <- rep(0, ncol(X))
beta_select[which(PPI > 0.5)] <- 1
which(PPI > 0.5)

# check the beta estiamtion 
beta_esti <- bsrmm$tembeta_save
beta_esti <- beta_esti[which(PPI > 0.5), (nburnin+2):(nburnin+niter+1)]
beta_esti[is.na(beta_esti)] <- 0
sum_beta <- colSums(beta_esti)

# save result for CrI and Mean 
beta_esti <- sweep(beta_esti, 1, STATS =  stand$Sx[which(PPI > 0.5)], FUN = "/")
beta_esti <- rbind(beta_esti, sum_beta)

# save the beta estimation result
beta_result <- data.frame(
  param = c(paste0("beta", which(PPI > 0.5)), "beta_sum"),
  mean = rowMeans(beta_esti),
  Lower_CrI = apply(beta_esti, 1, function(x) quantile(x, probs = 0.025)),
  Upper_CrI = apply(beta_esti, 1, function(x) quantile(x, probs = 0.975))
)

beta_result$ground_truth <- c(mydata$coeffs[which(PPI > 0.5),2], 0)

beta_result$mp <- proportion

# posterior mean for MVI
Y_mvi <- bsrmm$Y_mvi[(which(bsrmm$Ri==0)), (nburnin+2):(nburnin+niter+1)]
Y_true <- as.matrix(mydata$Y_true)
Y_true_train <- as.matrix(Y_true[train_ind, ])
Y_true_train <- Y_true_train[which(bsrmm$Ri==0), ]

mvi_result <- data.frame(
  mean = rowMeans(Y_mvi),
  Lower_CrI = apply(Y_mvi, 1, function(x) quantile(x, probs = 0.025)),
  Upper_CrI = apply(Y_mvi, 1, function(x) quantile(x, probs = 0.975)),
  ground_truth = Y_true_train
)

mvi_result$mp <- proportion

# result 
result <- list()
result$beta_result <- beta_result
result$mvi_result <- mvi_result

# display the PPI 
which(PPI > 0.5)

# save result
saveRDS(result, file = result_path)
