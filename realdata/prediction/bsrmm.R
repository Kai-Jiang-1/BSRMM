# set directory 
setwd("/Users/kjiang/Desktop/BSRMM/realdata/prediction")

# functions 
source("./function/bsrmmgibbs_v2.R")
source("./function/bsrmmbf.R")

# required library
library(truncnorm)
library(mvnfast)
library(MASS)
library(coda)
library(MCMCpack)

## ------------------------
## inital parameter
## ------------------------
outcome_name <- Sys.getenv("outcome")
outcome_name

a0 <- as.numeric(Sys.getenv("a0"))
a0

mydata_path <- Sys.getenv("mydata")
mydata_path

Q_matrix <- Sys.getenv("Q_matrix")
Q_matrix

# temporary parameter
# outcome_name <- "Y_thi"
# a0 <- -6
# mydata_path <- "mydata_relative_abundance"
# Q_matrix <- "q_corr_par"

## ------------------------
## import data
## ------------------------
mydata <- readRDS(paste0("~/Desktop/BSRMM/realdata/data_setup/", mydata_path, ".rds"))
ls(mydata)
length(mydata)

# extract predictor 
predictor <- (mydata$X)
nrow(predictor)
ncol(predictor)
sum(rowSums(predictor))
which(predictor == 0)

# log transform
predictor <- log(predictor)

# scale and center predictor
stand <- list()
stand$mux <- colMeans(predictor)
stand$sx <- apply(predictor, 2, sd)
predictor <- apply(predictor, 2, function(x) (x - mean(x))/sd(x))

# extract outcome 
outcome <- mydata[[outcome_name]]

# log-transformation for outcome 
# for zero, we keep them as zero. O.W. we perform log-transformation
outcome[which(outcome != 0)] <- log(outcome[which(outcome != 0)])

# create dt by including outcome, statudy group and missing index 
dt <- data.frame(outcome = outcome, group = mydata$group)
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

outcome_train <- outcome[train_ind]
outcome_test <- outcome[-train_ind]

dt_train <- dt[train_ind,]
dt_test <- dt[-train_ind,]

# scale and center the predictor 
X <- as.matrix(predictor_train)

# set Y as the outcome_train 
Y <- as.matrix(outcome_train)

# preparation for the gibbs sampling 
n <- nrow(X)
p <- ncol(X)
c <- 100

# construct the T matrix, because we cannot assign a matrix to T due to R's limitation 
# we name it as N matrix 
N <- rbind(diag(x = 1, nrow = p, ncol = p), matrix(data = c, nrow = 1, ncol = p))

# number of sampling procedure 
nburnin <- 10000
niter <- 20000

# initial number of predictors
nop <- floor(n/2)

# obtain the Q matrix 
Q <- mydata[[Q_matrix]]
print(sum(Q))
diag(Q) <- 0
print(sum(Q))

# initial hyperparameters 
tau <- 1
nu <- 0
omega <- 0
a <- rep(a0, p) # parameter for the ising pripr, to controlling the number od predictors

predict_result = TRUE
display = TRUE
bayestmetab = TRUE
# debug(bsrmmgibbs)
bsrmm <- bsrmmgibbs(nburnin, niter, p, nop, Y, X, N, a, Q, n, tau, nu, omega, seed, predict_result, bayestmetab, display)

# select features
PPI <- colMeans(bsrmm$gamma[(nburnin+2):(nburnin+niter+1),]) # posterior probabilities of inclusion
which(PPI > 0.5)
sum(PPI > 0.5)

# check the intercept
beta_0 <- mean(bsrmm$beta_0[(nburnin+2):(nburnin+niter+1)]) 
beta_0

# save the results
if (sum(PPI > 0.5) == 0) { # no variable is selectied
  Y_pred <- mean(bsrmm$beta_0[(nburnin+2):(nburnin+niter+1)]) # + stand$muy

  mse <- mean((Y_pred - dt_test$outcome[which(dt_test$Ri!=0)])^2)

  result <- list(
    mse = mse
  )
} else {  
  # step 1:calculate beta 
  Y_mvi <- as.matrix(rowMeans(bsrmm$Y_mvi[, (nburnin+2):(nburnin+niter+1)])) # using imputed value
  Xri <- X[, which(PPI > 0.5)] # select the predictors
  Tri <- N[, which(PPI > 0.5)] # select the predictors

  # Lambda = I - J/(n+1)
  Lambda <- diag(1, nrow = n, ncol = n) - matrix(1, nrow = n, ncol = n)/(n + 1)
  
  # Ari = X^T (I - J/(n+1)) X + sigma^(-2) tau^(-2) T^T T
  Ari <- t(Xri) %*% Lambda %*% Xri + tau^(-2) * t(Tri) %*% Tri
  invAri <- chol2inv(chol(Ari))
 
  # beta = Ari^(-1) X^T (I - J/(n+1)) Y
  beta <- invAri %*% t(Xri) %*% Lambda %*% Y_mvi
  
  # step 2: scale and center the test data
  X_test <- as.matrix(predictor_test)

  # step 3: calculate the predicted value, preidcted value is based on the log(outcome)
  Y_pred <- X_test[, which(PPI > 0.5)] %*% beta + mean(bsrmm$beta_0[(nburnin+2):(nburnin+niter+1)]) # + stand$muy

  # step 4: calculate the mse
  mse <- mean((Y_pred[which(dt_test$Ri!= 0), ] - as.matrix(dt_test$outcome)[which(dt_test$Ri!=0), ])^2)
  
  # result 
  result <- list(
    mse = mse
  )
  }

print(result$mse)

# save the result
# result path need to changed later on server
result_path <- paste0("./result/bsrmm_", mydata_path, "_", outcome_name, "_", Q_matrix, "_", a0, ".rds")
saveRDS(result, result_path)
