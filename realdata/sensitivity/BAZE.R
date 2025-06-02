# set directory 
setwd("~/Desktop/BSRMM/realdata/sensitivity/")

# function
source("./function/gibbsgamma.R")
source("./function/BayesFactor.R")

# required library
library(truncnorm)

## ------------------------
## inital parameter
## ------------------------
outcome_name <- Sys.getenv("outcome")
outcome_name

a0 <- as.numeric(Sys.getenv("a0"))
a0

mydata_path <- Sys.getenv("mydata")
mydata_path

imputation_method <- Sys.getenv("imputation_method")
imputation_method

Q_matrix <- Sys.getenv("Q_matrix")
Q_matrix

# temporary parameter
# outcome_name <- "Y_thi"
# a0 <- -6
# mydata_path <- "mydata_relative_abundance"
# imputation_method <- "mean"
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

# log trans
predictor <- log(predictor)

# scale and center predictor
stand <- list()
stand$mux <- colMeans(predictor)
stand$Sx <- apply(predictor, 2, sd)
predictor <- apply(predictor, 2, function(x) (x - mean(x))/sd(x))

# extract outcome 
outcome <- mydata[[outcome_name]]

# create dt by including outcome, statudy group and missing index 
dt <- data.frame(outcome = outcome, group = mydata$group)
dt$Ri <- 0
dt$Ri[which(dt$outcome != 0)] <- 1

# log-transformation for outcome 
# if imputation method is low, we assign 0.5 observed minimum
# if imputation method is mean, we assign the mean of observed value
if (imputation_method == "low") {
  outcome[which(outcome == 0)] <- 0.5 * sort(unique(outcome))[2] 
} else {
  outcome[which(outcome == 0)] <- mean(outcome[which(outcome != 0)])
}
outcome <- log(outcome)
dt$outcome <- outcome

# scale and center the outcome before split 
stand$muy <- mean(outcome)
stand$Sy <- sd(outcome)
outcome <- (outcome - stand$muy)/stand$Sy

# set parameters
seed <- 123

X <- as.matrix(predictor)
Y <- as.matrix(outcome)

n <- nrow(X)
p <- ncol(X)
c <- 100
N <- rbind(diag(x = 1, nrow = p, ncol = p), matrix(data = c, nrow = 1, ncol = p))

nburnin <- 10000
niter <- 20000

# initial number of predictors
nop <- floor(n/2)

# obtain the Q matrix 
Q <- mydata[[Q_matrix]]

# initial hyperparameters 
tau <- 1
nu <- 0
omega <- 0
a <- rep(a0, p) # parameter for the ising pripr, to controlling the number od predictors

baze <- gibbsgamma(nburnin, niter, p, nop, Y, X, N, a, Q, n, tau, nu, omega, seed, TRUE, stand, TRUE)

# selected features 
PPI <- colMeans(baze$gamma[(nburnin+2):(nburnin+niter+1),]) # posterior probabilities of inclusion
which(PPI > 0.5)
sum(PPI > 0.5)

# TPR, FPR, TNR, FNR and MSE
if (sum(PPI > 0.5) == 0) {
  nvs <- 0 

  nvs_idx <- which(PPI > 0.5)

  result <- list(
    nvs = nvs,
    a0 = a0,
    nvs_idx = nvs_idx
  )
  
} else {
  nvs <- sum(PPI > 0.5)

  nvs_idx <- which(PPI > 0.5)

  # result 
  result <- list(
    nvs = nvs,
    a0 = a0,
    nvs_idx = nvs_idx
  )
  }

# save the result
# result path need to changed later on server
result_path <- paste0("./result/baze_", mydata_path, "_", outcome_name, "_", imputation_method, "_", Q_matrix, "_", a0, ".rds")
saveRDS(result, result_path)
