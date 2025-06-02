# set directory 
setwd("~/Desktop/BSRMM/realdata/sensitivity/")

# functions 
source("./function/bsrmmgibbs.R")
source("./function/bsrmmbf.R")

# required library
library(truncnorm)
library(mvnfast)
library(MASS)
library(coda)
library(MCMCpack)
library(phyloseq)

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

# log trans
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

# set parameters
seed <- 123

X <- as.matrix(predictor)
Y <- as.matrix(outcome)

# preparation for the gibbs sampling 
n <- nrow(X)
p <- ncol(X)
c <- 100
N <- rbind(diag(x = 1, nrow = p, ncol = p), matrix(data = c, nrow = 1, ncol = p))

# number of sampling procedure 
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

predict_result = TRUE
display = TRUE
bayestmetab = TRUE
# debug(bsrmmgibbs)
bsrmm <- bsrmmgibbs(nburnin, niter, p, nop, Y, X, N, a, Q, n, tau, nu, omega, seed, predict_result, bayestmetab, display)

# select features
PPI <- colMeans(bsrmm$gamma[(nburnin+2):(nburnin+niter+1),]) # posterior probabilities of inclusion
which(PPI > 0.5)
sum(PPI > 0.5)
# 10 128 159 171 180 212 279 322
# save the results
if (sum(PPI > 0.5) == 0) { # no variable is selectied
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

# create the foreast plot 

if (sum(PPI > 0.5) > 0) {
  beta_est <- bsrmm$tembeta_save[which(PPI > 0.5), ((nburnin+2):(nburnin+niter+1))]
}

# save the result
# result path need to changed later on server
result_path <- paste0("./result/bsrmm_", mydata_path, "_", outcome_name, "_", Q_matrix, "_", a0, ".rds")
saveRDS(result, result_path)