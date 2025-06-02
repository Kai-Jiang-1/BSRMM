#' BRSMM (Bayesian sparse regression for microbiome-metabolite data integration)
#' Kai Jiang
#' Description: In this project, we would like perfrom MCMC diagnoistics for the BRSMM model 

# set directory 
rm(list = ls())
setwd("/Users/kjiang/Desktop/BSRMM/simulation/")

# functions 
source("./function/bsrmmgibbs.R")
source("./function/bsrmmbf.R")

# required library
library(truncnorm)
library(mvnfast)
library(MASS)
library(coda)
library(MCMCpack)
library(knitr)
library(kableExtra)

# parameter setting
imputation_method <- "low"
missing_type <- "mnar0.33"
correlated <- "cor"
n <- 300
p <- 1000
snr <- 1
proportion <- 0.40
jobid <- 22

## ------------------------
## data/result path
## ------------------------
if (correlated == "cor") {
  mydata_path <- paste0("./raw_data/", imputation_method, "_", missing_type, "_cor_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
  result_path <- paste0("./output/bsrmm_", imputation_method, "_imputation_", missing_type, "_cor_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
} else {
  mydata_path <- paste0("./raw_data/", imputation_method, "_", missing_type, "_ind_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
  result_path <- paste0("./output/bsrmm_", imputation_method, "_imputation_", missing_type, "_ind_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
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
set.seed(1996)

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
display = TRUE
bayestmetab = TRUE

# step 1: run four chains for the model
all_chains_bsrmm_output <- list()
num_chains <- 4
seed_chains <- c(2, 22, 222, 2222)
for (i in 1:num_chains) {
  seed <- seed_chains[i]
  cat(paste0("Running Chain ", i, "...\n"))
  
  all_chains_bsrmm_output[[i]] <- bsrmmgibbs(nburnin, niter, p, nop, Y, X, N, a, Q, n, tau, nu, omega, seed, predict_result, bayestmetab, display)
  
  # name the list 
  names(all_chains_bsrmm_output)[i] <- paste0("Chain_", i)
}


total_gamma_sum <- rep(0, p)
for (chain_idx in 1:num_chains) {
  current_bsrmm_output <- all_chains_bsrmm_output[[chain_idx]]

  print(which(colMeans(current_bsrmm_output$gamma[ (nburnin + 2) : (nburnin + niter + 1), ]) > 0.5))
  total_gamma_sum <- total_gamma_sum + colMeans(current_bsrmm_output$gamma[ (nburnin + 2) : (nburnin + niter + 1), ])
}

PPI <- total_gamma_sum / num_chains

which(PPI > 0.5)

# step 2: prepare the result for diagnostics 
all_chains_scalar_samples <- list()
all_chains_beta_samples <- list()

for (chain_idx in 1:num_chains) {
  current_bsrmm_output <- all_chains_bsrmm_output[[chain_idx]]

  # sigma2
  sigma2_samples <- current_bsrmm_output$sigma2[(nburnin + 2):(nburnin + niter + 1)]
  
  # beta_0
  beta_0_samples <- as.numeric(current_bsrmm_output$beta_0[(nburnin + 2):(nburnin + niter + 1)])

  # alpha
  alpha_samples <- current_bsrmm_output$alpha[(nburnin + 2):(nburnin + niter + 1)]
  
  # beta
  tembeta_samples <- current_bsrmm_output$tembeta_save[, (nburnin + 2):(nburnin + niter + 1)]
  tembeta_samples[is.na(tembeta_samples)] <- 0
  tembeta_samples <- tembeta_samples[which(PPI > 0.5), ]

  mcmc_data_matrix_current <- data.frame(
    sigma2 = sigma2_samples,
    alpha = alpha_samples,
    beta0 = beta_0_samples
  )
  
  # create t_tembeta_samples as data.frame, each column name based on beta_(which ppi > 0.5)
  t_tembeta_samples <- data.frame(t(tembeta_samples))
  colnames(t_tembeta_samples) <- paste0("beta ", which(PPI > 0.5))
  
  all_chains_scalar_samples[[chain_idx]] <- coda::mcmc(mcmc_data_matrix_current)
  all_chains_beta_samples[[chain_idx]] <- coda::mcmc(t_tembeta_samples)
  
  }

# Combine samples from all chains into mcmc.list objects
mcmc_scalar_list <- coda::as.mcmc.list(all_chains_scalar_samples)
mcmc_beta_list <- coda::as.mcmc.list(all_chains_beta_samples)

# step 3: MCMC diagnostics
# Summary Statistics
summary(mcmc_scalar_list)
summary(mcmc_beta_list)

# Convergence Diagnostics (Gelman-Rubin statistic)
gelman.diag(mcmc_scalar_list, autoburnin = FALSE)
gelman.diag(mcmc_beta_list, autoburnin = FALSE)

# save at the data frame 
rhat_scalar <- data.frame(
  parameter = c("sigma2", "alpha", "beta 0"),
  rhat = gelman.diag(mcmc_scalar_list, autoburnin = FALSE)$psrf[, 1],
  upper = gelman.diag(mcmc_scalar_list, autoburnin = FALSE)$psrf[, 2]
)

rhat_scalar[, 2:3] <- round(rhat_scalar[, 2:3], 3)

rhat_beta <- data.frame(
  parameter = colnames(mcmc_beta_list[[1]]),
  rhat = gelman.diag(mcmc_beta_list, autoburnin = FALSE)$psrf[, 1],
  upper = gelman.diag(mcmc_beta_list, autoburnin = FALSE)$psrf[, 2]
)

rhat_beta[, 2:3] <- round(rhat_beta[, 2:3], 3)

rhat <- rbind(rhat_scalar, rhat_beta)
row.names(rhat) <- NULL

cat(
  kable(rhat, "latex", booktabs = TRUE, row.names = FALSE, caption = "Gelman-Rubin Statistic") %>%
    kable_styling(latex_options = c("striped", "hold_position", "scale_down")) %>%
    row_spec(0, bold = TRUE))

pdf("/Users/kjiang/Library/CloudStorage/OneDrive-InsideMDAnderson/manuscript/traceplot_bsrmm.pdf", width = 12, height = 10)
par(mfrow = c(2, 3)) 
traceplot(mcmc_scalar_list)
traceplot(mcmc_beta_list)
dev.off()

# pick one chain for prediction result
bsrmm1 <- all_chains_bsrmm_output[[1]]

# plot between our imputed missing value vs the ground truth
Y_mvi <- as.matrix(rowMeans(bsrmm1$Y_mvi[, (nburnin+2):(nburnin+niter+1)]))
Y_mvi <- Y_mvi[which(bsrmm1$Ri==0), ]

Y_true <- as.matrix(mydata$Y_true)
Y_true_train <- as.matrix(Y_true[train_ind, ])
Y_true_test <- as.matrix(Y_true[-train_ind, ])
Y_true_train <- Y_true_train[which(bsrmm1$Ri==0), ]

# plot the imputed value vs the ground truth by using ggplot 
library(ggplot2)

ggplot(data = data.frame(x = Y_true_train, y = Y_mvi), aes(x = x, y = y)) +
  geom_point() +
  labs(x = "Ground truth values", y = "Imputed values from BSRMM") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1.5) +
  theme_classic(base_size = 16) 

ggsave("/Users/kjiang/Library/CloudStorage/OneDrive-InsideMDAnderson/manuscript/imputed_vs_true_bsrmm.pdf", width = 8, height = 6)
