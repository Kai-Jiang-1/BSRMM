# set directory 
rm(list = ls())
setwd("/Users/kjiang/Desktop/BSRMM/")

source("./function/gen_simdata_ind.R")

## ------------------------
## Path
## ------------------------
save_data_path <- paste0("./data/")

# parameter
imputation_method <- "low"
missing_type <- "mnar0.33"
correlated <- "cor"
n <- 300
p <- 1000
snr <- 1
proportion <- 0.30
jobid <- 1

set.seed(jobid)
## ------------------------
## Data Simulation Function
## ------------------------
gen_missing_value <-
  function(proportion = proportion,
           n = n,
           p = p,
           snr = snr,
           type = missing_type) {
    
    # Initialization gamma 
    gammatrue <- rep(0, p)
    
    # index of covariate
    true_index <- c(1:3, 6:8)
    gammatrue[true_index] <- 1
    
    b <- rep(0, p)
    b[c(1:3, 6:8)] <- as.numeric(c(1,-0.8, 0.6,-1.5,-0.5, 1.2))
    
    # Standard deviation of X
    sigmaX <- 1
    
    # Noise variance
    sigma <- mean(abs(b[b != 0])) / snr  # 0.09333
    
    theta <- rep(0, p)
    theta[1:5] <- log(0.5 * p)
    # pirnt out theta
    print("Theta values:")
    print(theta[1:10])

    # generate simulated data
    sim_data <-
      gen_simdata_ind(n, p, gammatrue, b, theta, sigmaX, sigma)
    
    # extract data from sim_data
    Xorg <- sim_data$X

    # print out colMeans for Xorg[, 1:10]
    print("Column mean of Xorg[, 1:10]:")
    print(colMeans(Xorg[, 1:10]))

    temp <- exp(2 * Xorg)
    
    # Z <- apply(temp, 2, function(x) x / sum(x)) # convert to compositional data by column
    
    Z <-
      t(apply(temp, 1, function(x)
        x / sum(x))) # convert to compositional data by row
    
    X <- log(Z) # predictors
    
    beta <- sim_data$beta
    epsilon <- sim_data$epsilon
    
    Y_true <- X %*% b + epsilon
    
    # generate missing value
    if (type == "mnar") {
      cutoff = proportion
      lod <- quantile(Y_true, probs = c(cutoff)) # set lod
      Y_miss <- Y_true
      if (imputation_method == "low") {
        Y_miss[Y_miss < lod] <- log(0.5) + lod
        print(sum(Y_miss == min(Y_miss)) / length(Y_miss))
      } else if (imputation_method == "mean") {
        Y_miss[Y_miss < lod] <- NA
        input_value <- mean(Y_miss, na.rm = TRUE)
        Y_miss[which(is.na(Y_miss))] <- input_value
        print(sum(Y_miss == input_value) / length(Y_miss))
      }
      data_type <- "mnar"
      
    } else if (type == "mar") {
      Y_miss <- Y_true
      mar_idx <-
        sample(nrow(Y_miss), floor(nrow(Y_miss) * proportion))
      lod <- min(Y_miss)
      if (imputation_method == "low") {
        Y_miss[mar_idx] <- log(0.5) + lod
        print(sum(Y_miss == min(Y_miss)) / length(Y_miss))
      } else if (imputation_method == "mean") {
        Y_miss[mar_idx] <- NA
        input_value <- mean(Y_miss, na.rm = TRUE)
        Y_miss[which(is.na(Y_miss))] <- input_value
        print(sum(Y_miss == input_value) / length(Y_miss))
      }
      data_type <- "mar"
      
    } else if (type == "mnar0.33") {
      cutoff = proportion / 3
      lod <- quantile(Y_true, probs = c(cutoff))
      Y_miss <- Y_true
      if (imputation_method == "low") {
        Y_miss[Y_miss < lod] <- log(0.5) + lod
        remaining_indices <-
          setdiff(1:nrow(Y_true), which(Y_true < lod))
        mar_idx <-
          sample(remaining_indices, floor((proportion - cutoff) * n))
        Y_miss[mar_idx] <- log(0.5) + lod
        print(sum(Y_miss == min(Y_miss)) / length(Y_miss))
      } else if (imputation_method == "mean") {
        Y_miss[Y_miss < lod] <- NA
        remaining_indices <-
          setdiff(1:nrow(Y_true), which(Y_true < lod))
        mar_idx <- 
          sample(remaining_indices, floor((proportion - cutoff) * n))
        Y_miss[mar_idx] <- NA
        input_value <- mean(Y_miss, na.rm = TRUE)
        Y_miss[which(is.na(Y_miss))] <- input_value
        print(sum(Y_miss == input_value) / length(Y_miss))
      }
      data_type <- "mnar0.33"
      
    } else if (type == "mnar0.67") {
      cutoff = (2 * proportion) / 3
      lod <- quantile(Y_true, probs = c(cutoff))
      Y_miss <- Y_true
      if (imputation_method == "low") {
        Y_miss[Y_miss < lod] <- log(0.5) + lod
        remaining_indices <-
          setdiff(1:nrow(Y_true), which(Y_true < lod))
        mar_idx <-
          sample(remaining_indices, floor((proportion - cutoff) * n))
        Y_miss[mar_idx] <- log(0.5) + lod
        print(sum(Y_miss == min(Y_miss)) / length(Y_miss))
      } else if (imputation_method == "mean") {
        Y_miss[Y_miss < lod] <- NA
        remaining_indices <-
          setdiff(1:nrow(Y_true), which(Y_true < lod))
        mar_idx <- 
          sample(remaining_indices, floor((proportion - cutoff) * n))
        Y_miss[mar_idx] <- NA
        input_value <- mean(Y_miss, na.rm = TRUE)
        Y_miss[which(is.na(Y_miss))] <- input_value
        print(sum(Y_miss == input_value) / length(Y_miss))
      }
      data_type <- "mnar0.67"
      
    }
    
    coeffs <- cbind(gammatrue, beta)

    mydata <-
      list(
        coeffs = coeffs,
        Y_true = Y_true,
        X = X,
        Y_miss = Y_miss,
        data_type = data_type
      )
    
    return(mydata)
  }

## ------------------------
## Data Generation
## ------------------------
mydata <-
  gen_missing_value(
    snr = snr,
    proportion = proportion,
    n = n,
    p = p,
    type = missing_type
  )

# Which is true predictors
which(mydata$coeffs[,1] == 1)

# Total number of true predictors
sum(mydata$coeffs[,1] == 1)

saveRDS(mydata, paste0(save_data_path, imputation_method, "_", mydata$data_type, "_ind_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS"))

