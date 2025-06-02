## ------------------------
## Data Simulation Function (Independent covariates)
## ------------------------
gen_missing_value_ind <-
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
    
    # generate simulated data
    sim_data <-
      gen_simdata_ind(n, p, gammatrue, b, theta, sigmaX, sigma)
    
    # extract data from sim_data
    Xorg <- sim_data$X
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
