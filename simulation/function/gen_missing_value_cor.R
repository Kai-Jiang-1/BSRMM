## ------------------------
## Data Simulation Function (Depedent covariates)
## ------------------------
gen_missing_value_cor <- 
  function(proportion = proportion,
           n = n,
           p = p,
           snr = snr, 
           type = missing_type) {
    
    # Initialization gamma 
    gammatrue <- rep(0, p)
    true_index <- c(seq(180, 400, by = 20), seq(580, 800, by = 20))
    gammatrue[true_index] <- 1
    
    # Initialization beta
    b <- rep(0, p)
    b1 <- as.numeric(c(0.8800, -1.3900, 1.0400, 1.2100, -1.8600, -1.3400, 1.7600, -0.9900, 0.6900, -0.5400, 1.3500, -0.8100))
    b2 <- as.numeric(c(-1.4100, -1.1500, 0.5100, -1.9500, 1.9300, -0.8500, -1.6600, 1.4800, 1.8700, 0.7200, 0.6700, -0.1600))
    b[true_index[seq(1, length(true_index), by=2)]] <- b1
    b[true_index[seq(2, length(true_index), by=2)]] <- b2
    
    # Standard deviation of X
    sigmaX <- rep(1, p)
    sigma <- mean(c(abs(b1), abs(b2))) / snr 
    
    Xcor <- diag(0,p)
    
    # First set of conditions
    for(i in seq(180, 380, by = 20)) {
      for(j in seq(i + 20, 400, by = 20)) {
        Xcor[i, j] <- 0.75 - 0.0015 * abs(i - j)
      }
    }
    
    # Second set of conditions
    for(i in seq(580, 780, by = 20)) {
      for(j in seq(i + 20, 800, by = 20)) {
        Xcor[i, j] <- 0.75 - 0.0015 * abs(i - j)
      }
    }
    
    # Third set of conditions
    for(i in 445:459) {
      for(j in (i + 1):460) {
        Xcor[i, j] <- 0.4 - 0.02 * abs(i - j)
      }
    }
    
    # Fourth set of conditions
    for(i in 945:959) {
      for(j in (i + 1):960) {
        Xcor[i, j] <- 0.4 - 0.02 * abs(i - j)
      }
    }
    Xcor <- Xcor + t(Xcor) + diag(1, p)
    
    theta <- rep(0, p)
    theta[c(seq(180, 400, by=20), seq(580, 800, by=20))] <- log(0.5*p)
    theta[c(45:60, 445:460, 945:960)] <- log(0.25*p)
    
    # generate simulated data
    sim_data <-
      gen_simdata_cor(n, gammatrue, b, theta, sigmaX, Xcor, sigma)
    
    # extract data from sim_data
    Xorg <- sim_data$X
    temp <- exp(2 * Xorg)
    
    Z <-
      t(apply(temp, 1, function(x)
        x / sum(x))) # convert to compositional data 
    
    X <- log(Z) # predictors
    
    beta <- sim_data$beta
    epsilon <- sim_data$epsilon
    
    Y_true <- X %*% beta + epsilon
    
    # generate missing value
    if (type == "mnar") { # missing below the LOD
      cutoff = proportion
      lod <- quantile(Y_true, probs = c(cutoff)) # set lod
      Y_miss <- Y_true
      if (imputation_method == "low") {
        Y_miss[Y_miss < lod] <- log(0.5) + lod
        print(sum(Y_miss == min(Y_miss)) / length(Y_miss))
      } else if (imputation_method == "mean"){
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
      } else if (imputation_method == "mean"){
        Y_miss[mar_idx] <- NA
        input_value <- mean(Y_miss, na.rm = TRUE)
        Y_miss[which(is.na(Y_miss))] <- input_value
        print(sum(Y_miss == input_value) / length(Y_miss))
      }
      data_type <- "mar"
      
    } else if (type == "mnar0.33") { # combine with 1/3 mnar and 2/3 mar 
      cutoff = proportion / 3
      lod <- quantile(Y_true, probs = c(cutoff))
      Y_miss <- Y_true
      if (imputation_method == "low") {
        Y_miss[Y_miss < lod] <- log(0.5) + lod
        remaining_indices <-
          setdiff(1:nrow(Y_true), which(Y_true < lod))
        mar_idx <-
          sample(remaining_indices, floor((proportion - cutoff)*n))
        Y_miss[mar_idx] <- log(0.5) + lod
        print(sum(Y_miss == min(Y_miss)) / length(Y_miss))
      } else if (imputation_method == "mean") {
        Y_miss[Y_miss < lod] <- NA
        remaining_indices <-
          setdiff(1:nrow(Y_true), which(Y_true < lod))
        mar_idx <- 
          sample(remaining_indices, floor((proportion - cutoff)*n))
        Y_miss[mar_idx] <- NA
        input_value <- mean(Y_miss, na.rm = TRUE)
        Y_miss[which(is.na(Y_miss))] <- input_value
        print(sum(Y_miss == input_value) / length(Y_miss))
      }
      data_type <- "mnar0.33"
      
    } else if (type == "mnar0.67") { # combine with 2/3 mnar and 1/3 mar
      cutoff = (2 * proportion) / 3
      lod <- quantile(Y_true, probs = c(cutoff))
      Y_miss <- Y_true
      if (imputation_method == "low") {
        Y_miss[Y_miss < lod] <- log(0.5) + lod
        remaining_indices <-
          setdiff(1:nrow(Y_true), which(Y_true < lod))
        mar_idx <-
          sample(remaining_indices, floor((proportion - cutoff)*n))
        Y_miss[mar_idx] <- log(0.5) + lod
        print(sum(Y_miss == min(Y_miss)) / length(Y_miss))
      } else if (imputation_method == "mean") {
        Y_miss[Y_miss < lod] <- NA
        remaining_indices <-
          setdiff(1:nrow(Y_true), which(Y_true < lod))
        mar_idx <- 
          sample(remaining_indices, floor((proportion - cutoff)*n))
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