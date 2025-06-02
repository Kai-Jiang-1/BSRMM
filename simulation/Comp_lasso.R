# set directory 
setwd("/Users/kjiang/Desktop/BSRMM/simulation/")

# import package 
library(Compack)

## ------------------------
## inital parameter
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
    result_path <- paste0("./output/complasso_", imputation_method, "_imputation_", missing_type, "_cor_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
} else {
    mydata_path <- paste0("./raw_data/", imputation_method, "_", missing_type, "_ind_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
    result_path <- paste0("./output/complasso_", imputation_method, "_imputation_", missing_type, "_ind_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS")
}

# import data 
mydata <- readRDS(mydata_path)

# extract predictor and outcome 
predictor <- mydata$X
outcome <- mydata$Y_miss
dt <- data.frame(outcome)

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

# check if the sum of each row is 1
# we need to exp transform the data to use the package
rowSums(exp(predictor_train))

# fit the cv model to obtian the lam.min
comp <- cv.compCL(y = outcome_train, Z = exp(predictor_train), Zc = NULL, intercept = TRUE)

# obtain the number of varaible selected
which(abs(coef(comp, s = "lam.min")[1:p]) > 0)
beta_select <- rep(0,ncol(predictor))
beta_select[which(abs(coef(comp, s = "lam.min")[1:p]) > 0)] <- 1

# TPR, FPR, TNR, FNR and MSE
if (sum((abs(coef(comp, s = "lam.min")[1:p]) > 0)) == 0) {
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
    y_pred <- mean(dt_test$outcome)
    
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
        observed_data_mse = observed_data_mse,
        coef_l2 = coef_l2
        )
} else if (sum((abs(coef(comp, s = "lam.min")[1:p]) > 0)) > 0) {
    # confusion matrix
    true_beta <- mydata$coeffs[,1]
    t <- table(beta_select, true_beta)
    
    TPR <- t[2,2]/(sum(t[,2]))
    FNR <- 1 - TPR
    TNR <- t[1,1]/(sum(t[,1]))
    FPR <- 1 - TNR

    TP <- t[2,2]

    FN <- t[1,2]

    TN <- t[1,1]

    FP <- t[2,1]

    # predict y test 
    Y_pred <- predict(comp, exp(predictor_test))

    # mse calculated by 
    observed_data_mse <- mean((Y_pred[which(dt_test$Ri!= 0), ] - as.matrix(dt_test$outcome)[which(dt_test$Ri!=0), ])^2)

    # obtain the F1 score
    F1_score <- TP/(TP + 0.5*(FP + FN))

    # l2 loss 
    beta_hat <- coef(comp, s= "lam.1se")[1:p]
    coef_l2 = sqrt(sum((mydata$coeffs[,2] - beta_hat)^2))

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
        observed_data_mse = observed_data_mse,
        coef_l2 = coef_l2
        )
}

result$nrmse <- nrmse

saveRDS(result, file = result_path)
