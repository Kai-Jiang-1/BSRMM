# set directory 
setwd("/Users/kjiang/Desktop/BSRMM/realdata/prediction")

# import package 
library(Compack)

## ------------------------
## inital parameter
## ------------------------
outcome_name <- Sys.getenv("outcome")
outcome_name

# imputation method
imputation_method <- Sys.getenv("imputation_method")
imputation_method

mydata_path <- Sys.getenv("mydata")
mydata_path

# temporary parameter
# outcome_name <- "Y_thi"
# a0 <- -6
# mydata_path <- "mydata_relative_abundance"
# imputation_method <- "low"
# Q_matrix <- "q_corr_par"

## ------------------------
## import data
## ------------------------
mydata <- readRDS(paste0("~/Desktop/BSRMM/realdata/data_setup/", mydata_path, ".rds"))
ls(mydata)

# extract predictor 
predictor <- (mydata$X)
nrow(predictor)
ncol(predictor)

# the package requires each rowsum is 1
rowSums(predictor)

# for zero, we assign 0.5 * observed minimum of each column 
predictor <- apply(predictor, 2, function(x) {
  min_nonzero <- min(x[x > 0])
  impute_nonzero <- 0.5 * min_nonzero
  x[x == 0] <- impute_nonzero
  return(x)
})

# redo the compostional 
predictor <- t(apply(predictor, 1, function(x) {
  x <- x / sum(x)
  return(x)
}))
sum(rowSums(predictor))

# extract outcome
outcome <- mydata[[outcome_name]]

# create dt by including outcome, statudy group and missing index 
dt <- data.frame(outcome = outcome, group = mydata$group)
dt$Ri <- 0
dt$Ri[which(dt$outcome != 0)] <- 1

if (imputation_method == "low") {
  outcome[which(outcome == 0)] <- 0.5 * sort(unique(outcome))[2] 
} else {
  outcome[which(outcome == 0)] <- mean(outcome[which(outcome != 0)])
}
outcome <- log(outcome)
dt$outcome <- outcome

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

# fit the cv model to obtian the lam.min
comp <- cv.compCL(y = outcome_train, Z = predictor_train, Zc = NULL, intercept = TRUE)

# obtain the number of varaible selected
p <- ncol(predictor)
which(abs(coef(comp, s = "lam.min")[1:p]) > 0)


# TPR, FPR, TNR, FNR and MSE
if (sum((abs(coef(comp, s = "lam.min")[1:p]) > 0)) == 0) {
  # calcualte mse based on the test data 
  Y_pred <- mean(outcome)
  
  mse <- mean((Y_pred - dt_test$outcome[which(dt_test$Ri!=0)])^2)
  
  result <- list(
    mse = mse
  )
  
} else if (sum((abs(coef(comp, s = "lam.min")[1:p]) > 0)) > 0) {
    Y_pred <- predict(comp, predictor_test)

    mse <- mean((Y_pred[which(dt_test$Ri!= 0), ] - as.matrix(dt_test$outcome)[which(dt_test$Ri!=0), ])^2)
    
    result <- list(
      mse = mse
    )
}

print(result$mse)

result_path <- paste0("./result/complasso_", mydata_path, "_", outcome_name, "_", imputation_method, ".rds")
saveRDS(result, result_path)
