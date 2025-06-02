# set directory 
setwd("~/Desktop/BSRMM/realdata/sensitivity/")

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
# mydata_path <- "mydata_relative_abundance"
# outcome_name <- "Y_thi"
# imputation_method <- "mean"

## ------------------------
## import data
## ------------------------
mydata <- readRDS(paste0("~/Desktop/BSRMM/realdata/data_setup/", mydata_path, ".rds"))
ls(mydata)

# extract predictor 
predictor <- (mydata$X)
nrow(predictor)
ncol(predictor)

# for the real data, we need to make sure that rowSum is 1
# to make sure it is good for the package
sum(rowSums(predictor))

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

# set parameters
set.seed(123)

# fit the cv model to obtian the lam.min
comp <- cv.compCL(y = outcome, Z = predictor, Zc = NULL, intercept = TRUE)

# obtain the number of varaible selected
p <- ncol(predictor)
which(abs(coef(comp, s = "lam.min")[1:p]) > 0)

# TPR, FPR, TNR, FNR and MSE
if (sum((abs(coef(comp, s = "lam.min")[1:p]) > 0)) == 0) {
  nvs <- 0 

  nvs_idx <- which(abs(coef(comp, s = "lam.min")[1:p]) > 0)

  result <- list(
    nvs = nvs,
    nvs_idx = nvs_idx
  )
  
} else if (sum((abs(coef(comp, s = "lam.min")[1:p]) > 0)) > 0) {
    nvs <- sum((abs(coef(comp, s = "lam.min")[1:p]) > 0))
    
    nvs_idx <- which(abs(coef(comp, s = "lam.min")[1:p]) > 0)

    result <- list(
      nvs = nvs,
      nvs_idx = nvs_idx
    )
}

result_path <- paste0("./result/complasso_", mydata_path, "_", outcome_name, "_", imputation_method, ".rds")
saveRDS(result, result_path)
