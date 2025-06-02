# set working directory 
setwd("/Users/kjiang/Desktop/BSRMM/simulation/")

# function
source("./function/gen_simdata_ind.R")
source("./function/gen_missing_value_ind.R")

## ------------------------
## Path
## ------------------------
save_data_path <- paste0("./raw_data/")

## ------------------------
## initial parameter
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

# imputation method
imputation_method <- Sys.getenv("imputation_method")
imputation_method

# # temporary parameter
# jobid <- 1
# missing_type <- "mnar0.33"
# n <- 300
# p <- 1000
# snr <- 1
# proportion <- 0.30
# imputation_method <- "low"

## ------------------------
## Data Generation
## ------------------------
set.seed(jobid)

mydata <-
  gen_missing_value_ind(
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

