# set working directory 
setwd("/Users/kjiang/Desktop/BSRMM/simulation/")

# function
source("./function/gen_simdata_cor.R")
source("./function/gen_missing_value_cor.R")

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
# we have mar, mnar, mnar0.33 (1/3 mnar and 2/3 mar), and mnar0.67 (2/3 mnar and 1/3 mar)

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
# jobid <- 22
# missing_type <- "mnar0.33"
# n <- 300
# p <- 1000
# snr <- 1
# proportion <- 0.40
# imputation_method <- "low"

## ------------------------
## Data Generation
## ------------------------
set.seed(jobid)

mydata <-
  gen_missing_value_cor(
    proportion = proportion,
    n = n,
    p = p,
    snr = snr,
    type = missing_type
  )

# Which is true predictors
which(mydata$coeffs[,1] == 1)

# Total numver of true predictors
sum(mydata$coeffs[,1] == 1)

saveRDS(mydata, paste0(save_data_path, imputation_method, "_", mydata$data_type, "_cor_", n, "_", p, "_", snr, "_", proportion, "_", jobid, ".RDS"))
