# BSRMM

**Author**: Kai Jiang  
**Reference**: Jiang K, Saha S, Peterson CB. Bayesian sparse regression for microbiome-metabolite data integration.

## Overview:

This repository provides the implementation of **Bayesian Sparse Regression for Microbiome-Metabolite Data Integration (BSRMM)**. This method enables joint modeling of high-dimensional microbiome predictors and metabolomic outcomes. The BSRMM model is implemented in R and is designed for for variable selection and outcome missing value imputation.

## Directory Structure:

- `simulation/` – Contains scripts for generate simulated datasets and run the BSRMM model on them.
- `realdata/` – Contains scripts and preprocessed datasets for real-world applications.

## Key functions:
- `bsrmm`: The main function to fit the BSRMM model.

- `Comp_lasso` and `BAZE`: Compared methods. 

## Example: 
Under the `simulation/` directory, you can find a script named as example to run all functions from simualtioned data to evluate the perforamnce. 
