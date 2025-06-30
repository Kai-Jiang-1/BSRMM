# BSRMM

**Author**: Kai Jiang  
**Reference**: Jiang K, Saha S, Peterson CB. Bayesian sparse regression for microbiome-metabolite data integration.

## Overview

This repository provides the implementation of **Bayesian Sparse Regression for Microbiome-Metabolite Data Integration (BSRMM)**. This method enables joint modeling of high-dimensional microbiome predictors and metabolomic outcomes. The BSRMM model is implemented in R and is designed for variable selection and missing value imputation for outcomes.

## Directory Structure

- `simulation/` – Contains scripts to generate simulated datasets and run the BSRMM model.
- `realdata/` – Contains scripts and preprocessed datasets for real-world applications.

## Key Functions

- `bsrmm`: Main function to fit the BSRMM model.
- `Comp_lasso` and `BAZE`: Benchmark methods for comparison.

## Example

In the `simulation/` directory, you can find an example script that runs all functions, from data simulation to model evaluation.
