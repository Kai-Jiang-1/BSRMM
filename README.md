README
================
Kai Jiang

# BSRMM

Bayesian Reduced Rank Regression for Microbiome-Metabolite Data
Integration

**Author**: Kai Jiang  
**Reference**: Jiang K, Saha S, Peterson CB. Bayesian sparse regression
for microbiome-metabolite data integration.

## Overview

This repository provides the implementation of **Bayesian Sparse
Regression for Microbiome-Metabolite Data Integration (BSRMM)**. This
method enables joint modeling of high-dimensional microbiome predictors
and metabolomic outcomes. The BSRMM model is implemented in R and is
designed for variable selection and missing value imputation for
outcomes.

## Directory Structure

- `code/` – Contains scripts to generate simulated datasets and run the
  BSRMM model.
- `data/` – Contains simulated datasets.
- `function` - Contains the implementation of the BSRMM model and
  simulation functions.
- `plot` - Contains MCMC diagnostic plots.
- `result` - Contains the results of the BSRMM model and diagnostic
  R-hat values.

## Example
