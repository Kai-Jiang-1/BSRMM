# BSRMM

**Author**: Kai Jiang  
**Reference**: Jiang K, Saha S, Peterson CB. Bayesian sparse regression for microbiome-metabolite data integration.

## Overview:

This repository provides the implementation of **Bayesian Sparse Regression for Microbiome-Metabolite Data Integration (BSRMM)**. This method enables joint modeling of high-dimensional microbiome predictors and metabolomic outcomes. It allows for variable selection and outcome missing value imputation under a fully Bayesian framework.

## Directory Structure:

- `simulation/` – Contains scripts for generate simulated datasets and run the BSRMM model on them.
- `realdata/` – Contains scripts and preprocessed datasets for real-world applications.

## Key functions:
\begin{itemize}

\item `bsrmm`: The main function to fit the BSRMM model.

\item `Comp_lasso` and `BAZE`: Compared methods. 

\end{itemize}
