
# Math5472 Project

This repository contains the final project for the course **MATH 5472**. The project implements five algorithms and encapsulates them into an R package named **Math5472Project**.

## Reference

The implementation is based on the paper:
**Matrix Completion and Low-Rank SVD via Fast Alternating Least Squares** .

## Implemented Functions

### 1. `rank_restricted_soft_svd`
Implements Algorithm 2.1 from the reference paper.

### 2. `rank_restricted_softImpute_ALS_3.1`
Implements Algorithm 3.1 from the reference paper.

### 3. `rank_restricted_softImpute_ALS_5.1`
Implements Algorithm 5.1 from the reference paper.

### 4. `als_matrix_completion`
Implements a modified version of Algorithm 5.2 from the reference paper, incorporating an additional regularization parameter $(\lambda\)$.

### 5. `soft_impute`
Implements the original algorithm by Mazumder et al. (2010), as described in steps (2)–(4) of their method. It uses `rank_restricted_soft_svd` with $(\lambda = 0)$ to perform SVD for $(\hat{X})$.

## How to Reproduce the Results in the Report

### Step 1: Install the Package

Open RStudio and install the package using the following command:
```R
# Install the package from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("dd828/Math5472_project")
```

### Step 2: Load the Package

Load the package into your R session:
```R
# Load the package
library(Math5472Project)
```

## Notes

- Ensure that all required dependencies are installed in your R environment.
- For additional details about each function, please refer to the documentation provided within the package.








