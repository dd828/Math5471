
# Math5472 Project

This repository contains the final project for the course **MATH 5472**. The project implements five algorithms and two experiment functions and encapsulates them into an R package named **Math5472Project**.

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
Implements the original algorithm by Mazumder et al. (2010), as described in steps (2)–(4) of reference paper. It uses `rank_restricted_soft_svd` with $(\lambda = 0)$ to perform SVD for $(\hat{X})$.

### 6. `Timing_experiments`
This function performs timing experiments comparing three matrix completion algorithms: ALS (Alternating Least Squares),
Soft-Impute, and Soft-Impute-ALS. It generates a synthetic low-rank matrix, introduces noise and missing values,
applies the algorithms, and evaluates their performance in terms of computation time and relative error.

### 7. `RMSE_experiments`
This function compares the performance of three matrix completion algorithms—ALS, Soft-Impute, and Soft-Impute-ALS—
 based on the Root Mean Squared Error (RMSE). It generates a synthetic low-rank matrix, introduces noise and missing
 values, applies the specified algorithm(s), and computes the RMSE for both observed and full matrices.
 
## How to Reproduce the Results in the Report?

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

### Step 3: Perform timing experiments
Type the parameters you want to reproduce.
```{r}
timing_experiments(m = 100, n =50, r = 10, rank = 5, lambda = 15, missing_rate = 0.93, max_iter = 1000, tol = 1e-8, k = 1, sigma = 0.1)
```

### Step4: Perform RMSE experiments
Tyep teh parameters you want to reproduce.
```{r}
RMSE_experiments(m = 100, n = 50, r=10,rank = 5, noise_sd = 0.1, missing_rate = 0.2, lambda = 0.1, max_iter = 2000, method = "all")

```
## Notes

- Ensure that all required dependencies are installed in your R environment.
- For additional details about each function, please refer to the documentation provided within the package.








