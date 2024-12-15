#' RMSE Experiments for Matrix Completion Algorithms
#'
#' This function compares the performance of three matrix completion algorithms—ALS, Soft-Impute, and Soft-Impute-ALS—
#' based on the Root Mean Squared Error (RMSE). It generates a synthetic low-rank matrix, introduces noise and missing
#' values, applies the specified algorithm(s), and computes the RMSE for both observed and full matrices.
#'
#' @param m Integer. Number of rows in the matrix.
#' @param n Integer. Number of columns in the matrix.
#' @param r Integer. Operator rank, the dimension of the low-rank approximation used in algorithms.
#' @param rank Integer. Target rank of the matrix (ground truth rank of the synthetic matrix).
#' @param noise_sd Numeric. Standard deviation of the noise added to the low-rank matrix.
#' @param missing_rate Numeric. Proportion of elements in the matrix that are set to missing (NA).
#' @param lambda Numeric. Regularization parameter for the algorithms.
#' @param max_iter Integer. Maximum number of iterations for the algorithms.
#' @param method Character. Which method to use for matrix completion. Options are "all", "ALS", "softImpute", or "softImputeALS". Default is "all".
#' @param tol Numeric. Convergence tolerance for stopping the algorithms.
#'
#'  @return A list containing:
#'   - `RMSE_results`: A named list with the RMSE values for each method (`observed_rmse` and `full_rmse`).
#'   - The function also prints the RMSE results and plots the reconstructed matrices for visual comparison.
#' @details
#' RMSE is computed on both the observed entries (those that were not missing) and the full matrix (including missing values).
#'
#' @examples
#' # Example usage:
#' RMSE_experiments(
#'   m = 100,
#'   n = 50,
#'   r = 5,
#'   rank = 5,
#'   noise_sd = 0.1,
#'   missing_rate = 0.2,
#'   lambda = 0.1,
#'   max_iter = 2000,
#'   method = "all"
#' )
#'
#' @import ggplot2
#' @export
RMSE_experiments <- function(m, n, r,rank, noise_sd, missing_rate, lambda, max_iter, tol=1e-6,method = "all") {

  # Simulate ground truth matrix
  set.seed(42)
  A_true <- matrix(rnorm(m * rank), m, rank)
  B_true <- matrix(rnorm(n * rank), n, rank)
  X_true <- A_true %*% t(B_true)  # Low-rank matrix

  # Add noise
  X_noisy <- X_true + matrix(rnorm(m * n, sd = noise_sd), m, n)

  # Introduce missing values
  X_missing <- X_noisy
  missing_indices <- sample(length(X_missing), size = missing_rate * length(X_missing))
  X_missing[missing_indices] <- NA

  # Initialize results list
  results <- list()

  if (method == "all" || method == "ALS") {
    # Apply ALS algorithm
    als_results <- als_matrix_completion(X_missing, r = r, lambda = lambda, max_iter = max_iter,tol=tol)
    X_reconstructed_als <- als_results$A %*% t(als_results$B)

    # Evaluate accuracy for ALS
    observed_indices <- which(!is.na(X_missing), arr.ind = TRUE)
    observed_rmse_als <- sqrt(mean((X_reconstructed_als[observed_indices] - X_noisy[observed_indices])^2))
    full_rmse_als <- sqrt(mean((X_reconstructed_als - X_true)^2))

    results$ALS <- list(observed_rmse = observed_rmse_als, full_rmse = full_rmse_als, X_reconstructed = X_reconstructed_als)
  }

  if (method == "all" || method == "softImpute") {
    # Apply Soft-Impute algorithm
    soft_impute_results <- soft_impute(X_missing, r =  min(m,n), lambda = lambda, max_iter = max_iter,tol=tol)
    X_reconstructed_soft_impute <- soft_impute_results$M_hat

    # Evaluate accuracy for Soft-Impute
    observed_rmse_soft_impute <- sqrt(mean((X_reconstructed_soft_impute[observed_indices] - X_noisy[observed_indices])^2))
    full_rmse_soft_impute <- sqrt(mean((X_reconstructed_soft_impute - X_true)^2))

    results$SoftImpute <- list(observed_rmse = observed_rmse_soft_impute, full_rmse = full_rmse_soft_impute, X_reconstructed = X_reconstructed_soft_impute)
  }

  if (method == "all" || method == "softImputeALS") {
    # Apply Soft-Impute-ALS algorithm
    soft_impute_als_results <- rank_restricted_softImpute_ALS_5.1(X_missing, r = r, lambda = lambda, max_iter = max_iter,tol=tol)
    X_reconstructed_soft_impute_als <- soft_impute_als_results$A %*% t(soft_impute_als_results$B)

    # Evaluate accuracy for Soft-Impute-ALS
    observed_rmse_soft_impute_als <- sqrt(mean((X_reconstructed_soft_impute_als[observed_indices] - X_noisy[observed_indices])^2))
    full_rmse_soft_impute_als <- sqrt(mean((X_reconstructed_soft_impute_als - X_true)^2))

    results$SoftImputeALS <- list(observed_rmse = observed_rmse_soft_impute_als, full_rmse = full_rmse_soft_impute_als, X_reconstructed = X_reconstructed_soft_impute_als)
  }

  # Print RMSE results
  for (method_name in names(results)) {
    cat(paste0("RMSE for ", method_name, ":\n"))
    cat("  RMSE on observed entries:", results[[method_name]]$observed_rmse, "\n")
    cat("  RMSE on full matrix:", results[[method_name]]$full_rmse, "\n")
    cat("\n")
  }

  # Optional: Visualization of true vs reconstructed matrices for each method
  par(mfrow = c(3, 2))
  for (method_name in names(results)) {
    image(results[[method_name]]$X_reconstructed, main = paste0(method_name, " Reconstructed Matrix"))
    image(X_true, main = "Original Matrix")
  }

}

# Example usage:
# matrix_completion_with_comparison(m = 100, n = 50, r = 5, noise_sd = 0.1, missing_rate = 0.2, lambda = 0.1, max_iter = 2000, max_iter_svd = 100, method = "all")

