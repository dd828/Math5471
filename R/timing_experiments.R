

#' Timing Experiments for Matrix Completion Algorithms
#'
#' This function performs timing experiments comparing three matrix completion algorithms: ALS (Alternating Least Squares),
#' Soft-Impute, and Soft-Impute-ALS. It generates a synthetic low-rank matrix, introduces noise and missing values,
#' applies the algorithms, and evaluates their performance in terms of computation time and relative error.
#'
#' @param m Integer. Number of rows in the matrix.
#' @param n Integer. Number of columns in the matrix.
#' @param r Integer. Operator rank, the dimension of the low-rank approximation used in algorithms.
#' @param rank Integer. Target rank of the matrix (ground truth rank of the synthetic matrix).
#' @param lambda Numeric. Regularization parameter for penalizing large values in the matrix factorization.
#' @param missing_rate Numeric. Proportion of elements in the matrix that are set to missing (NA).
#' @param max_iter Integer. Maximum number of iterations for the optimization algorithms.
#' @param tol Numeric. Convergence tolerance for stopping the algorithms.
#' @param k Integer. An integer for recording cumulative time and objective loss every `k` steps. Default is 1.
#' @param sigma Numeric. Standard deviation of the noise added to the synthetic low-rank matrix.
#'
#' @return A ggplot2 object showing the relative error as a function of computation time for the three algorithms.
#'
#' @details
#' Relative error and computation times are recorded for each algorithm. The results are plotted, showing the
#' trade-offs between time and accuracy.
#'
#' @examples
#' # Example usage:
#' timing_experiments(
#'   m = 100,
#'   n = 100,
#'   r = 5,
#'   rank = 5,
#'   lambda = 0.1,
#'   missing_rate = 0.2,
#'   max_iter = 100,
#'   tol = 1e-4,
#'   k = 10,
#'   sigma = 0.1
#' )
#'
#' @import ggplot2
#' @export
timing_experiments <- function(m, n, r, rank, lambda, missing_rate, max_iter, tol, k, sigma) {


  set.seed(42)

  # Ground truth matrix
  A_true <- matrix(rnorm(m * rank), m, rank)
  B_true <- matrix(rnorm(n * rank), n, rank)
  X_true <- A_true %*% t(B_true)  # Low-rank matrix

  # Add noise
  noise_sd <- sigma
  X_noisy <- X_true + matrix(rnorm(m * n, sd = noise_sd), m, n)

  # Introduce missing values
  X_obs <- X_noisy
  missing_indices <- sample(length(X_obs), size = missing_rate * length(X_obs))
  X_obs[missing_indices] <- NA

  # Run ALS
  als_results <- als_matrix_completion(X_obs, r, lambda, max_iter = max_iter, tol = tol, k = k)

  # Run Soft-Impute
  soft_impute_results <- soft_impute(X_obs, r, lambda, max_iter = max_iter, tol = tol, k = k)

  # Run Soft-Impute-ALS
  soft_impute_als_results <- rank_restricted_softImpute_ALS_5.1(X_obs, r, lambda, max_iter = max_iter, tol = tol, k = k)

  # Prepare data for plotting
  relative_als <- abs(als_results$objectives - als_results$obj_min) / als_results$obj_min
  relative_soft_impute <- abs(soft_impute_results$objectives - soft_impute_results$obj_min) / soft_impute_results$obj_min
  relative_soft_impute_als <- abs(soft_impute_als_results$objectives - soft_impute_als_results$obj_min) / soft_impute_als_results$obj_min

  df_als <- data.frame(Time = als_results$times, RelativeError = relative_als, Method = "ALS")
  df_soft_impute <- data.frame(Time = soft_impute_results$times, RelativeError = relative_soft_impute, Method = "Soft-Impute")
  df_soft_impute_als <- data.frame(Time = soft_impute_als_results$times, RelativeError = relative_soft_impute_als, Method = "Soft-Impute-ALS")

  df <- rbind(df_als, df_soft_impute, df_soft_impute_als)

  # Plot the data
  library(ggplot2)

  ggplot(df, aes(x = Time, y = RelativeError, color = Method)) +
    geom_line() +
    geom_point() +
    scale_y_log10() +
    labs(
      title = paste0("(", m, ", ", n, ") ", missing_rate * 100, "% NAs \u03bb=", lambda, " rank=", rank, " r=", r),
      x = "Time in Seconds",
      y = "Relative Objective "
    ) +
    theme_minimal()
}

