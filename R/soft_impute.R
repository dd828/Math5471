#' Soft-Impute Algorithm for Matrix Completion
#'
#' Implements the Soft-Impute algorithm for matrix completion. The algorithm uses Rank Restricted Soft SVD
#'algorithm with lambda=0 to do svd for X_hat.
#'
#'
#' @param X A numeric matrix of size \code{m x n} with missing values (use \code{NA} for missing entries).
#' @param r An integer specifying the target rank for the low-rank approximation.
#' @param lambda A numeric value for the regularization parameter, controlling the soft-thresholding applied to singular values.
#' @param max_iter An integer specifying the maximum number of iterations for the Soft-Impute algorithm. Default is 100.
#' @param max_iter_svd An integer specifying the maximum number of iterations for the rank-restricted SVD solver. Default is 100.
#' @param tol A numeric value for the convergence tolerance. The algorithm stops if the change in the objective
#'   loss between iterations is less than \code{tol}. Default is \code{1e-6}.
#' @param k An integer specifying the interval for recording cumulative time and objective loss.
#'   Metrics are recorded every \code{k} steps or at convergence. Default is 1.
#'
#' @return A list containing the following elements:
#'   \item{M_hat}{The completed matrix of size \code{m x n}, with missing values imputed.}
#'   \item{times}{A numeric vector of cumulative times (in seconds) recorded at intervals or at convergence.}
#'   \item{objectives}{A numeric vector of objective loss values recorded at intervals or at convergence.}
#'   \item{obj_min}{The minimum objective loss achieved during the iterations.}
#'
#' @details
#' The algorithm alternates between estimating the missing entries of the matrix \code{X} and performing singular value
#' decomposition (SVD) on the imputed matrix. Singular values are soft-thresholded using the regularization parameter \code{lambda}.
#' Convergence is determined by the change in the objective loss between iterations or reaching the maximum number of iterations.
#'
#' The function relies on the \code{rank_restricted_soft_svd} function to perform SVD with a rank restriction.
#'
#' @examples
#' # Example usage:
#' set.seed(123)
#' X <- matrix(rnorm(100), nrow = 10)
#' X[sample(1:100, 20)] <- NA # Introduce missing values
#' result <- soft_impute(X, r = 3, lambda = 0.1)
#' print(result$M_hat) # Imputed matrix
#' plot(result$times, result$objectives, type = "o", main = "Objective Loss vs Time")
#'
#' @export
soft_impute <- function(X, r, lambda, max_iter = 100, max_iter_svd = 100, tol = 1e-6, k = 1) {
  m <- nrow(X)
  n <- ncol(X)

  # Mask of observed entries
  M <- !is.na(X)

  # Step 1: Initialize M_hat
  M_hat <- X
  M_hat[!M] <- 0
  obj_old <- Inf

  times <- numeric()
  objectives <- numeric()
  cumulative_time <- 0

  for (iter in 1:max_iter) {
    start_time <- Sys.time()

    # Step 2: Update X_hat and compute SVD
    R_obs <- matrix(0, m, n)
    R_obs[M] <- (X - M_hat)[M]
    X_hat <- M_hat + R_obs

    svd_res <- rank_restricted_soft_svd(X_hat, r, max_iter = max_iter_svd, lambda = 0, verbose = FALSE)
    U <- svd_res$U
    D <- svd_res$D
    V <- svd_res$V

    # Apply soft-thresholding to singular values
    D_soft <- pmax(diag(D) - lambda, 0)
    D_soft_mat <- diag(D_soft)

    # Update M_hat
    M_hat <- U %*% D_soft_mat %*% t(V)

    # Compute objective function
    obj_new <- 0.5 * sum(((X - M_hat)[M])^2) + lambda * sum(D_soft)

    # Record cumulative time and objective loss every k steps or at convergence
    end_time <- Sys.time()
    cumulative_time <- cumulative_time + as.numeric(difftime(end_time, start_time, units = "secs"))
    if (iter %% k == 0 || iter == max_iter || abs(obj_old - obj_new)  < tol) {
      times <- c(times, cumulative_time)
      objectives <- c(objectives, obj_new)
    }

    if (abs(obj_old - obj_new) < tol) {
      cat("Soft_Impute Converged at iteration:", iter, "with objective loss:", obj_new, "\n")
      break
    }

    if (iter == max_iter) {
      cat("Soft_Impute Reached maximum iterations:", max_iter, "with objective loss:", obj_new, "\n")
    }

    obj_old <- obj_new
  }

  return(list(M_hat = M_hat, times = times, objectives = objectives, obj_min = min(objectives)))
}








