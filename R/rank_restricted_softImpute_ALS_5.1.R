


#' Rank-Restricted Soft-Impute ALS Algorithm (Algorithm 5.1)
#'Implement the algorithm 5.1 in the reference paper.
#' @title Matrix Completion Using Soft-Impute-ALS
#' @param X A numeric matrix with missing values (`NA`) to be imputed.
#' @param r An integer specifying the target rank for approximation.
#' @param lambda A numeric value for the regularization parameter.
#' @param max_iter An integer for the maximum number of iterations. Default is 100.
#' @param tol A numeric tolerance for convergence. Default is `1e-6`.
#' @param k An integer for recording cumulative time and objective loss every `k` steps. Default is 1.
#' @return A list containing:
#' \describe{
#'   \item{A}{The left singular factor matrix of rank `r`.}
#'   \item{B}{The right singular factor matrix of rank `r`.}
#'   \item{times}{Cumulative time at each recorded iteration.}
#'   \item{objectives}{Objective loss at each recorded iteration.}
#'   \item{obj_min}{The minimum objective loss achieved.}
#' }
#' @examples
#' # Example usage
#' X <- matrix(c(1, NA, 3, 4, 5, NA), nrow = 2)
#' result <- rank_restricted_softImpute_ALS_5.1(X, r = 2, lambda = 0.1)
#' @export
rank_restricted_softImpute_ALS_5.1 <- function(X, r, lambda, max_iter = 100, tol = 1e-6, k = 1) {
  M <- !is.na(X)
  X_filled <- X
  X_filled[!M] <- 0

  m <- nrow(X)
  n <- ncol(X)
  # Step 1:  initialization of U

  M_hat <- X
  M_hat[!M] <- 0
  U=svd(M_hat)$u[,1:r]
  D=diag(sqrt(svd(M_hat)$d[1:r]))
  V=svd(M_hat)$v[,1:r]
  A <- U %*% D
  B <- V %*% D
  old_ABt <- A %*% t(B)
  obj_old <- Inf


  times <- numeric()
  objectives <- numeric()
  cumulative_time <- 0

  for (iter in 1:max_iter) {
    start_time <- Sys.time()

  # Step 2: Update B
    ABt <- A %*% t(B)
    R_obs <- matrix(0, m, n)
    R_obs[M] <- (X - ABt)[M]
    X_star <- ABt + R_obs
    D2_lI_inv <- solve(t(A) %*% A + lambda * diag(r))
    B <- t(X_star) %*% A %*% D2_lI_inv

  # Step 3: Update A
    ABt <- A %*% t(B)
    R_obs <- matrix(0, m, n)
    R_obs[M] <- (X - ABt)[M]
    X_star <- ABt + R_obs
    D2_lI_inv <- solve(t(B) %*% B + lambda * diag(r))
    A <- X_star %*% B %*% D2_lI_inv

  # Step 4: Check convergence
    ABt <- A %*% t(B)
    obj_new <- 0.5 * sum(((X - ABt)[M])^2) + (lambda / 2) * (sum(A^2) + sum(B^2))

  # Record cumulative time and objective loss every k steps or at convergence
    end_time <- Sys.time()
    cumulative_time <- cumulative_time + as.numeric(difftime(end_time, start_time, units = "secs"))
    if (iter %% k == 0 || iter == max_iter || abs(obj_old - obj_new) < tol) {
      times <- c(times, cumulative_time)
      objectives <- c(objectives, obj_new)
    }

    if (abs(obj_old - obj_new) < tol) {
      cat("Soft-Impute-ALS Converged at iteration:", iter, "with  objective loss:", obj_new, "\n")
      break
    }

    if (iter == max_iter) {
      cat("Soft-Impute-ALS Reached maximum iterations:", max_iter, "with  objective loss:", obj_new, "\n")
    }

    obj_old <- obj_new
  }

  return(list(A = A, B = B, times = times, objectives = objectives, obj_min = min(objectives)))
}


