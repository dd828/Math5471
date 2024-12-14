#' Alternating Least Squares (ALS) Matrix Completion
#'
#' Implements a modified Alternating Least Squares (ALS) algorithm for matrix completion,
#' based on Algorithm 5.2 in the referenced paper, with the addition of a regularization parameter (\code{lambda}).
#'
#' @param X A numeric matrix of size \code{m x n} with missing values (use \code{NA} for missing entries).
#' @param r An integer specifying the target rank for the low-rank approximation.
#' @param lambda A numeric value for the regularization parameter, controlling the ridge penalty applied to the factor matrices.
#' @param max_iter An integer specifying the maximum number of ALS iterations. Default is 100.
#' @param tol A numeric value for the convergence tolerance. The algorithm stops if the change in the objective
#'   loss between iterations is less than \code{tol}. Default is \code{1e-6}.
#' @param k An integer specifying the interval for recording cumulative time and objective loss.
#'   Metrics are recorded every \code{k} steps or at convergence. Default is 1.
#'
#' @return A list containing the following elements:
#'   \item{A}{The factor matrix of size \code{m x r} corresponding to the rows of the input matrix.}
#'   \item{B}{The factor matrix of size \code{n x r} corresponding to the columns of the input matrix.}
#'   \item{times}{A numeric vector of cumulative times (in seconds) recorded at intervals or at convergence.}
#'   \item{objectives}{A numeric vector of objective loss values recorded at intervals or at convergence.}
#'   \item{obj_min}{The minimum objective loss achieved during the iterations.}
#'
#'
#'
#' @examples
#' # Example usage:
#' set.seed(123)
#' X <- matrix(rnorm(100), nrow = 10)
#' X[sample(1:100, 20)] <- NA # Introduce missing values
#' result <- als_matrix_completion(X, r = 3, lambda = 0.1)
#' print(result$A) # Row factors
#' print(result$B) # Column factors
#' plot(result$times, result$objectives, type = "o", main = "Objective Loss vs Time")
#'
#' @export
als_matrix_completion<- function(X, r, lambda, max_iter = 100, tol = 1e-6, k = 1) {
  m <- nrow(X)
  n <- ncol(X)

  # Create a mask for observed entries
  M <- !is.na(X)

  # Step 1: Initialize A and B
  M_hat <- X
  M_hat[!M] <- 0
  U=svd(M_hat)$u[,1:r]
  D=diag(sqrt(svd(M_hat)$d[1:r]))
  V=svd(M_hat)$v[,1:r]
  A <- U %*% D
  B <- V %*% D


  I_r <- diag(lambda, r, r)  # ridge terms
  obj_old <- Inf

  times <- numeric()
  objectives <- numeric()
  cumulative_time <- 0

  for (iter in 1:max_iter) {
    start_time <- Sys.time()

  # Step 2: Update A, row by row
    for (i in 1:m) {
      obs_i <- which(M[i, ])
      if (length(obs_i) > 0) {
        B_obs <- B[obs_i, , drop = FALSE]
        BtB <- t(B_obs) %*% B_obs
        rhs <- t(B_obs) %*% X[i, obs_i]
        A[i, ] <- solve(BtB + I_r, rhs)
      }
    }

  # Step 3: Update B, column by column
    for (j in 1:n) {
      obs_j <- which(M[, j])
      if (length(obs_j) > 0) {
        A_obs <- A[obs_j, , drop = FALSE]
        AtA <- t(A_obs) %*% A_obs
        rhs <- t(A_obs) %*% X[obs_j, j]
        B[j, ] <- solve(AtA + I_r, rhs)
      }
    }


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
      cat("ALS Converged at iteration:", iter, "with objective loss:", obj_new, "\n")
      break
    }

    if (iter == max_iter) {
      cat("ALS Reached maximum iterations:", max_iter, "with objective loss:", obj_new, "\n")
    }

    obj_old <- obj_new
  }

  return(list(A = A, B = B, times = times, objectives = objectives, obj_min = min(objectives)))
}
