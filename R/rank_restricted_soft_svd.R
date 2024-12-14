#' Rank-Restricted Soft SVD
#'
#' Implements the Rank-Restricted Soft Singular Value Decomposition (SVD) algorithm for low-rank approximation
#' with regularization, based on Algorithm 2.1 in the referenced paper.
#'
#' @param X A numeric matrix of size \code{m x n}.
#' @param r An integer specifying the target rank for the low-rank approximation.
#' @param lambda A numeric value for the regularization parameter, controlling the shrinkage applied to the factor matrices.
#' @param max_iter An integer specifying the maximum number of iterations for the algorithm. Default is 100.
#' @param tol A numeric value for the convergence tolerance. The algorithm stops if the change in the objective
#'   loss between iterations is less than \code{tol}. Default is \code{1e-6}.
#' @param verbose A logical flag indicating whether to print convergence messages during the iterations. Default is \code{TRUE}.
#'
#' @return A list containing the following elements:
#'   \item{U}{The left singular matrix of size \code{m x r}.}
#'   \item{D}{The diagonal matrix of shrunk singular values of size \code{r x r}.}
#'   \item{V}{The right singular matrix of size \code{n x r}.}
#'
#' @details
#' The algorithm alternates between solving for the left and right singular vectors \code{U} and \code{V}, and the
#' singular values \code{D}, while applying a soft-thresholding operation to regularize the singular values.
#' The rank restriction is imposed to ensure that only the top \code{r} singular values and corresponding singular
#' vectors are used in the approximation.
#'
#' Convergence is determined by monitoring the change in the objective loss, which includes both the reconstruction
#' error and the regularization term.
#'
#' @references
#' Reference the original paper or technical documentation for Algorithm 2.1, where this method is described in detail.
#'
#' @examples
#' # Example usage:
#' set.seed(123)
#' X <- matrix(rnorm(100 * 50), 100, 50)
#' result <- rank_restricted_soft_svd(X, r = 5, lambda = 0.1)
#' U <- result$U
#' D <- result$D
#' V <- result$V
#' print(dim(U)) # Dimensions of U
#' print(dim(D)) # Dimensions of D
#' print(dim(V)) # Dimensions of V
#'
#' @export
rank_restricted_soft_svd<- function(X, r, lambda, max_iter = 100, tol = 1e-6, verbose = TRUE) {

  m <- nrow(X)
  n <- ncol(X)

  # Step 1: Random Initialize A = U D.
  U_0 <- qr.Q(qr(matrix(rnorm(m * r), m, r)))
  D_0 <- diag(1, r, r)
  A <- U_0 %*% D_0

  obj_old <- Inf
  for (iter in 1:max_iter) {
    # Step 2: Given A, solve for B
    if (iter == 1) {
      U <- U_0
      D <- D_0
    } else {

    }
    D2_lI_inv <- solve(D^2 + lambda * diag(r))
    B_t <- D2_lI_inv %*% D %*% t(U) %*% X  # size r x n

    # Step 3: Update V, D and B
    B=t(B_t)
    B_D <- B %*% D
    svd_BD <- svd(B_D)
    V_tilde <- svd_BD$u
    D_tilde <- diag(sqrt(svd_BD$d), r, r)
    V <- V_tilde
    D <- D_tilde
    B <- V %*% D

    # Step 4: Given B, solve for A

    D2_lI_inv <- solve(D^2 + lambda * diag(r))
    A_tilde <- X %*% V %*% D %*% D2_lI_inv

    # Step 5: Update U, D, A
    A_tildeD <- A_tilde %*% D
    svd_AtD <- svd(A_tildeD)
    U_tilde <- svd_AtD$u
    D_tilde <- diag(sqrt(svd_AtD$d), r, r)
    U <- U_tilde
    D <- D_tilde
    A <- U %*% D

    # Step 6: Check convergence:
    ABt <- A %*% t(B)
    obj_new <- 0.5 * sum(((X - ABt))^2) + (lambda / 2) * (sum(A^2) + sum(B^2))
    if (iter>=2 & abs(obj_old - obj_new)  < tol) {
      if (verbose) {
        cat("Converged at iteration:", iter,"with  objective loss:",  obj_new , "\n")
      }
      break
    }


    # If maximum iterations reached, print message
    if (iter == max_iter) {
      if (verbose) {
        cat("Reached maximum iterations:", max_iter,"with objective loss:", obj_new, "\n")
      }
    }

    obj_old <- obj_new
  }

  # Step 7: Compute M = X V and then SVD, apply soft thresholding
  M <- X %*% V
  svd_M <- svd(M)
  U_M <- svd_M$u
  D_M <- svd_M$d
  R_M <- svd_M$v
  D_soft <- pmax(D_M - lambda, 0)
  U_final <- U_M
  V_final <- V %*% R_M
  D_final <- diag(D_soft, length(D_soft), length(D_soft))

  return(list(U = U_final[, 1:r, drop = FALSE],
              D = D_final[1:r, 1:r, drop = FALSE],
              V = V_final[, 1:r, drop = FALSE]))
}


