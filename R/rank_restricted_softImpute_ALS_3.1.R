#' Soft-Impute-ALS Algorithm (Algorithm 3.1)
#'
#' Implements the Soft-Impute-ALS Algorithm  (Algorithm 3.1 in the referenced paper).
#' @param X A numeric matrix of size \code{m x n} with missing values (use \code{NA} for missing entries).
#' @param r An integer specifying the target rank for the low-rank approximation.
#' @param lambda A numeric value for the regularization parameter, controlling the shrinkage applied to the factor matrices.
#' @param max_iter An integer specifying the maximum number of ALS iterations. Default is 100.
#' @param tol A numeric value for the convergence tolerance. The algorithm stops if the relative change in
#'   the objective loss between iterations is less than \code{tol}. Default is \code{1e-6}.
#'
#' @return A list containing the following elements:
#'   \item{U}{The left singular matrix of size \code{m x r}.}
#'   \item{D}{The diagonal matrix of shrunk singular values of size \code{r x r}.}
#'   \item{V}{The right singular matrix of size \code{n x r}.}
#'   \item{A}{The factor matrix of size \code{m x r}, representing the row factors.}
#'   \item{B}{The factor matrix of size \code{n x r}, representing the column factors.}
#'   \item{objectives}{A numeric vector of the objective loss values at each iteration.}
#'
#' @details
#' Convergence is determined by monitoring the relative change in the objective loss, which includes the
#' reconstruction error and the regularization penalty.
#'
#'
#' @examples
#' # Example usage:
#' set.seed(123)
#' X <- matrix(rnorm(100 * 50), 100, 50)
#' X[sample(1:5000, 100)] <- NA  # Introduce missing values
#' result <- rank_restricted_softImpute_ALS_3.1(X, r = 5, lambda = 0.1)
#' U <- result$U
#' D <- result$D
#' V <- result$V
#' A <- result$A
#' B <- result$B
#' print(dim(U)) # Dimensions of U
#' print(dim(D)) # Dimensions of D
#' print(dim(V)) # Dimensions of V
#' plot(result$objectives, type = "o", main = "Objective Loss vs Iterations")
#'
#' @export
rank_restricted_softImpute_ALS_3.1<- function(X, r, lambda, max_iter = 100, tol = 1e-6) {

  # Create a mask for observed entries; true:non-misssing; False:missing
  M <- !is.na(X)
  X_filled <- X #copy of X
  X_filled[!M] <- 0  # replace NAs with 0 for convenience

  m <- nrow(X)
  n <- ncol(X)


  # Step 1: initialization of U
  M_hat <- X
  M_hat[!M] <- 0
  U=svd(M_hat)$u[,1:r]
  D=diag(sqrt(svd(M_hat)$d[1:r]))
  V=svd(M_hat)$v[,1:r]
  A <- U %*% D
  B <- V %*% D
  old_ABt <- A %*% t(B)
  obj_old <- Inf
  objectives <- numeric()
  for (iter in 1:max_iter) {
  # Step2:  Given A, solve for B

    ABt <- A %*% t(B)

    R_obs <- matrix(0, m, n)
    R_obs[M] <- (X - ABt)[M]
    X_star <- ABt + R_obs #stored as sparse plus low-rank.
    D2_lI_inv <- solve(D^2 + lambda * diag(r))
    B_t <- D2_lI_inv %*% D %*% t(U) %*% X_star
  # Step 3: Update V ,D,  B
    B=t(B_t)
    B_D <- B %*% D
    svd_BD <- svd(B_D)
    U_B <- svd_BD$u
    D_B <- diag(sqrt(svd_BD$d), r, r)
    V <- U_B
    D <- D_B
    B <- V %*% D

  # Step 4: Given B, solve for A
    ABt <- A %*% t(B)
    R_obs <- matrix(0, m, n)
    R_obs[M] <- (X - ABt)[M]
    X_star <- ABt + R_obs

    D2_lI_inv <- solve(D^2 + lambda * diag(r))
    A= X_star%*%V%*%D%*%D2_lI_inv
  # Step 5: Update U, D, A
    A_D <- A %*% D
    svd_AD <- svd(A_D)
    U_A <- svd_AD$u
    D_A <- diag(sqrt(svd_AD$d), r, r)
    U <- U_A
    D <- D_A
    A=U%*%D

  # Step 6: Check convergence using relative objective error
    ABt <- A %*% t(B)
    obj_new <- 0.5 * sum(((X - ABt)[M])^2) + (lambda / 2) * (sum(A^2) + sum(B^2))
    objectives <- c(objectives, obj_new)
    if (iter>=2 & abs(obj_old - obj_new) < tol) {
      cat("Converged at iteration:", iter,"with objective loss:",  obj_new, "\n")
      break
    }


    # If maximum iterations reached, print message
    if (iter == max_iter) {
      cat("Reached maximum iterations:", max_iter,"with objective loss:",  obj_new, "\n")
    }

    obj_old <- obj_new
  }

  # Step 7: Compute M = X V and then SVD, apply soft thresholding
  M <-  X_star %*% V
  svd_M <- svd(M)
  U_final <- svd_M$u
  D_sigma <- svd_M$d
  R_final <- svd_M$v
  D_soft <- pmax(D_sigma - lambda, 0)

  # Update U, V, D:
  U_final <- U_final
  V_final <- V %*% R_final
  D_final <- diag(D_soft, length(D_soft), length(D_soft))

  # Keep only the top r components:
  U_final <- U_final[, 1:r, drop = FALSE]
  D_final <- D_final[1:r, 1:r, drop = FALSE]
  V_final <- V_final[, 1:r, drop = FALSE]
  A_final=U_final%*% sqrt(D_final)
  B_final=V_final%*%sqrt(D_final)
  return(list(U = U_final, D = D_final, V = V_final,A=A_final,B=B_final,objectives = objectives))
}

