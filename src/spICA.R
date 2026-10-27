whiten_data <- function(X, na.action = na.fail) {
  # Whitening exactly as in from JADE:::FOBI
  X <- na.action(X)
  X <- as.matrix(X)

  X_centered <- sweep(X, 2, colMeans(X), "-")

  COV <- crossprod(X_centered) / nrow(X_centered)

  EVD <- eigen(COV, symmetric = TRUE)
  U <- EVD$vectors
  Lambda <- EVD$values

  COV_inv_sqrt <- U %*% diag(1 / sqrt(Lambda)) %*% t(U)

  X_whitened <- X_centered %*% COV_inv_sqrt

  return(list(
    Y = X_whitened,
    W0 = COV_inv_sqrt
  ))
}

compute_bkl <- function(k, l, kmat, Y) {
  # For n distinct locations, computes the fourth order cross-moment
  # matrices B_f^{kl}.
  n <- nrow(Y)
  vkl <- Y[, k] * Y[, l]
  w <- as.vector(kmat %*% vkl)
  return(as.matrix(crossprod(Y * w, Y) / n))
}

compute_ckl <- function(k, l, kmat, lcov, Y, Sf) {
  n <- nrow(Y)
  p <- ncol(Y)
  Bkl <- compute_bkl(k, l, kmat, Y)

  u <- lcov[, k, drop = FALSE]
  v <- lcov[, l, drop = FALSE]
  L <- if (k == l) 2 * tcrossprod(u) else tcrossprod(u, v) + tcrossprod(v, u)

  Ckl <- if (k == l) Bkl - L / Sf - diag(p) * Sf else Bkl - L / Sf
  return(as.matrix(Ckl))
}

# --- HELPER FUNCTION FOR IRREGULAR COORDINATES ---
spectral_normalize_kernels <- function(kernels) {
  lapply(kernels, function(kmat) {
    # Calculate the largest eigenvalue (spectral radius)
    # distance/contiguity matrices are naturally symmetric
    max_eigen <- max(eigen(kmat, symmetric = TRUE, only.values = TRUE)$values)

    # Prevent division by zero if the matrix is completely empty
    if (max_eigen == 0) {
      max_eigen <- 1
    }

    return(kmat / max_eigen)
  })
}
# -------------------------------------------------

spjade <- function(x, kernel_list, Sf_list = NULL, eps = 1e-06, maxiter = 100) {
  x <- as.matrix(x)
  p <- ncol(x)
  K <- length(kernel_list)

  x_mu <- colMeans(x)

  whitening <- whiten_data(x)
  Y <- whitening$Y
  W0 <- whitening$W0

  lcovs <- lapply(kernel_list, function(W) {
    as.matrix(crossprod(Y, W %*% Y) / nrow(Y))
  })

  C_array <- array(0, dim = c(p, p, p * p * K))
  idx <- 1
  for (i in 1:K) {
    kmat <- kernel_list[[i]]
    lcov <- lcovs[[i]]
    Sf <- Sf_list[[i]]

    for (k in 1:p) {
      for (l in 1:p) {
        C_array[,, idx] <- compute_ckl(k, l, kmat, lcov, Y, Sf)
        idx <- idx + 1
      }
    }
  }

  num_matrices <- dim(C_array)[3]

  jd_res <- tryCatch(
    suppressWarnings(JADE::frjd(C_array, maxiter = maxiter, eps = eps)),
    error = function(e) {
      message("frjd failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(jd_res)) {
    stop("Joint diagonalization failed.")
  }

  W_est <- crossprod(jd_res$V, W0)

  x_centered <- sweep(x, 2, x_mu, "-")
  s <- x_centered %*% t(W_est)

  colnames(s) <- paste0("spJADE_", 1:p)

  res <- list(
    s = s,
    coords = NULL,
    w = W_est,
    x_mu = x_mu,
    d = jd_res$D,
    C_array = C_array
  )

  class(res) <- "sbss"

  return(res)
}

spfobi <- function(x, kernels, eps = 1e-06, maxiter = 100) {
  x <- as.matrix(x)

  K <- length(kernels)
  n <- nrow(x)
  p <- ncol(x)

  x_mu <- colMeans(x)

  whitening <- whiten_data(x)
  Y <- whitening$Y
  W0 <- whitening$W0

  B_array <- array(0, dim = c(p, p, K))

  for (i in 1:K) {
    kmat <- kernels[[i]]
    B_sum <- matrix(0, p, p)
    for (k in 1:p) {
      B <- compute_bkl(k, k, kmat, Y)
      B_sum <- B_sum + B
    }

    B_array[,, i] <- B_sum
  }

  num_matrices <- dim(B_array)[3]

  jd_res <- tryCatch(
    suppressWarnings(JADE::frjd(B_array, maxiter = maxiter, eps = eps)),
    error = function(e) {
      message("frjd failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(jd_res)) {
    stop("Joint diagonalization failed.")
  }

  W_est <- crossprod(jd_res$V, W0)

  x_centered <- sweep(x, 2, x_mu, "-")
  s <- x_centered %*% t(W_est)

  colnames(s) <- paste0("spFOBI_", 1:p)

  res <- list(
    s = s,
    coords = NULL,
    w = W_est,
    x_mu = x_mu,
    d = jd_res$D,
    B_array = B_array
  )

  class(res) <- "sbss"

  return(res)
}
