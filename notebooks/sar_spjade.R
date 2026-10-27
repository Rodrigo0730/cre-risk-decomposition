library(sf)
library(spdep)
library(spatialreg)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dbscan)
library(Matrix)
library(RSpectra)
library(moments)
library(knitr)
library(glasso)
library(igraph)
library(SpatialBSS)
library(JADE)
source("src/spICA.R")

data <- st_read("data/processed/londonCRE.gpkg", quiet = TRUE)

set.seed(123) # optional, ensures you get the same sample each run
data <- data[sample(nrow(data), 5000), ]

cre_27700 <- st_transform(data, crs = 27700)
coords <- st_coordinates(cre_27700)
n <- nrow(coords)

### Create sparse spatial weights for linear spatial model fitting later
knn_list <- knearneigh(coords, k = 10)
nb_knn <- knn2nb(knn_list)
listw_sparse <- nb2listw(nb_knn, style = "W")

conditional_log_transform <- function(x, var_name, kurtosis_threshold = 10) {
  # Calculate the kurtosis of the variable
  k_val <- moments::kurtosis(x, na.rm = TRUE)

  # 1. Test if transformation is needed (Extreme Kurtosis)
  if (k_val > kurtosis_threshold) {
    cat(sprintf(
      "Variable '%s' has extreme kurtosis (%.2f). Applying log transformation.\n",
      var_name,
      k_val
    ))

    # 2. Test for zeros or negative values to prevent -Inf from log()
    min_val <- min(x, na.rm = TRUE)
    if (min_val <= 0) {
      # Shift the data so the minimum value becomes strictly positive (1) before log
      shift_constant <- abs(min_val) + 1
      cat(sprintf(
        " -> Contains zeros or negatives. Shifting by %.2f before log.\n",
        shift_constant
      ))
      return(log(x + shift_constant))
    } else {
      # Data is strictly positive, safe to apply standard log
      return(log(x))
    }
  } else {
    # Kurtosis is within normal bounds, no transformation needed
    cat(sprintf(
      "Variable '%s' has stable kurtosis (%.2f). No transformation applied.\n",
      var_name,
      k_val
    ))
    return(x)
  }
}

# Apply the test and conditional transformation to all target variables
cre_27700$Trans_Value <- conditional_log_transform(cre_27700$Value, "Value")
cre_27700$Trans_Volatility <- conditional_log_transform(
  cre_27700$Volatility,
  "Volatility"
)
cre_27700$Trans_Liquidity <- conditional_log_transform(
  cre_27700$Liquidity,
  "Liquidity"
)

# Standardize the processed data for BSS extraction
X <- scale(cbind(
  cre_27700$Trans_Value,
  cre_27700$Trans_Volatility,
  cre_27700$Trans_Liquidity
))
colnames(X) <- c("Value", "Volatility", "Liquidity")

knn_econ <- dbscan::kNN(X, k = 5)
k_neighbors <- ncol(knn_econ$id)

i <- rep(1:n, times = k_neighbors)
j <- as.vector(knn_econ$id)

Adjacency <- sparseMatrix(
  i = i,
  j = j,
  x = 1,
  dims = c(n, n),
  symmetric = FALSE
)

# 3. Make the matrix symmetric (if Property A is like B, B is like A)
Adjacency <- Adjacency + t(Adjacency)
Adjacency[Adjacency > 0] <- 1
diag(Adjacency) <- 0

# 4. Generate multiple network lags to avoid eigenvalue degeneracy in spJADE
g <- graph_from_adjacency_matrix(Adjacency, mode = "undirected")
dist_mat <- distances(g)

data_driven_kernels_raw <- list(
  as(dist_mat == 1, "dgCMatrix") * 1, # Lag 1: Direct economic neighbors
  as(dist_mat == 2, "dgCMatrix") * 1, # Lag 2: 2-hop economic neighbors
  as(dist_mat == 3, "dgCMatrix") * 1 # Lag 3: 3-hop economic neighbors
)

spectral_normalize_sparse <- function(kernels) {
  lapply(kernels, function(kmat) {
    max_eigen <- RSpectra::eigs_sym(kmat, k = 1, which = "LA")$values
    if (max_eigen == 0) {
      max_eigen <- 1
    }
    return(kmat / max_eigen)
  })
}


data_driven_kernel_list <- spectral_normalize_sparse(data_driven_kernels_raw)
Sf_list_data <- lapply(data_driven_kernel_list, function(kmat) {
  sum(kmat) / nrow(kmat)
})


phi_values <- c(100, 250, 500, 1000)

kernel_list_raw <- lapply(phi_values, function(phi) {
  cutoff_dist <- 5 * phi
  nn <- dbscan::frNN(coords, eps = cutoff_dist)

  i <- rep(1:n, lengths(nn$id))
  j <- unlist(nn$id)
  d <- unlist(nn$dist)

  weights <- exp(-d / phi)
  W_sparse <- sparseMatrix(
    i = i,
    j = j,
    x = weights,
    dims = c(n, n),
    symmetric = FALSE
  )
  diag(W_sparse) <- 1
  return(as(W_sparse, "dgCMatrix"))
})


kernel_list <- spectral_normalize_sparse(kernel_list_raw)
kernel_list <- c(list(Matrix::Diagonal(n)), kernel_list)

Sf_list <- lapply(kernel_list, function(kmat) sum(kmat) / nrow(kmat))
kernel_list_dense <- lapply(kernel_list[-1], as.matrix)

# 2. Extract Independent Components from the raw standardized data
fobi_res <- FOBI(X)
jade_res <- JADE(X)
sbss_res <- sbss(x = X, kernel_list = data_driven_kernel_list)
spjade_res <- spjade(
  x = X,
  kernel_list = data_driven_kernel_list,
  Sf_list = Sf_list_data
)
spfobi_res <- spfobi(x = X, kernels = data_driven_kernel_list)

# Bind extracted components back to the spatial dataframe
cre_27700 <- cre_27700 %>%
  mutate(
    spjade_IC1 = spjade_res$s[, 1],
    spjade_IC2 = spjade_res$s[, 2],
    spjade_IC3 = spjade_res$s[, 3],
    spfobi_IC1 = spfobi_res$s[, 1],
    spfobi_IC2 = spfobi_res$s[, 2],
    spfobi_IC3 = spfobi_res$s[, 3],
    jade_IC1 = jade_res$S[, 1],
    jade_IC2 = jade_res$S[, 2],
    jade_IC3 = jade_res$S[, 3],
    fobi_IC1 = fobi_res$S[, 1],
    fobi_IC2 = fobi_res$S[, 2],
    fobi_IC3 = fobi_res$S[, 3]
  )

knn_test <- knearneigh(coords, k = 15)
nb_test <- knn2nb(knn_test)
listw_test <- nb2listw(nb_test, style = "W", zero.policy = TRUE)

get_spatial_q_stat <- function(data, prefix, listw_fit, listw_test) {
  df_clean <- st_drop_geometry(data)
  cols <- grep(
    paste0("^", prefix, "_IC"),
    names(df_clean),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(cols) == 0) {
    return(rep(NA, 3))
  }

  q_stats <- sapply(cols, function(col_name) {
    x <- df_clean[[col_name]]

    cat("\n--- Processing:", col_name, "---\n")
    cat("Mean of x:", mean(x), "\n")
    cat("Var of x:", var(x), "\n")

    # Step A: Filter out linear spatial autocorrelation
    fit <- spatialreg::errorsarlm(x ~ 1, listw = listw_fit, method = "LU")
    cat("GMerrorsar lambda:", fit$lambda, "\n")

    # Step B & C: Extract residuals and square them
    resids <- as.numeric(residuals(fit))
    cat("Residual mean:", mean(resids), "\n")
    cat("Residual variance:", var(resids), "\n")

    resids_sq <- resids^2
    cat("Squared residuals mean:", mean(resids_sq), "\n")
    cat("Squared residuals variance:", var(resids_sq), "\n")

    # Step D: Compute Moran's I
    moran_test <- spdep::moran.test(
      resids_sq,
      listw_test,
      randomisation = FALSE,
      zero.policy = TRUE
    )

    cat("Full moran.test estimate:\n")
    print(moran_test$estimate)
    cat("Names of estimate:", names(moran_test$estimate), "\n")

    # FIXED: element [1] is the Moran's I statistic
    return(as.numeric(moran_test$estimate[1]))
  })

  ordered_q_stats <- sort(q_stats, decreasing = TRUE)
  return(ordered_q_stats[1:3])
}
prefixes <- c("spjade", "spfobi", "jade", "fobi")
method_names <- c("spJADE", "spFOBI", "JADE", "FOBI")

# Compute the rigorous Spatial Q-statistic for all methods
sq_morans_list <- lapply(prefixes, function(p) {
  get_spatial_q_stat(
    cre_27700,
    p,
    listw_fit = listw_sparse,
    listw_test = listw_test
  )
})
sq_morans_matrix <- do.call(rbind, sq_morans_list)

# Construct the final evaluation table
results_squared_morans <- data.frame(
  Methodology = method_names,
  IC1_Volatility = sq_morans_matrix[, 1],
  IC2_Volatility = sq_morans_matrix[, 2],
  IC3_Volatility = sq_morans_matrix[, 3]
)

# Output as a clean Markdown table
kable(
  results_squared_morans,
  format = "markdown",
  digits = 4,
  col.names = c(
    "Methodology",
    "Component 1 (I)",
    "Component 2 (I)",
    "Component 3 (I)"
  ),
  caption = "Spatial autocorrelation of SQUARED RESIDUALS (Volatility Clustering). Higher values indicate the algorithm successfully grouped structural shocks geographically."
)
