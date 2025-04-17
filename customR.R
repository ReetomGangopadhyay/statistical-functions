mymvr = function(Xk, Yk, zx=FALSE, ridged = 0) {
  
  n <- nrow(Yk)
  d <- ncol(Yk)
  p <- ncol(Xk)
  Y <- Yk
  
  if(is.null(p)) { # for case of simple linear regression
    p <- 1
  }
  
  # standardizes predictors if needed
  if(zx == TRUE) {
    Xk <- scale(Xk)
  }
  
  X <- cbind(1, Xk)
  
  # pseudo-determinant of design matrix
  dx <- sqrt(abs(det(t(X) %*% X)))
  
  xtxi = solve(t(X) %*% X + ridged^2 * diag(p+1))
  ##xtxi <- solve(t(X) %*% X + ridged * diag(p + 1))
  
  #betahat <- solve(t(X) %*% X + ridged^2 %*% diag(p+1)) %*% t(X) %*% Y
  betahat <- solve(t(X) %*% X + ridged^2 * diag(p+1)) %*% t(X) %*% Y
  
  Yhat <- X %*% betahat
  resid <- Y - Yhat
  
  # define hat matrix
  H <- X %*% solve(t(X) %*% X + ridged^2 * diag(p+1)) %*% t(X)
  lev <- diag(H)
  
  R <- cor(Xk)
  # variance inflation factors
  # vif = 1 / (1- R^2) r2 is correlation between xj and the other xks
  # mci = 'multicollinearity index' = sqrt(vif)
  vif <- diag(solve(R))
  mci <- sqrt(vif) # close to 1 good, bigger than 3* bad.
  
  if (d == 1) {
    SST <- var(Y) * (n - 1) # If Y is a matrix then variance of Y will not be a number
    SSE <- sum(resid^2)
    SSM <- sum((Yhat - mean(Y))^2)
    
    MST <- SST / (n - 1)
    MSE <- SSE / (n - p - 1)
    MSM <- SSM / (p)
    sig2 <- MSE
    
    Fstat <- MSM / MSE
    pval <- pf(Fstat, p, n - p - 1, lower.tail = FALSE)
    
    # var(eps) = sig^2 * (I - H)
    sresid <- resid / (sqrt(sig2) * (sqrt(1 - lev)))
    
    # var(betahat) = sig^2 * inv(t(X) %*% X)
    sebhat <- sqrt(sig2 * solve(t(X) %*% X))
    
    r2 <- 1 - SSE / SST
    r2adj <- 1 - MSE / MST
  }
  
  if (d > 1) {
    SST <- t(Y) %*% (diag(n) - (1/n) * matrix(1, n, n)) %*% Y
    SSE <- t(Y) %*% (diag(n) - H) %*% Y
    # SSM <- t(Y) %*% (H - (1/n) * matrix(1, n, n)) %*% Y
    SSM <- SST - SSE
    
    sig2 <- SSE / (n - p - 1)
    
    MST <- SST / (d * (n - 1))
    MSE <- SSE / (d * (n - p - 1))
    MSM <- SSM / (d * p)
    
    r2 <- c(1 - (sum(diag(SSE)) / sum(diag(SST))), 1 - (det(SSE) / det(SST)))
    r2adj <- c(1 - (sum(diag(MSE)) / sum(diag(MST))), 1 - (det(MSE) / det(MST)))
    
    sebhat <- matrix(sqrt(diag(kronecker(sig2, xtxi))), p+1, d)
    seeps <- matrix(sqrt(diag(kronecker(sig2, diag(n) - H))), n, d)
    
    sresid <- resid / seeps
    
    # pillai roy lowry wilks modified wilks
    
    g <- p + 1
    Fstat_pi <- sum(diag(solve(SST) %*% SSM))
    Fstat_roy <- max(eigen(solve(SSE) %*% SSM)$values)
    Fstat_lowry <- sum(diag(solve(SSE) %*% SSM))
    Fstat_wilks <- (n - g - 0.5 * (d + 2)) * log(det(SST) / det(SSE))
    Fstat_modwilks <- ((n - g - d + 1) / (d * (p - 1))) * log(det(SST) / det(SSE))
    
    Fstat <- c(Fstat_pi, Fstat_roy, Fstat_lowry, Fstat_wilks, Fstat_modwilks)
    
    pval_pi <- pf(Fstat_pi, d * g, n - g - 2, lower.tail = FALSE) 
    pval_roy <- pf(Fstat_roy, d, n - g - d - 1, lower.tail = FALSE)
    pval_lowry <- pf(Fstat_lowry, d * (g - 1), n - g - d + 1, lower.tail = FALSE)
    pval_wilks <- pchisq(Fstat_wilks, d * (g - 1), lower.tail = FALSE)
    pval_modwilks <- pf(Fstat_modwilks, d * (g - 1), n - g - d + 1, lower.tail = FALSE)
    
    pval <- c(pval_pi, pval_roy, pval_lowry, pval_wilks, pval_modwilks)
  }
  
  
  results <- list('detx' = dx,
                  'betahat' = betahat,
                  'sebetahat' = sebhat,
                  'resid' = resid,
                  'sresid' = sresid,
                  'r2' = r2,
                  'r2adj' = r2adj,
                  'Fstat' = Fstat,
                  'pval' = pval,
                  'mci' = mci)
  
  return(results)
}



# Function to calculate Shannon Entropy and Gini Impurity
calculate_metrics <- function(variable, base = exp(1)) {
  # Get proportions of each class
  proportions <- prop.table(table(variable))
  
  # Shannon Entropy
  entropy <- -sum(proportions * log(proportions, base = base))
  
  # Gini Impurity
  gini <- 1 - sum(proportions^2)
  
  return(list(Entropy = entropy, Gini = gini))
}



# Mahalanobis distance function
mymaha <- function(X, mu, sighat) {
  n = dim(X)[1]  # Number of observations
  d = dim(X)[2]  # Number of dimensions
  
  dm = c()  # To store distances
  siginv = solve(sighat)  # Inverse of covariance matrix
  for (k in 1:n) {
    myx = X[k, ]
    dm = c(dm, sqrt((myx - mu) %*% siginv %*% (myx - mu)))
  }
  return(dm)
}

# Henze-Zirkler test statistic function
hz_statistic <- function(dm, eta_1, eta_2) {
  n <- length(dm)
  
  # Compute w_k = (1/2) * exp(-0.5 * D_i_mean)
  w_k <- (1 / 2) * exp(-0.5 * dm^2)
  
  # Henze-Zirkler statistic
  hz <- n + 2 * eta_1 * sum(w_k) + eta_2 * (sum(w_k^2) - n)
  return(hz)
}


# Function to calculate HZ statistic for a given sample
calculate_hz <- function(sample_data, eta1, eta2) {
  # Ensure sample_data is a numeric matrix
  sample_data <- as.matrix(sample_data)
  
  n <- nrow(sample_data)
  d <- ncol(sample_data)
  
  mu <- colMeans(sample_data)
  Sigma <- cov(sample_data)
  
  # Compute Mahalanobis distances using mymaha
  D_i_mean <- (mymaha(sample_data, mu, Sigma))^2
  
  # Compute w_k = (1/2) * exp(-0.5 * D_i_mean)
  w_k <- (1 / 2) * exp(-0.5 * D_i_mean)
  
  # Compute the HZ test statistic
  hz_statistic <- n + 2 * eta1 * sum(w_k) + eta2 * (sum(w_k^2) - n)
  
  return(round(hz_statistic, 4))
}



initialize_centroids <- function(k, min_vals, max_vals) {
  centroids <- matrix(NA, nrow = k, ncol = ncol(X_scaled))
  for (i in 1:ncol(X_scaled)) {
    centroids[, i] <- runif(k, min_vals[i], max_vals[i])
  }
  return(centroids)
}

lightning_macqueen <- function(X, k, max_iter = 100) {
  # Initialize centroids w/ my goofy function
  centroids <- initialize_centroids(k, min_vals, max_vals)
  
  # Initialize cluster assignments
  clusters <- rep(0, nrow(X))
  
  for (iter in 1:max_iter) {
    changed <- FALSE
    
    for (i in 1:nrow(X)) {
      # Assign each point to the nearest centroid
      distances <- apply(centroids, 1, function(c) sqrt(sum((X[i, ] - c)^2))) # Euclidean distance
      new_cluster <- which.min(distances)
      
      # Update cluster assignment
      if (clusters[i] != new_cluster) {
        clusters[i] <- new_cluster
        changed <- TRUE
      }
      
      # Update centroid one by one as macqueen wants
      cluster_points <- X[clusters == new_cluster, , drop = FALSE]
      centroids[new_cluster, ] <- colMeans(cluster_points)
    }
    
    # kaboom if no change
    if (!changed) break
  }
  
  return(list(clusters = clusters, centroids = centroids))
}

tow_mater <- function(clusters, k) { # counts individuals in each cluster
  sapply(1:k, function(cluster) sum(clusters == cluster))
}

## NOT SURE IF YOU WANT MY DOO-HICKEY FUNCTION HER,FOR Q2 its not general.
evil_perc <- function(K, distances, labels) {
  # Find the indices of the K nearest neighbors
  mykranks <- which(rank(distances, ties.method = "first") <= K)
  
  nearest_labels <- labels[mykranks]
  
  # Identify the most frequent department
  most_frequent_dept <- names(sort(table(nearest_labels), decreasing = TRUE))[1]
  
  # Count the percentage of neighbors in the most frequent department
  most_frequent_count <- sum(nearest_labels == most_frequent_dept)
  percentage <- (most_frequent_count / K) * 100
  
  return(list(percentage = percentage, department = most_frequent_dept))
}

natural_entropy <- function(data_column) {
  prob <- table(data_column) / length(data_column)
  -sum(prob * log(prob), na.rm = TRUE)
}
