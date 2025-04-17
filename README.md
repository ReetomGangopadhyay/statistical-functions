# Custom R Functions Library

This repository provides a collection of custom R functions for multivariate regression, clustering, statistical diagnostics, and information-theoretic metrics. These utilities are designed for both educational and applied statistical modeling tasks.

---

## Function Overview

### `mymvr(Xk, Yk, zx = FALSE, ridged = 0)`

Performs multiple or multivariate linear regression with optional standardization and ridge regularization.

**Arguments**:
- `Xk`: Predictor matrix.
- `Yk`: Response vector or matrix.
- `zx`: Whether to standardize predictors (default: `FALSE`).
- `ridged`: Ridge penalty coefficient (default: `0`).

**Returns**:
- `detx`: Pseudo-determinant of the design matrix.
- `betahat`: Estimated coefficients.
- `sebetahat`: Standard errors of `betahat`.
- `resid`: Residuals.
- `sresid`: Standardized residuals.
- `r2`, `r2adj`: Coefficient of determination and adjusted R².
- `Fstat`, `pval`: F-statistics and p-values for model tests.
- `mci`: Multicollinearity indices (sqrt(VIF)).

---

### `calculate_metrics(variable, base = exp(1))`

Computes entropy and Gini impurity for a categorical variable.

**Arguments**:
- `variable`: Vector of class labels.
- `base`: Logarithmic base for entropy (default: natural log).

**Returns**:
- `Entropy`: Shannon entropy.
- `Gini`: Gini impurity.

---

### `mymaha(X, mu, sighat)`

Computes Mahalanobis distances of observations in `X` from mean `mu` using covariance matrix `sighat`.

**Arguments**:
- `X`: Data matrix.
- `mu`: Mean vector.
- `sighat`: Covariance matrix.

**Returns**:
- A vector of Mahalanobis distances.

---

### `hz_statistic(dm, eta_1, eta_2)`

Calculates the Henze-Zirkler statistic for multivariate normality.

**Arguments**:
- `dm`: Vector of Mahalanobis distances.
- `eta_1`, `eta_2`: Constants based on dimensionality.

**Returns**:
- Henze-Zirkler test statistic.

---

### `calculate_hz(sample_data, eta1, eta2)`

Wrapper for computing the Henze-Zirkler statistic directly from sample data.

**Arguments**:
- `sample_data`: Matrix or data frame of multivariate observations.
- `eta1`, `eta2`: Constants.

**Returns**:
- Rounded Henze-Zirkler statistic.

---

### `initialize_centroids(k, min_vals, max_vals)`

Generates `k` random centroids within feature bounds.

**Arguments**:
- `k`: Number of clusters.
- `min_vals`, `max_vals`: Vectors of feature-wise bounds.

**Returns**:
- A `k × d` matrix of initial centroids.

---

### `lightning_macqueen(X, k, max_iter = 100)`

Performs online (MacQueen-style) k-means clustering.

**Arguments**:
- `X`: Data matrix.
- `k`: Number of clusters.
- `max_iter`: Maximum number of iterations.

**Returns**:
- `clusters`: Final cluster assignments.
- `centroids`: Final cluster centroids.

---

### `tow_mater(clusters, k)`

Counts the number of items in each cluster.

**Arguments**:
- `clusters`: Cluster assignment vector.
- `k`: Number of clusters.

**Returns**:
- A vector of cluster counts.

---

### `evil_perc(K, distances, labels)`

Computes the percentage of K-nearest neighbors belonging to the most common label.

**Arguments**:
- `K`: Number of neighbors.
- `distances`: Distance vector.
- `labels`: Class labels of data points.

**Returns**:
- `percentage`: Percent in dominant class.
- `department`: Dominant class label.

---

### `natural_entropy(data_column)`

Computes natural (log base *e*) entropy of a categorical variable.

**Arguments**:
- `data_column`: Vector of categorical data.

**Returns**:
- Scalar entropy value.

---

## Notes

- All functions assume numeric inputs where appropriate.
- Use `mymvr()` for both univariate and multivariate linear modeling.
- Ridge penalty in `mymvr()` uses `ridged^2` in design matrix regularization.
- `lightning_macqueen()` simulates a MacQueen-style update for centroids.

---

## Use Cases

- Regression diagnostics and regularized inference.
- Evaluating multicollinearity.
- Clustering for unsupervised learning.
- Entropy-based classification metrics.
- Testing multivariate normality via HZ-statistic.

---

### `generate_tree_hierarchy(mytreeindex)`

Randomly generates a subset of variables (a "branch") from a given index set. Useful for decision tree and ensemble modeling.

**Arguments**:
- `mytreeindex`: A vector of variable indices.

**Returns**:
- A random subset of `mytreeindex` without replacement.

---

### `calculate_branch_metrics(predictor, target)`

Computes Gini impurity and information gain for a potential decision tree split.

**Arguments**:
- `predictor`: A categorical predictor variable.
- `target`: The target class variable (for which impurity is measured).

**Returns**:
- `E1`: Weighted Gini impurity after the split.
- `IG1`: Information gain (requires prior global impurity `E0` to be defined in scope).

>  Note: Assumes `E0` (entropy or Gini of parent node) is defined in the global environment.

---

### `calculate_zone_percentages(K, distances, labels)`

Computes the percentage distribution of the `K` nearest neighbors across different categorical labels (e.g., climate zones).

**Arguments**:
- `K`: Number of neighbors.
- `distances`: Vector of distances to neighbors.
- `labels`: Vector of class labels corresponding to each observation.

**Returns**:
- A named vector containing the percentage of each label among the nearest `K` neighbors.

---



## Example Usage

```r
# Fit multivariate regression
result <- mymvr(X, Y, zx = TRUE, ridged = 0.1)
print(result$betahat)

# Compute entropy and Gini impurity
calculate_metrics(c("A", "A", "B", "C", "A", "B"))

# Run clustering
res <- lightning_macqueen(X_scaled, k = 3)
table(res$clusters)
