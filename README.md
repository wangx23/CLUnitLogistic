# CLUnitLogistic

Clustered unit-level logistic models for small area estimation.

The package fits a unit-level logistic mixed model with area-specific random
effects in which the area-level regression coefficients are grouped into an
unknown number of clusters. Clustering is estimated with an MCP fusion penalty
solved by an ADMM-type algorithm, the number of clusters is selected by BIC, and
the fitted model is combined with population covariates to build model-based
estimators of finite population proportions.

## Installation

```r
devtools::install("CLUnitLogistic")
```

## Exported functions

- `Clust_nested_logistic_bicw()` — fit the clustered logistic model over a grid
  of fusion-penalty values (warm-started path) and return the per-lambda fits
  with a summary data frame.
- `Est_area_clust()` — run the penalty path, select the number of clusters by
  BIC, and produce area-level proportion estimators from population covariates.

See `?Clust_nested_logistic_bicw` and `?Est_area_clust` for details.

## Example

```r
library(CLUnitLogistic)

lamvec <- exp(seq(-2, 0.5, length.out = 30)) ## lambda value

res_estY <- Est_area_clust(indexy, y = y, x = x, xunit, index_d = 1:2, 
                           lamvec = lamvec)
res_estY$bic #bic path
table(res_estY$out$cluster) #estimated cluster
res_estY$out$betac # estimated beta
```

`x` is the sampled design matrix (include a column of 1s for an intercept);
`xunit` is the population design matrix whose first column is the area index and
whose remaining columns match the order of `x`. `index_d` gives the columns of
`x` whose coefficients are area-specific and subject to clustering.
