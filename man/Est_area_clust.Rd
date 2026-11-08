\name{Est_area_clust}
\alias{Est_area_clust}
\title{Model-Based Small Area Estimators from a Clustered Logistic Model}
\description{
Fits the clustered nested-error logistic model over a penalty path, selects the
number of clusters by a BIC criterion, and uses the selected model together with
population-level covariates to construct model-based estimators of the finite
population proportion in each area.
}
\usage{
Est_area_clust(area, y, x, xunit, index_d, lamvec,
  initial0 = NULL, vv0 = rep(0, length(unique(area))),
  var_v0 = max(0.001, var(vv0)), nu = 1, gam = 3, lam0 = 0.01,
  maxiter = 300, inner_loop = 1, tol = 1e-3)
}
\arguments{
\item{area}{Vector of area (domain) indices for the sampled units.}
\item{y}{Binary response vector (0/1) for the sampled units.}
\item{x}{Design matrix for the sampled units (include a column of 1s for an
intercept).}
\item{xunit}{Population-level design matrix whose first column is the area index
and whose remaining columns match the order of \code{x} (include the intercept
column if the model has one).}
\item{index_d}{Integer vector giving the columns of \code{x} whose coefficients
are area-specific and clustered; remaining columns are common effects.}
\item{lamvec}{Numeric vector of fusion-penalty values (the regularization path).}
\item{initial0}{Optional list of starting values; computed automatically when
\code{NULL}.}
\item{vv0}{Starting values for the area random effects.}
\item{var_v0}{Starting value for the random-effect variance.}
\item{nu}{ADMM augmented-Lagrangian parameter.}
\item{gam}{MCP concavity parameter.}
\item{lam0}{Ridge parameter used in the initialization step.}
\item{maxiter}{Maximum number of outer iterations.}
\item{inner_loop}{Number of inner coordinate updates per outer iteration.}
\item{tol}{Convergence tolerance.}
}
\value{
A list with components
\item{bic}{Data frame of the penalty path with the BIC column added.}
\item{refit}{The \code{lme4} refit at the selected cluster structure, or
\code{NULL} if the refit failed.}
\item{out}{The selected single-lambda fit.}
\item{est}{Matrix with two columns of area-level proportion estimators:
a synthetic estimator (\code{estimator1}) and a bias-corrected / survey-adjusted
estimator (\code{estimator2}), one row per area.}
\item{lam}{The penalty grid that was supplied.}
\item{initial00}{The initialization used for the path.}
}
\details{
The number of clusters is chosen as the minimizer of
\code{2 * nll + log(m * p) * log(m) * Khat * p}, where \code{m} is the number of
areas and \code{p = length(index_d)}. The selected model is then refit with
\code{lme4::glmer} on the estimated cluster structure when possible, and the
fitted coefficients and predicted random effects are combined with the
population covariates in \code{xunit} to produce the area estimators.
}
\seealso{\code{\link{Clust_nested_logistic_bicw}}, \code{\link{datexample}}}
\examples{
data(datexample)

lamvec <- exp(seq(-2, 0.5, length.out = 30)) ## lambda grid

\donttest{
res_estY <- Est_area_clust(indexy, y = y, x = x, xunit,
                           index_d = 1:2, lamvec = lamvec)
res_estY$bic              # BIC path
table(res_estY$out$cluster)  # estimated cluster membership
head(res_estY$est)        # area-level estimators
}
}
