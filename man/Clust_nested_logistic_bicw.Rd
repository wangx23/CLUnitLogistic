\name{Clust_nested_logistic_bicw}
\alias{Clust_nested_logistic_bicw}
\alias{Clust_nested_logistic}
\title{Clustered Nested-Error Logistic Regression over a Penalty Path}
\description{
Fits the unit-level logistic mixed model with area-specific random effects and
clustered regression coefficients over a grid of fusion-penalty values.
\code{Clust_nested_logistic} fits the model for a single penalty value;
\code{Clust_nested_logistic_bicw} runs the full warm-started path and returns
the fit and a summary data frame for every value in \code{lamvec}.
}
\usage{
Clust_nested_logistic_bicw(indexy, y, x, index_d, lamvec,
  initial0 = NULL, vv0 = rep(0, length(unique(indexy))),
  var_v0 = max(0, var(vv0)), nu = 1, gam = 3, lam0 = 0.01, k0 = 10,
  maxiter = 500, inner_loop = 1, tol = 1e-3)

Clust_nested_logistic(indexy, y, x, index_d, lam, initial0,
  vv0 = rep(0, length(unique(indexy))), var_v0 = max(0.001, var(vv0)),
  nu = 1, gam = 3, lam0 = 0.01, maxiter = 500, inner_loop = 1, tol = 1e-3)
}
\arguments{
\item{indexy}{Vector of area (domain) indices, one entry per sampled unit.}
\item{y}{Binary response vector (0/1).}
\item{x}{Design matrix of covariates (include a column of 1s for an intercept).}
\item{index_d}{Integer vector giving the columns of \code{x} whose coefficients
are area-specific and subject to clustering; the remaining columns are treated
as common (shared) effects.}
\item{lamvec}{Numeric vector of fusion-penalty values (the regularization path).}
\item{lam}{A single fusion-penalty value (for \code{Clust_nested_logistic}).}
\item{initial0}{Optional list of starting values (\code{betam}, and \code{eta}
when there are common effects). If \code{NULL}, starting values are computed
with an EM + PQL initialization.}
\item{vv0}{Starting values for the area random effects.}
\item{var_v0}{Starting value for the random-effect variance.}
\item{nu}{ADMM augmented-Lagrangian parameter.}
\item{gam}{MCP concavity parameter.}
\item{lam0}{Ridge parameter used in the initialization step.}
\item{k0}{Number of clusters used to seed the initialization.}
\item{maxiter}{Maximum number of outer iterations.}
\item{inner_loop}{Number of inner coordinate updates per outer iteration.}
\item{tol}{Convergence tolerance.}
}
\value{
\code{Clust_nested_logistic_bicw} returns a list with components
\item{out}{A list of per-lambda fits, each the output of
\code{Clust_nested_logistic} (containing \code{betam}, \code{betac},
\code{rand}, \code{var_v}, \code{cluster}, \code{mu}, \code{nll}, etc.).}
\item{df}{Data frame with columns \code{lambda}, \code{Khat} (estimated number
of clusters), and the log-likelihood values \code{nll} and \code{nll1}.}
\item{initial00}{The initialization used for the path.}
}
\details{
The path is warm-started: the clustered fit at one penalty value seeds the next.
When the estimated number of clusters collapses to one, the path stops early.
}
\seealso{\code{\link{Est_area_clust}}, \code{\link{datexample}}}
\examples{
data(datexample)

lamvec <- exp(seq(-2, 0.5, length.out = 30)) ## lambda grid

\donttest{
res <- Clust_nested_logistic_bicw(indexy, y = y, x = x,
                                  index_d = 1:2, lamvec = lamvec)
res$df                      # estimated number of clusters along the path
table(res$out[[1]]$cluster) # cluster membership at the first lambda
}
}
