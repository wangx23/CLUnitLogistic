\name{datexample}
\docType{data}
\alias{datexample}
\alias{indexy}
\alias{y}
\alias{x}
\alias{xunit}
\title{Example Data for Clustered Unit-Level Logistic Estimation}
\description{
A simulated example data set illustrating the inputs expected by
\code{\link{Est_area_clust}} and \code{\link{Clust_nested_logistic_bicw}}.
Loading it creates four objects: the sampled data (\code{indexy}, \code{y},
\code{x}) and the population-level covariates (\code{xunit}).
}
\usage{data(datexample)}
\format{
Loading \code{datexample} creates the following objects:
\describe{
\item{\code{indexy}}{Numeric vector of length 4505 giving the area (domain)
index of each sampled unit.}
\item{\code{y}}{Numeric 0/1 response vector of length 4505.}
\item{\code{x}}{Numeric matrix with 4505 rows and 4 columns
(\code{x1}, \code{x2}, \code{x3}, \code{x4}), the sampled design matrix.}
\item{\code{xunit}}{Numeric matrix with 45072 rows and 5 columns
(\code{area}, \code{x1}, \code{x2}, \code{x3}, \code{x4}), the population-level
design matrix whose first column is the area index and whose remaining columns
match the order of \code{x}.}
}
}
\details{
In the accompanying examples, \code{index_d = 1:2} marks \code{x1} and \code{x2}
as the area-specific (clustered) coefficients, while \code{x3} and \code{x4} are
treated as common effects.
}
\examples{
data(datexample)
dim(x)
dim(xunit)
table(table(indexy))   # sample sizes per area
}
\keyword{datasets}
