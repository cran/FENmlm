% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MiscFuns.R
\name{logLik.femlm}
\alias{logLik.femlm}
\title{Extracts the log-likelihood}
\usage{
\method{logLik}{femlm}(object, ...)
}
\arguments{
\item{object}{An object of class \code{femlm}. Typically the result of a \code{\link[FENmlm]{femlm}} estimation.}

\item{...}{Not currently used.}
}
\value{
It returns a numeric scalar.
}
\description{
This function extracts the log-likelihood from a \code{\link[FENmlm]{femlm}} estimation.
}
\details{
This function extracts the log-likelihood based on the model fit. You can have more information on the likelihoods in the details of the function \code{\link[FENmlm]{femlm}}.
}
\examples{

# simple estimation on iris data, clustering by "Species"
res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
            Petal.Width | Species, iris)

nobs(res)
logLik(res)


}
\seealso{
\code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{AIC.femlm}}, \code{\link[FENmlm]{BIC.femlm}}, \code{\link[FENmlm]{nobs.femlm}}.
}
\author{
Laurent Berge
}
