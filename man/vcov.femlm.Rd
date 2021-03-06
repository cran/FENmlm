% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MiscFuns.R
\name{vcov.femlm}
\alias{vcov.femlm}
\title{Extract the variance/covariance of a femlm fit}
\usage{
\method{vcov}{femlm}(object, se = c("standard", "white", "cluster",
  "twoway", "threeway", "fourway"), cluster, dof_correction = FALSE,
  forceCovariance = FALSE, keepBounded = FALSE, ...)
}
\arguments{
\item{object}{A femlm object. Obtained using \code{\link[FENmlm]{femlm}}.}

\item{se}{Character scalar. Which kind of standard error should be computed: \dQuote{standard} (default), \dQuote{White}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}?}

\item{cluster}{A list of vectors. Used only if \code{se="cluster"}, \dQuote{se=twoway}, \dQuote{se=threeway} or \dQuote{se=fourway}. The vectors should give the cluster of each observation. Note that if the estimation was run using \code{cluster}, the standard error is automatically clustered along the cluster given in \code{\link[FENmlm]{femlm}}. For one-way clustering, this argument can directly be a vector (instead of a list). If the estimation has been done with cluster variables, you can give a character vector of the dimensions over which to cluster the SE.}

\item{dof_correction}{Logical, default is \code{FALSE}. Should there be a degree of freedom correction to the standard errors of the coefficients?}

\item{forceCovariance}{(Advanced users.) Logical, default is \code{FALSE}. In the peculiar case where the obtained Hessian is not invertible (usually because of collinearity of some variables), use this option force the covariance matrix, by using a generalized inverse of the Hessian. This can be useful to spot where possible problems come from.}

\item{keepBounded}{(Advanced users.) Logical, default is \code{FALSE}. If \code{TRUE}, then the bounded coefficients (if any) are treated as unrestricted coefficients and their S.E. is computed (otherwise it is not).}

\item{...}{Other arguments to be passed to \code{\link[FENmlm]{summary.femlm}}.

The computation of the VCOV matrix is first done in \code{\link[FENmlm]{summary.femlm}}.}
}
\value{
It returns a \eqn{N\times N} square matrix where \eqn{N} is the number of variables of the fitted model.
This matrix has an attribute \dQuote{type} specifying how this variance/covariance matrix has been commputed (i.e. was it created using White correction, or was it clustered along a specific factor, etc).
}
\description{
This function extracts the variance-covariance of estimated parameters from a model estimated with \code{\link[FENmlm]{femlm}}.
}
\examples{

# Load trade data
data(trade)

# We estimate the effect of distance on trade (with 3 fixed-effects)
est_pois = femlm(Euros ~ log(dist_km) + log(Year) | Origin + Destination +
                 Product, trade)

# "normal" VCOV
vcov(est_pois)

# "white" VCOV
vcov(est_pois, se = "white")

# "clustered" VCOV (with respect to the Origin factor)
vcov(est_pois, se = "cluster")

# "clustered" VCOV (with respect to the Product factor)
vcov(est_pois, se = "cluster", cluster = trade$Product)
# another way to make the same request:
vcov(est_pois, se = "cluster", cluster = "Product")

# Another estimation without cluster:
est_pois_simple = femlm(Euros ~ log(dist_km) + log(Year), trade)

# We can still get the clustered VCOV,
# but we need to give the cluster-vector:
vcov(est_pois_simple, se = "cluster", cluster = trade$Product)


}
\seealso{
\code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{summary.femlm}}, \code{\link[FENmlm]{confint.femlm}}, \code{\link[FENmlm]{resid.femlm}}, \code{\link[FENmlm]{predict.femlm}}, \code{\link[FENmlm]{getFE}}.
}
\author{
Laurent Berge
}
