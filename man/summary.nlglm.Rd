\name{summary.feNmlm}
\alias{summary.feNmlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of a \code{feNmlm} object. Computes different types of standard errors.
}
\description{
This function is similar to \code{print.feNmlm}. It provides the table of coefficients along with other information on the fit of the estimation. It can compute different types of standard errors. The new variance covariance matrix is an object returned.
}
\usage{
\method{summary}{feNmlm}(object,sd=c("standard","white","cluster","twoway"),
         cluster,dof_correction=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A feNmlm object.
}
  \item{sd}{
Character scalar. Which kind of standard error should be prompted: \dQuote{standard} (default), \dQuote{White}, or \dQuote{cluster}?
}
  \item{cluster}{
A list of vectors. Used only if \code{sd = "cluster"} or \code{sd="twoway"}. The vectors should give the cluster of each observation. Note that if the estimation was run using \code{dummy}, the standard error is automatically clustered along the cluster given in \code{feNmlm}.
}
  \item{dof_correction}{
Logical. Should a finite sample correcton be applied? (Default is \code{TRUE}.)
}
  \item{\dots}{
Currently unused.
}
}

\value{
The same values as a \code{feMmlm} object plus:
\item{vcov}{The variance-covariance matrix whose type is the one requested by the user.}
}

\author{
Laurent Berge
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{


#The data
n = 100
x = rnorm(n,1,5)**2
y = rnorm(n,-1,5)**2
z = rpois(n,x*y)
base = data.frame(x,y,z)

#Comparing the results of a 'linear' function
est0L = feNmlm(z~0,base,~log(x)+log(y),family="poi")
est0NL = feNmlm(z~a*log(x)+b*log(y),base,start = list(a=0,b=0),
					family="poisson",	linear.fml=~1)

# Displaying the summary
summary(est0L,sd="white")
myWhiteVcov = summary(est0L,sd="white")$vcov


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
