\name{print.feNmlm}
\alias{print.feNmlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A print facility for \code{feNmlm} objects. It can compute different types of standard errors.
}
\description{
This function is very similar to usual \code{summary} functions as it provides the table of coefficients along with other information on the fit of the estimation.
}
\usage{
\method{print}{feNmlm}(x, sd = c("standard", "white","cluster","twoway"),cluster, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A feNmlm object.
}
  \item{sd}{
Character scalar. Which kind of standard error should be prompted: \dQuote{standard} (default), \dQuote{White}, or \dQuote{cluster}?
}
	\item{cluster}{
A list of vectors. Used only if \code{sd = "cluster"} or \code{sd="twoway"}. The vectors should give the cluster of each observation. Note that if the estimation was run using \code{dummy}, the standard error is automatically clustered along the cluster given in \code{feNmlm}.
}
  \item{\dots}{
Currently unused.
}
}

\author{
Laurent Berge
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{feNmlm}}.
}
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
					family="poisson", linear.fml=~1)

print(est0L)
print(est0NL)

#Generating a non-linear relation
z2 = rpois(n,x + y)
base$z2 = z2

#Using a non-linear form
est1L = feNmlm(z2~0,base,~log(x)+log(y),family="poi")
est1NL = feNmlm(z2~log(a*x + b*y),base,start = list(a=1,b=2),family="poisson")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
