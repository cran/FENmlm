\name{getFE}
\alias{getFE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the Fixed-Effects from a \code{feNmlm} estimation.
}
\description{
This function retrives the fixed effects from a feNmlm estimation. It is useful only when there are more than one cluster.
}
\usage{
getFE(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A feMmlm object.
}
}

\value{
A list containig the vectors of the fixed effects.
}

\author{
Laurent Berge
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

# Bilateral network
nb = 20
n = nb**2
k = nb
id1 = factor(rep(1:k,each=n/k))
id2 = factor(rep(1:(n/k),times=k))
d = rep(rnorm(k)**2,each=n/k)
x = rnorm(n,1,5)**2
y = rnorm(n,-1,5)**2
z = rpois(n,x*y+rnorm(n,sd = 3)**2)
base = data.frame(x,y,z,id1,id2)

# We want to use the ID's of each observation as a variable: we use the option dummy
est_poisson = femlm(z~log(x)+log(y),base,family="poisson",dummy=c("id1","id2"))

# To get the FE:
myFE = getFE(est_poisson)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
