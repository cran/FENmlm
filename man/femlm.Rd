\name{femlm}
\alias{femlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fixed effects maximum likelihood models.
}
\description{
This function estimates maximum likelihood models (e.g., Poisson or Logit) and is efficient to handle fixed effects (i.e. cluster variables). It further allows for nonlinear right hand sides.
}
\usage{
femlm(linear.fml, data, dummy, linear.start = 0,
      useHessian = TRUE, opt_method = c("nlminb", "optim"),
      debug = FALSE, family = c("poisson", "negbin", "logit"),
      opt.control=list(),optim.method="BFGS",...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{linear.fml}{
A formula. The linear formula to be estimated.
}
  \item{family}{
Character scalar. It should provide the family. Currently \code{family="poisson"}, \code{family="negbin"} and \code{family="logit"} are implemented. Note that the log link is used by default.
}
  \item{data}{
A data.frame containing the necessary variables to run the model. The variables of the non-linear right hand side of the formula are identified with this data.frame names. Note that no NA is allowed.
}
  \item{start}{
A list. Starting values for the non-linear parameters. ALL the parameters are to be named and given a staring value. Example: \code{start=list(a=1,b=5,c=0)}.
Though, there is an exception: if all parameters are to be given the same starting value, use \code{start.init}. Yet this is not recommended.
}
  \item{dummy}{
Character vector. The name/s of a/some variable/s within the dataset. These variables should contain the identifier of each observation (e.g., think of it as a panel identifier).
}
  \item{linear.start}{
Numeric named vector. The starting values of the linear part. If it is just a numeric scalar, all coefficients are set to \code{linear.start}.
}
  \item{useHessian}{
Logical. (Only if omptimization method is \code{optim}). Should the Hessian be computed in the optimization stage? Default is \code{TRUE}.
}
  \item{opt_method}{
Character scalar. Which optimization method should be used. Either \code{nlminb} or \code{optim}. Default is \code{nlminb}.
}
  \item{opt.control}{
List of elements to be passed to the optimization method (nlminb or optim).
}
  \item{optim.method}{
Character scalar. If \code{opt_method="optim"}, it is the algorithm to be used by optim (default is "BFGS"). See \code{optim} help pages for detail.
}
  \item{debug}{
Logical. If \code{TRUE} then the log-likelihood as well as all parameters are printed at each iteration. Default is \code{FALSE}.
}
  \item{\dots}{
Not currently used.
}
}

\value{
An \code{feNmlm} object.
\item{coef}{The coefficients.}
\item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
\item{loglik}{The loglikelihood.}
\item{iterations}{Number of iterations of the algorithm.}
\item{n}{The number of observations.}
\item{k}{The number of parameters of the model.}
\item{call}{The call.}
\item{nonlinear.fml}{The nonlinear formula of the call. It also contains the dependent variable.}
\item{linear.formula}{The linear formula of the call.}
\item{ll_null}{Log-likelyhood of the null model}
\item{pseudo_r2}{The adjusted pseudo R2.}
\item{naive.r2}{The R2 as if the expected predictor was the linear predictor in OLS.}
\item{message}{The convergence message from the optimization procedures.}
\item{sq.cor}{Squared correlation between the dependent variable and its expected value as given by the optimization.}
\item{expected.predictor}{The expected predictor is the expected value of the dependent variable.}
\item{cov.unscaled}{The variance covariance matrix of the parameters.}
\item{sd}{The standard error of the parameters.}
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

# Results of the Poisson..
est_poisson = femlm(z~log(x)+log(y),base,family="poisson")
# .. and of the Negative Binomial
est_negbin = femlm(z~log(x)+log(y),base,family="negbin")

# Displaying the results
est_poisson
est_negbin

# Changing the way the standard errors are computed:
summary(est_poisson,sd="white")
summary(est_negbin,sd="white")

#
# Now with dummies
#

# Bilateral network
nb = 20
n = nb**2
k = nb
id1 = factor(rep(1:k,each=n/k))
id2 = factor(rep(1:(n/k),times=k))
x = rnorm(n,1,5)**2
y = rnorm(n,-1,5)**2
z = rpois(n,x*y+rnorm(n,sd = 3)**2)
base = data.frame(x,y,z,id1,id2)

# We want to use the ID's of each observation as a variable: we use the option dummy
est_poisson = femlm(z~log(x)+log(y),base,family="poisson",dummy=c("id1","id2"))
# Displaying the results with twoway clustered santard-errors
print(est_poisson,"t")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
