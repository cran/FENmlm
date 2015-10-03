\name{feNmlm}
\alias{feNmlm}

\title{
Fixed effects non-linear maximum likelihood models
}

\description{
This function estimates maximum likelihood models (e.g., Poisson or Logit) and is efficient to handle fixed effects (i.e. cluster variables). It further allows for nonlinear right hand sides.
}
\usage{
feNmlm(fml,data,linear.fml,start,lower,upper,
       env,dummy,start.init,nl.gradient,linear.start=0,
       jacobian.method=c("simple","Richardson"),useHessian=TRUE,
       d.hessian,opt_method=c("nlminb","optim"),debug=FALSE,
       family=c("poisson","negbin","logit"),
       opt.control=list(),optim.method="BFGS",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fml}{
A formula. This formula must provide the dependent variable as well as the non linear part of the right hand side (RHS). It can be for instance \code{y~a*log(b*x + c*x^3)}. If there is no non-linear part, the RHS of the formula should be 0: e.g. \code{y~0}.
}
  \item{data}{
A data.frame containing the necessary variables to run the model. The variables of the non-linear right hand side of the formula are identified with this data.frame names. Note that no NA is allowed.
}
  \item{family}{
Character scalar. It should provide the family. Currently \code{family="poisson"}, \code{family="negbin"} and \code{family="logit"} are implemented. Note that the log link is used by default.
}
  \item{linear.fml}{
A formula with no left hand side. This formula states the linear parameters (as the constant for instance). Putting linear parameters in this formula enhances A LOT the performance of the algorithm. Example: \code{linear.fml = ~ 1} to include only the constant, or \code{linear.fml = ~ z + factor(f)} for other variables along with the constant. Note that by default there is not any linear parameter (not even the constant).
}
  \item{start}{
A list. Starting values for the non-linear parameters. ALL the parameters are to be named and given a staring value. Example: \code{start=list(a=1,b=5,c=0)}.
Though, there is an exception: if all parameters are to be given the same starting value, use \code{start.init}. Yet this is not recommended.
}
  \item{lower}{
A list. The lower bound for each of the non-linear parameters that requires one. Example: \code{lower=list(b=0,c=0)}.
Beware, if the estimated parameter is at his lower bound, problems can be raised when computing the Jacobian or the Hessian. A proper setting of \code{lower} or by using \code{d.hessian} or \code{d.jacobian} can solve these issues. See details.
}
  \item{upper}{
A list. The upper bound for each of the non-linear parameters that requires one. Example: \code{upper=list(a=10,c=50)}.
Beware, if the estimated parameter is at his upper bound, problems can be raised when computing the Jacobian or the Hessian. A proper setting of \code{upper} or a proper use of  \code{d.hessian} can solve this issue. See details.
}
  \item{env}{
An environment. You can provide an environement in which the non-linear part will be evaluated. (May be useful for some particular non-linear functions.)
}
  \item{dummy}{
Character vector. The name/s of a/some variable/s within the dataset. These variables should contain the identifier of each observation (e.g., think of it as a panel identifier).
}
  \item{start.init}{
Numeric scalar. If the argument \code{start} is not provided, or only partially filled (i.e. there remain non-linear parameters with no starting value), then the starting value of all remaining non-linear parameters is set to \code{start.init}.
}
  \item{nl.gradient}{
A formula. The user can prodide a function that computes the gradient of the non-linear part. The formula should be of the form \code{~f0(a1,x1,a2,a2)}. The important point is that it should be able to be evaluated by: \code{eval(nl.gradient[[2]],env)} where \code{env} is the working environment of the algorithm (which contains all variables and parameters). The function should return a list or a data.frame whose names are the non-linear parameters.
}
  \item{linear.start}{
Numeric named vector. The starting values of the linear part.
}
  \item{jacobian.method}{
Character scalar. Provides the method used to numerically compute the jacobian. Can be either \code{"simple"} or \code{"Richardson"}. Default is \code{"simple"}.
See the help of numDeriv for more information.
}
  \item{useHessian}{
Logical. (Only if omptimization method is \code{optim}). Should the Hessian be computed in the optimization stage? Default is \code{TRUE}.
}
  \item{d.hessian}{
Numeric scalar. It provides an argument to the function \code{hessian} of the package \code{numDeriv}. It defines the step used to compute the hessian. The default being 0.1, it can lead to problems when some parameters are at their lower or upper bound. See details for more information.
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
\details{
When the paramters are at their lower or upper bound, there can be problems when computing the Hessian. This is because the values of the paremeters are shifted to compute numerically the hessian. The defaults of those steps are 0.1 (see the help pages of \code{\link[numDeriv]{hessian}}). Thus, in the case where the non-linear part CANNOT be estimated when the parameter is beyond its bound, the hessian will not be possibly computed numerically. Thus the most straightforward way to circumvent this problem is to either rise the lower (resp. lower the upper) bound by more than 0.1, or to set \code{d.hessian} to a lower value (while slighlty rising/lowering the bound).
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

est0NL_hess = feNmlm(z~a*log(x)+b*log(y),base,start = list(a=0,b=0),
                     family="poisson", linear.fml=~1, useHessian=TRUE)

#Generating a non-linear relation
z2 = rpois(n,x + y)
base$z2 = z2

#Using a non-linear form
est1L = feNmlm(z2~0,base,~log(x)+log(y),family="poi")
est1NL = feNmlm(z2~log(a*x + b*y),base,start = list(a=1,b=2),family="poisson")
est1NL_hess = feNmlm(z2~log(a*x + b*y),base,start = list(a=1,b=2),
                     family="poisson",useHessian=TRUE)

#Using a custom Jacobian
myGrad = function(a,x,b,y){
	#Custom Jacobian
	s = a*x+b*y
	data.frame(a = x/s, b = y/s)
}

est1NL_grad = feNmlm(z2~log(a*x + b*y), base, start = list(a=1,b=2),
                     family="poisson", nl.gradient = ~myGrad(a,x,b,y))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
