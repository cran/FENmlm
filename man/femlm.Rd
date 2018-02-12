% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/femlm.R
\name{femlm}
\alias{femlm}
\title{Fixed effects maximum likelihood models}
\usage{
femlm(fml, data, family = c("poisson", "negbin", "logit", "gaussian"), NL.fml,
  cluster, useAcc = TRUE, start, lower, upper, env, start.init, offset,
  nl.gradient, linear.start = 0, jacobian.method = c("simple",
  "Richardson"), useHessian = TRUE, opt.control = list(), cores = 1,
  debug = FALSE, theta.init, ...)
}
\arguments{
\item{fml}{A formula. This formula gives the linear formula to be estimated (it is similar to a \code{lm} formula), for example: \code{fml = z~x+y}. To include cluster variables, you can 1) either insert them in this formula using a pipe (e.g. \code{fml = z~x+y|cluster1+cluster2}), or 2) either use the argment \code{cluster}. You can add a non-linear element in this formula by using the argment \code{NL.fml}. If you want to estimate only a non-linear formula without even the intercept, you can use \code{fml = z~0} in combination with \code{NL.fml}.}

\item{data}{A data.frame containing the necessary variables to run the model. The variables of the non-linear right hand side of the formula are identified with this \code{data.frame} names. Note that no \code{NA} is allowed in the variables to be used in the estimation.}

\item{family}{Character scalar. It should provide the family. The possible values are "poisson" (Poisson model with log-link, the default), "negbin" (Negative Binomial model with log-link), "logit" (LOGIT model with log-link), "gaussian" (Gaussian model).}

\item{NL.fml}{A formula. If provided, this formula represents the non-linear part of the right hand side (RHS). Note that contrary to the \code{fml} argument, the coefficients must explicitely appear in this formula. For instance, it can be \code{~a*log(b*x + c*x^3)}, where \code{a}, \code{b}, and \code{c} are the coefficients to be estimated. Note that only the RHS of the formula is to be provided, and NOT the left hand side.}

\item{cluster}{Character vector. The name/s of a/some variable/s within the dataset to be used as clusters. These variables should contain the identifier of each observation (e.g., think of it as a panel identifier).}

\item{useAcc}{Default is \code{TRUE}. Whether an acceleration algorithm (Irons and Tuck iterations) should be used to otbain the cluster coefficients when there are two clusters.}

\item{start}{A list. Starting values for the non-linear parameters. ALL the parameters are to be named and given a staring value. Example: \code{start=list(a=1,b=5,c=0)}. Though, there is an exception: if all parameters are to be given the same starting value, you can use the argument \code{start.init}.}

\item{lower}{A list. The lower bound for each of the non-linear parameters that requires one. Example: \code{lower=list(b=0,c=0)}. Beware, if the estimated parameter is at his lower bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.}

\item{upper}{A list. The upper bound for each of the non-linear parameters that requires one. Example: \code{upper=list(a=10,c=50)}. Beware, if the estimated parameter is at his upper bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.}

\item{env}{An environment. You can provide an environement in which the non-linear part will be evaluated. (May be useful for some particular non-linear functions.)}

\item{start.init}{Numeric scalar. If the argument \code{start} is not provided, or only partially filled (i.e. there remain non-linear parameters with no starting value), then the starting value of all remaining non-linear parameters is set to \code{start.init}.}

\item{offset}{A formula. An offset can be added to the estimation. It should be a formula of the form (for example) ~0.5*x**2. This offset is linearily added to the elements of the main formula 'fml'. Note that when using the argument 'NL.fml', you can directly add the offset there.}

\item{nl.gradient}{A formula. The user can prodide a function that computes the gradient of the non-linear part. The formula should be of the form \code{~f0(a1,x1,a2,a2)}. The important point is that it should be able to be evaluated by: \code{eval(nl.gradient[[2]], env)} where \code{env} is the working environment of the algorithm (which contains all variables and parameters). The function should return a list or a data.frame whose names are the non-linear parameters.}

\item{linear.start}{Numeric named vector. The starting values of the linear part. Note that you can}

\item{jacobian.method}{Character scalar. Provides the method used to numerically compute the jacobian of the non-linear part. Can be either \code{"simple"} or \code{"Richardson"}. Default is \code{"simple"}. See the help of \code{\link[numDeriv]{jacobian}} for more information.}

\item{useHessian}{Logical. Should the Hessian be computed in the optimization stage? Default is \code{TRUE}.}

\item{opt.control}{List of elements to be passed to the optimization method \code{\link[stats]{nlminb}}.}

\item{cores}{Integer, default is 1. Number of threads to be used (accelerates the algorithm via the use of openMP routines). This is particularly efficient for the negative binomial and logit models, less so for Gaussian and Poisson likelihoods (unless for large datasets).}

\item{debug}{Logical. If \code{TRUE} then the log-likelihood as well as all parameters are printed at each iteration. Default is \code{FALSE}.}

\item{theta.init}{Positive numeric scalar. The starting value of the dispersion parameter if \code{family="negbin"}. By default, the algorithm uses as a starting value the theta obtained from the model with only the intercept.}

\item{...}{Not currently used.}
}
\value{
An \code{femlm} object.
\item{coef}{The coefficients.}
\item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
\item{loglik}{The loglikelihood.}
\item{iterations}{Number of iterations of the algorithm.}
\item{n}{The number of observations.}
\item{k}{The number of parameters of the model.}
\item{call}{The call.}
\item{NL.fml}{The nonlinear formula of the call.}
\item{linear.formula}{The linear formula of the call.}
\item{ll_null}{Log-likelihood of the null model (i.e. with the intercept only).}
\item{pseudo_r2}{The adjusted pseudo R2.}
\item{naive.r2}{The R2 as if the expected predictor was the linear predictor in OLS.}
\item{message}{The convergence message from the optimization procedures.}
\item{sq.cor}{Squared correlation between the dependent variable and its expected value as given by the optimization.}
\item{hessian}{The Hessian of the parameters.}
\item{expected.predictor}{The expected predictor is the expected value of the dependent variable.}
\item{cov.unscaled}{The variance-covariance matrix of the parameters.}
\item{bounds}{Whether the coefficients were upper or lower bounded. -- This can only be the case when a non-linear formula is included and the arguments 'lower' or 'upper' are provided.}
\item{isBounded}{The logical vector that gives for each coefficient whether it was bounded or not. This can only be the case when a non-linear formula is included and the arguments 'lower' or 'upper' are provided.}
\item{se}{The standard-error of the parameters.}
\item{scores}{The matrix of the scores (first derivative for each observation).}
\item{family}{The ML family that was used for the estimation.}
\item{resids}{The difference between the dependent variable and the expected predictor.}
\item{dummies}{The sum of the cluster coefficients for each observation.}
\item{clusterNames}{The names of each cluster.}
\item{id_dummies}{The list (of length the number of clusters) of the cluser identifiers for each observation.}
\item{clusterSize}{The size of each cluster.}
\item{obsRemoved}{In the case there were clusters and some observations were removed because of only 0/1 outcome wirhin a cluster, it gives the row numbers of the observations that were removed.}
\item{clusterRemoved}{In the case there were clusters and some observations were removed because of only 0/1 outcome wirhin a cluster, it gives the list (for each cluster) of the clustr identifiers that were removed.}
\item{theta}{In the case of a negative binomial estimation: the overdispersion parameter.}
}
\description{
This function estimates maximum likelihood models (e.g., Poisson or Logit) and is efficient to handle any number of fixed effects (i.e. cluster variables). It further allows for nonlinear in parameters right hand sides.
}
\details{
This function estimates maximum likelihood models where the conditional expectations are as follows:

Gaussian likelihood:
\deqn{E(Y|X)=X\beta}{E(Y|X) = X*beta}
Poisson and Negative Binomial likelihoods:
\deqn{E(Y|X)=\exp(X\beta)}{E(Y|X) = exp(X*beta)}
where in the Negative Binomial there is the parameter \eqn{\theta}{theta} used to model the variance as \eqn{\mu+\mu^2/\theta}{mu+mu^2/theta}, with \eqn{\mu}{mu} the conditional expectation.
Logit likelihood:
\deqn{E(Y|X)=\frac{\exp(X\beta)}{1+\exp(X\beta)}}{E(Y|X) = exp(X*beta) / (1 + exp(X*beta))}

When there are one or more clusters, the conditional expectation can be written as:
\deqn{E(Y|X) = h(X\beta+\sum_{k}\sum_{m}\gamma_{m}^{k}\times C_{im}^{k}),}
where \eqn{h(.)} is the function corresponding to the likelihood function as shown before. \eqn{C^k} is the matrix associated to cluster \eqn{k} such that \eqn{C^k_{im}} is equal to 1 if observation \eqn{i} is of category \eqn{m} in cluster \eqn{k} and 0 otherwise.

When there are non linear in parameters functions, we can schematically split the set of regressors in two:
\deqn{f(X,\beta)=X^1\beta^1 + g(X^2,\beta^2)}
with first a linear term and then a non linear part expressed by the function g. That is, we add a non-linear term to the linear terms (which are \eqn{X*beta} and the cluster coefficients). It is always better (more efficient) to put into the argument \code{NL.fml} only the non-linear in parameter terms, and add all linear terms in the \code{fml} argument.
}
\examples{

#
# Linear examples
#

# Load trade data
data(trade)

# We estimate the effect of distance on trade => we account for 3 cluster effects
# 1) Poisson estimation
est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
# alternative formulation giving the same results:
# est_pois = femlm(Euros ~ log(dist_km), trade, cluster = c("Origin", "Destination", "Product"))

# 2) Log-Log Gaussian estimation
est_gaus = femlm(log(Euros+1) ~ log(dist_km)|Origin+Destination+Product, trade, family="gaussian")

# 3) Negative Binomial estimation
est_nb = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade, family="negbin")

# Comparison of the results using the function res2table
res2table(est_pois, est_gaus, est_nb)
# Now using two way clustered standard-errors
res2table(est_pois, est_gaus, est_nb, se = "twoway")

# Comparing different types of standard errors
sum_white = summary(est_pois, se = "white")
sum_oneway = summary(est_pois, se = "cluster")
sum_twoway = summary(est_pois, se = "twoway")
sum_threeway = summary(est_pois, se = "threeway")

res2table(sum_white, sum_oneway, sum_twoway, sum_threeway)


#
# Non-linear examples
#

# Generating data for a simple example
n = 100
x = rnorm(n, 1, 5)**2
y = rnorm(n, -1, 5)**2
z = rpois(n, x*y) + rpois(n, 2)
base = data.frame(x, y, z)

# Comparing the results of a 'linear' function using a 'non-linear' call
est0L = femlm(z~log(x)+log(y), base)
est0NL = femlm(z~1, base, NL.fml = ~a*log(x)+b*log(y), start = list(a=0, b=0))
# we compare the estimates with the function res2table
res2table(est0L, est0NL)

# Generating a non-linear relation
z2 = rpois(n, x + y) + rpois(n, 1)
base$z2 = z2

# Using a non-linear form
est1NL = femlm(z2~0, base, NL.fml = ~log(a*x + b*y), start = list(a=1, b=2), lower = list(a=0, b=0))
# we can't estimate this relation linearily
# => closest we can do:
est1L = femlm(z2~log(x)+log(y), base)

res2table(est1L, est1NL)

# Using a custom Jacobian for the function log(a*x + b*y)
myGrad = function(a,x,b,y){
	# Custom Jacobian
	s = a*x+b*y
	data.frame(a = x/s, b = y/s)
}

est1NL_grad = femlm(z2~0, base, NL.fml = ~log(a*x + b*y), start = list(a=1,b=2),
                     nl.gradient = ~myGrad(a,x,b,y))


}
\references{
For models with multiple fixed-effects:

Gaure, Simen, 2013, "OLS with multiple high dimensional category variables", Computational Statistics & Data Analysis 66 pp. 8--18



On the unconditionnal Negative Binomial model:

Allison, Paul D and Waterman, Richard P, 2002, "Fixed-Effects Negative Binomial Regression Models", Sociological Methodology 32(1) pp. 247--265
}
\seealso{
See also \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
}
\author{
Laurent Berge
}