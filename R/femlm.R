# Commands to genereate the help files:
# load("data/trade.RData")
# roxygen2::roxygenise(roclets = "rd")
# devtools::install(args = "--no-multiarch", build_vignettes = TRUE) # To build the vignettes


#' Fixed effects maximum likelihood models
#'
#' This function estimates maximum likelihood models (e.g., Poisson or Logit) and is efficient to handle any number of fixed effects (i.e. cluster variables). It further allows for nonlinear in parameters right hand sides.
#'
#' @param fml A formula. This formula gives the linear formula to be estimated (it is similar to a \code{lm} formula), for example: \code{fml = z~x+y}. To include cluster variables, you can 1) either insert them in this formula using a pipe (e.g. \code{fml = z~x+y|cluster1+cluster2}), or 2) either use the argment \code{cluster}. To include a non-linear in parameters element, you must use the argment \code{NL.fml}.
#' @param NL.fml A formula. If provided, this formula represents the non-linear part of the right hand side (RHS). Note that contrary to the \code{fml} argument, the coefficients must explicitely appear in this formula. For instance, it can be \code{~a*log(b*x + c*x^3)}, where \code{a}, \code{b}, and \code{c} are the coefficients to be estimated. Note that only the RHS of the formula is to be provided, and NOT the left hand side.
#' @param data A data.frame containing the necessary variables to run the model. The variables of the non-linear right hand side of the formula are identified with this \code{data.frame} names. Note that no \code{NA} is allowed in the variables to be used in the estimation. Can also be a matrix.
#' @param family Character scalar. It should provide the family. The possible values are "poisson" (Poisson model with log-link, the default), "negbin" (Negative Binomial model with log-link), "logit" (LOGIT model with log-link), "gaussian" (Gaussian model).
#' @param cluster Character vector. The name/s of a/some variable/s within the dataset to be used as clusters. These variables should contain the identifier of each observation (e.g., think of it as a panel identifier).
#' @param na.rm Logical, default is \code{FALSE}. If the variables necessary for the estimation contain NAs and \code{na.rm = TRUE}, then all observations containing NA are removed prior to estimation and a warning message is raised detailing the number of observations removed.
#' @param useAcc Default is \code{TRUE}. Whether an acceleration algorithm (Irons and Tuck iterations) should be used to otbain the cluster coefficients when there are two or more clusters.
#' @param NL.start (For NL models only) A list of starting values for the non-linear parameters. ALL the parameters are to be named and given a staring value. Example: \code{NL.start=list(a=1,b=5,c=0)}. Though, there is an exception: if all parameters are to be given the same starting value, you can use the argument \code{NL.start.init}.
#' @param lower (For NL models only) A list. The lower bound for each of the non-linear parameters that requires one. Example: \code{lower=list(b=0,c=0)}. Beware, if the estimated parameter is at his lower bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.
#' @param upper (For NL models only) A list. The upper bound for each of the non-linear parameters that requires one. Example: \code{upper=list(a=10,c=50)}. Beware, if the estimated parameter is at his upper bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.
#' @param env (For NL models only) An environment. You can provide an environement in which the non-linear part will be evaluated. (May be useful for some particular non-linear functions.)
#' @param NL.start.init (For NL models only) Numeric scalar. If the argument \code{NL.start} is not provided, or only partially filled (i.e. there remain non-linear parameters with no starting value), then the starting value of all remaining non-linear parameters is set to \code{NL.start.init}.
#' @param offset A formula. An offset can be added to the estimation. It should be a formula of the form (for example) ~0.5*x**2. This offset is linearily added to the elements of the main formula 'fml'. Note that when using the argument 'NL.fml', you can directly add the offset there.
#' @param nl.gradient (For NL models only) A formula. The user can prodide a function that computes the gradient of the non-linear part. The formula should be of the form \code{~f0(a1,x1,a2,a2)}. The important point is that it should be able to be evaluated by: \code{eval(nl.gradient[[2]], env)} where \code{env} is the working environment of the algorithm (which contains all variables and parameters). The function should return a list or a data.frame whose names are the non-linear parameters.
#' @param linear.start Numeric named vector. The starting values of the linear part.
#' @param jacobian.method Character scalar. Provides the method used to numerically compute the jacobian of the non-linear part. Can be either \code{"simple"} or \code{"Richardson"}. Default is \code{"simple"}. See the help of \code{\link[numDeriv]{jacobian}} for more information.
#' @param useHessian Logical. Should the Hessian be computed in the optimization stage? Default is \code{TRUE}.
#' @param opt.control List of elements to be passed to the optimization method \code{\link[stats]{nlminb}}. See the help page of \code{\link[stats]{nlminb}} for more information.
#' @param cores Integer, default is 1. Number of threads to be used (accelerates the algorithm via the use of openMP routines). This is particularly efficient for the negative binomial and logit models, less so for the Gaussian and Poisson likelihoods (unless for very large datasets).
#' @param verbose Integer, default is 0. It represents the level of information that should be reported during the optimisation process. If \code{verbose=0}: nothing is reported. If \code{verbose=1}: the value of the coefficients and the likelihood are reported. If \code{verbose=2}: \code{1} + information on the computing tiime of the null model, the cluster coefficients and the hessian are reported.
#' @param theta.init Positive numeric scalar. The starting value of the dispersion parameter if \code{family="negbin"}. By default, the algorithm uses as a starting value the theta obtained from the model with only the intercept.
#' @param precision.cluster Precision used to obtain the fixed-effects (ie cluster coefficients). Defaults to \code{1e-5}. It corresponds to the maximum absolute difference allowed between two iterations. Argument \code{precision.cluster} cannot be lower than \code{10000*.Machine$double.eps}.
#' @param itermax.cluster Maximum number of iterations in the step obtaining the fixed-effects (only in use for 2+ clusters). Default is 10000.
#' @param itermax.deriv Maximum number of iterations in the step obtaining the derivative of the fixed-effects (only in use for 2+ clusters). Default is 5000.
#' @param ... Not currently used.
#'
#' @details
#' This function estimates maximum likelihood models where the conditional expectations are as follows:
#'
#' Gaussian likelihood:
#' \deqn{E(Y|X)=X\beta}{E(Y|X) = X*beta}
#' Poisson and Negative Binomial likelihoods:
#' \deqn{E(Y|X)=\exp(X\beta)}{E(Y|X) = exp(X*beta)}
#' where in the Negative Binomial there is the parameter \eqn{\theta}{theta} used to model the variance as \eqn{\mu+\mu^2/\theta}{mu+mu^2/theta}, with \eqn{\mu}{mu} the conditional expectation.
#' Logit likelihood:
#' \deqn{E(Y|X)=\frac{\exp(X\beta)}{1+\exp(X\beta)}}{E(Y|X) = exp(X*beta) / (1 + exp(X*beta))}
#'
#' When there are one or more clusters, the conditional expectation can be written as:
#' \deqn{E(Y|X) = h(X\beta+\sum_{k}\sum_{m}\gamma_{m}^{k}\times C_{im}^{k}),}
#' where \eqn{h(.)} is the function corresponding to the likelihood function as shown before. \eqn{C^k} is the matrix associated to cluster \eqn{k} such that \eqn{C^k_{im}} is equal to 1 if observation \eqn{i} is of category \eqn{m} in cluster \eqn{k} and 0 otherwise.
#'
#' When there are non linear in parameters functions, we can schematically split the set of regressors in two:
#' \deqn{f(X,\beta)=X^1\beta^1 + g(X^2,\beta^2)}
#' with first a linear term and then a non linear part expressed by the function g. That is, we add a non-linear term to the linear terms (which are \eqn{X*beta} and the cluster coefficients). It is always better (more efficient) to put into the argument \code{NL.fml} only the non-linear in parameter terms, and add all linear terms in the \code{fml} argument.
#'
#' To estimate only a non-linear formula without even the intercept, you must exclude the intercept from the linear formula by using, e.g., \code{fml = z~0}.
#'
#' The over-dispersion parameter of the Negative Binomial family, theta, is capped at 10,000. If theta reaches this high value, it means that there is no overdispersion.
#'
#' @return
#' An \code{femlm} object.
#' \item{coefficients}{The named vector of coefficients.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
#' \item{loglik}{The loglikelihood.}
#' \item{iterations}{Number of iterations of the algorithm.}
#' \item{n}{The number of observations.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{call}{The call.}
#' \item{fml}{The linear formula of the call.}
#' \item{ll_null}{Log-likelihood of the null model (i.e. with the intercept only).}
#' \item{pseudo_r2}{The adjusted pseudo R2.}
#' \item{message}{The convergence message from the optimization procedures.}
#' \item{sq.cor}{Squared correlation between the dependent variable and the expected predictor (i.e. fitted.values) obtained by the estimation.}
#' \item{hessian}{The Hessian of the parameters.}
#' \item{fitted.values}{The fitted values are the expected value of the dependent variable for the fitted model: that is \eqn{E(Y|X)}.}
#' \item{cov.unscaled}{The variance-covariance matrix of the parameters.}
#' \item{se}{The standard-error of the parameters.}
#' \item{scores}{The matrix of the scores (first derivative for each observation).}
#' \item{family}{The ML family that was used for the estimation.}
#' \item{residuals}{The difference between the dependent variable and the expected predictor.}
#' \item{sumFE}{The sum of the fixed-effects for each observation.}
#' \item{offset}{The offset formula.}
#' \item{NL.fml}{The nonlinear formula of the call.}
#' \item{bounds}{Whether the coefficients were upper or lower bounded. -- This can only be the case when a non-linear formula is included and the arguments 'lower' or 'upper' are provided.}
#' \item{isBounded}{The logical vector that gives for each coefficient whether it was bounded or not. This can only be the case when a non-linear formula is included and the arguments 'lower' or 'upper' are provided.}
#' \item{clusterNames}{The names of each cluster.}
#' \item{id_dummies}{The list (of length the number of clusters) of the cluster identifiers for each observation.}
#' \item{clusterSize}{The size of each cluster.}
#' \item{obsRemoved}{In the case there were clusters and some observations were removed because of only 0/1 outcome within a cluster, it gives the row numbers of the observations that were removed.}
#' \item{clusterRemoved}{In the case there were clusters and some observations were removed because of only 0/1 outcome within a cluster, it gives the list (for each cluster) of the clustr identifiers that were removed.}
#' \item{theta}{In the case of a negative binomial estimation: the overdispersion parameter.}
#'
#' @seealso
#' See also \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @references
#'
#' Berge, Laurent, 2018, "Efficient estimation of maximum likelihood models with multiple fixed-effects: the R package FENmlm." CREA Discussion Papers, 13 (\url{https://wwwen.uni.lu/content/download/110162/1299525/file/2018_13}).
#'
#' For models with multiple fixed-effects:
#'
#' Gaure, Simen, 2013, "OLS with multiple high dimensional category variables", Computational Statistics & Data Analysis 66 pp. 8--18
#'
#' On the unconditionnal Negative Binomial model:
#'
#' Allison, Paul D and Waterman, Richard P, 2002, "Fixed-Effects Negative Binomial Regression Models", Sociological Methodology 32(1) pp. 247--265
#'
#' @examples
#'
#' #
#' # Linear examples
#' #
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade => we account for 3 cluster effects
#' # 1) Poisson estimation
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#' # alternative formulation giving the same results:
#' # est_pois = femlm(Euros ~ log(dist_km), trade, cluster = c("Origin", "Destination", "Product"))
#'
#' # 2) Log-Log Gaussian estimation (with same clusters)
#' est_gaus = update(est_pois, log(Euros+1) ~ ., family="gaussian")
#'
#' # 3) Negative Binomial estimation
#' est_nb = update(est_pois, family="negbin")
#'
#' # Comparison of the results using the function res2table
#' res2table(est_pois, est_gaus, est_nb)
#' # Now using two way clustered standard-errors
#' res2table(est_pois, est_gaus, est_nb, se = "twoway")
#'
#' # Comparing different types of standard errors
#' sum_white = summary(est_pois, se = "white")
#' sum_oneway = summary(est_pois, se = "cluster")
#' sum_twoway = summary(est_pois, se = "twoway")
#' sum_threeway = summary(est_pois, se = "threeway")
#'
#' res2table(sum_white, sum_oneway, sum_twoway, sum_threeway)
#'
#'
#' #
#' # Example of Equivalences
#' #
#'
#' # equivalence with glm poisson
#' est_glm <- glm(Euros ~ log(dist_km) + factor(Origin) +
#'             factor(Destination) + factor(Product), trade, family = poisson)
#'
#' # coefficient estimates + Standard-error
#' summary(est_glm)$coefficients["log(dist_km)", ]
#' est_pois$coeftable
#'
#' # equivalence with lm
#' est_lm <- lm(log(Euros+1) ~ log(dist_km) + factor(Origin) +
#'             factor(Destination) + factor(Product), trade)
#'
#' # coefficient estimates + Standard-error
#' summary(est_lm)$coefficients["log(dist_km)", ]
#' summary(est_gaus, dof_correction = TRUE)$coeftable
#'
#' #
#' # Non-linear examples
#' #
#'
#' # Generating data for a simple example
#' n = 100
#' x = rnorm(n, 1, 5)**2
#' y = rnorm(n, -1, 5)**2
#' z1 = rpois(n, x*y) + rpois(n, 2)
#' base = data.frame(x, y, z1)
#'
#' # Estimating a 'linear' relation:
#' est1_L = femlm(z1 ~ log(x) + log(y), base)
#' # Estimating the same 'linear' relation using a 'non-linear' call
#' est1_NL = femlm(z1 ~ 1, base, NL.fml = ~a*log(x)+b*log(y), NL.start = list(a=0, b=0))
#' # we compare the estimates with the function res2table (they are identical)
#' res2table(est1_L, est1_NL)
#'
#' # Now generating a non-linear relation (E(z2) = x + y + 1):
#' z2 = rpois(n, x + y) + rpois(n, 1)
#' base$z2 = z2
#'
#' # Estimation using this non-linear form
#' est2_NL = femlm(z2~0, base, NL.fml = ~log(a*x + b*y),
#'                NL.start = list(a=1, b=2), lower = list(a=0, b=0))
#' # we can't estimate this relation linearily
#' # => closest we can do:
#' est2_L = femlm(z2~log(x)+log(y), base)
#'
#' # Difference between the two models:
#' res2table(est2_L, est2_NL)
#'
#' # Plotting the fits:
#' plot(x, z2, pch = 18)
#' points(x, fitted(est2_L), col = 2, pch = 1)
#' points(x, fitted(est2_NL), col = 4, pch = 2)
#'
#'
#' # Using a custom Jacobian for the function log(a*x + b*y)
#' myGrad = function(a,x,b,y){
#' 	s = a*x+b*y
#' 	data.frame(a = x/s, b = y/s)
#' }
#'
#' est2_NL_grad = femlm(z2~0, base, NL.fml = ~log(a*x + b*y),
#'                      NL.start = list(a=1,b=2), nl.gradient = ~myGrad(a,x,b,y))
#'
#'
femlm <- function(fml, data, family=c("poisson", "negbin", "logit", "gaussian"), NL.fml, cluster, na.rm = FALSE, useAcc=TRUE, NL.start, lower, upper, env, NL.start.init, offset, nl.gradient, linear.start=0, jacobian.method=c("simple", "Richardson"), useHessian=TRUE, opt.control=list(), cores = 1, verbose=0, theta.init, precision.cluster, itermax.cluster = 10000, itermax.deriv = 5000, ...){

	# use of the conjugate gradient in the gaussian case to get
	# the cluster coefficients
	accDeriv = TRUE

	jacobian.method <- match.arg(jacobian.method)
	family = match.arg(family)

	# Some settings (too complicated to be tweaked by the user)
	# Nber of evaluations of the NL part to be kept in memory
	# Default keeps the two last evaluations
	NLsave = 2

	# DOTS
	dots = list(...)

	# Future parameter in development: na.rm
	# na.rm = ifelse(is.null(dots$na.rm), FALSE, dots$na.rm)
	isNA_sample = FALSE

	ptm = proc.time()

	# DEPRECATED INFORMATION
	# I initially called the cluster dummies... I keep it for compatibility
	if(missing(cluster) && "dummy" %in% names(dots)) cluster = dots$dummy
	if("linear.fml" %in% names(dots)) stop("Argument 'linear.fml' is deprecated, now use 'fml' in combination with 'NL.fml'.")
	if(missing(NL.start) && "start" %in% names(dots)) {
		warning("Argument 'start' is deprecated.\nUse 'NL.start' instead.", immediate. = TRUE)
		NL.start = dots$start
	}
	if(missing(NL.start.init) && "start.init" %in% names(dots)) {
		warning("Argument 'start.init' is deprecated.\nUse 'NL.start.init' instead.", immediate. = TRUE)
		NL.start.init = dots$start.init
	}

	#
	# The clusters => controls + setup
	if(class(fml) != "formula") stop("The argument 'fml' must be a formula.")
	if(length(fml) != 3) stop("The formula must be two sided.\nEG: a~exp(b/x), or a~0 if there is no linear part.")

	FML = Formula::Formula(fml)
	n_rhs = length(FML)[2]

	if(n_rhs > 2){
		stop("The argument 'fml' cannot contain more than two parts separated by a pipe ('|').")
	}

	if(n_rhs == 2){
		if(missing(cluster) || length(cluster) == 0){
			cluster = formula(FML, lhs = 0, rhs = 2)
			fml = formula(FML, lhs = 1, rhs = 1)
		} else {
			stop("To add cluster variables: either include them in argument 'fml' using a pipe ('|'), either use the argument 'cluster'. You cannot use both!")
		}
	}

	# Other paramaters
	if(!is.null(dots$debug) && dots$debug) verbose = 2
	d.hessian = dots$d.hessian

	#
	# cores argument
	FENmlm_CORES = 1
	if(!length(cores) == 1 || !is.numeric(cores) || !(cores%%1) == 0 || cores < 0){
		stop("The argument 'cores' must be an integer greater than 0 and lower than the number of threads of the computer.")
	} else if(cores == 1){
		isMulticore = FALSE
	} else if(is.na(parallel::detectCores())){
		# This can happen...
		isMulticore = FALSE
		warning("The number of cores has been set to 1 because the function detectCores() could not evaluate the maximum number of nodes.")
	} else if(cores > parallel::detectCores()){
		stop("The argument 'cores' must be lower or equal to the number of possible threads (equal to ", parallel::detectCores(), ") in this computer.")
	} else {
		FENmlm_CORES = cores
		isMulticore = TRUE
	}

	famFuns = switch(family,
						  poisson = ml_poisson(),
						  negbin = ml_negbin(),
						  logit = ml_logit(),
						  gaussian = ml_gaussian())

	call = match.call(expand.dots = FALSE)
	# cleaning the call from update method
	# we drop 'hidden' arguments for a clean call
	call$"..." = NULL

	if(is.matrix(data)){
		if(is.null(colnames(data))){
			stop("If argument data is to be a matrix, its columns must be named.")
		}
		data = as.data.frame(data)
	}
	# The conversion of the data (due to data.table)
	if(!"data.frame" %in% class(data)){
		stop("The argument 'data' must be a data.frame or a matrix.")
	}
	if("data.table" %in% class(data)){
		# this is a local change only
		class(data) = "data.frame"
	}

	dataNames = names(data)

	# The LHS must contain only values in the DF
	namesLHS = all.vars(fml[[2]])
	if(!all(namesLHS %in% dataNames)) stop("Some elements on the LHS of the formula are not in the dataset:\n", paste0(namesLHS[!namesLHS %in% dataNames], collapse = ", "))

	# Now the nonlinear part:

	isNA_NL = FALSE # initialisation of NAs flag (FALSE is neutral)
	if(!missing(NL.fml) && !is.null(NL.fml)){

		if(!class(NL.fml) == "formula") stop("Argument 'NL.fml' must be a formula.")

		isNonLinear = TRUE
		nl.call = NL.fml[[length(NL.fml)]]

		allnames = all.vars(nl.call)
		nonlinear.params = allnames[!allnames %in% dataNames]
		nonlinear.varnames = allnames[allnames %in% dataNames]

		if(length(nonlinear.params) == 0){
			warning("As there is no parameter to estimate in argument 'NL.fml', this argument is ignored.\nIf you want to add an offset, use argument 'offset'.")
		}

		# Control for NAs
		if(anyNA(data[, nonlinear.varnames])){
			if(!na.rm){
				# Default
				varWithNA = nonlinear.varnames[which(apply(data[, nonlinear.varnames, FALSE], 2, anyNA))]
				text = show_vars_limited_width(varWithNA)
				stop("Some variables in 'NL.fml' contain NA. NAs are not supported, please remove them first (or use na.rm). FYI, the variables are:\n", text, call. = FALSE)
			} else {
				# If Na.rm => we keep track of the NAs
				isNA_NL = is.na(rowSums(data[, nonlinear.varnames]))
				isNA_sample = isNA_sample || TRUE
			}

		}

	} else {
		isNonLinear = FALSE
		nl.call = 0
		allnames = nonlinear.params = nonlinear.varnames = character(0)
	}

	# The dependent variable: lhs==left_hand_side
	lhs = as.vector(eval(fml[[2]], data))

	# creation de l'environnement
	if(missing(env)) env <- new.env()
	else stopifnot(class(env)=="environment")

	#
	# First check
	#

	isNA_y = FALSE # initialisation of NAs flag (FALSE is neutral)

	if(anyNA(lhs)){
		if(!na.rm){
			# Default behavior
			stop("The left hand side of the fomula has NA values. Please provide data without NA (or use na.rm).")
		} else {
			# If Na.rm => we keep track of the NAs
			isNA_y = is.na(lhs)
			isNA_sample = isNA_sample || TRUE
		}

		# We repeat twice the controls, but one for a cleaned y
		# It's a bit "ugly" but I don't recreate unecessary vectors this way
		# We check that the dep var is not a constant

		lhs_clean = lhs[!isNA_y]

		# we check the var is not a constant
		if(var(lhs_clean) == 0){
			stop("The dependent variable is a constant. Estimation cannot be done.")
		}
		if(family %in% c("poisson", "negbin") & any(lhs_clean<0)) stop("Negative values of the dependant variable \nare not allowed for the \"", family, "\" family.", sep="")
		if(family %in% c("logit") & !all(lhs_clean==0 | lhs_clean==1)) stop("The dependant variable has values different from 0 or 1.\nThis is not allowed with the \"logit\" family.")

	} else if(any(!is.finite(lhs))){
		stop("The dependent variable contains non-finite values.")
	} else {
		# regular controls when there is no NA

		# We check that the dep var is not a constant
		if(var(lhs) == 0){
			stop("The dependent variable is a constant. Estimation cannot be done.")
		}
		if(family %in% c("poisson", "negbin") & any(lhs<0)) stop("Negative values of the dependant variable \nare not allowed for the \"", family, "\" family.", sep="")
		if(family %in% c("logit") & !all(lhs==0 | lhs==1)) stop("The dependant variable has values different from 0 or 1.\nThis is not allowed with the \"logit\" family.")
	}

	#
	# Controls and setting of the linear part:
	#

	isLinear = FALSE
	linear.varnames = all.vars(fml[[3]])

	if(length(linear.varnames) > 0 || attr(terms(fml), "intercept") == 1){
		isLinear = TRUE
		linear.fml = fml
	}

	isNA_L = FALSE # initialisation of NAs flag (FALSE is neutral)
	if(isLinear){

		if(!all(linear.varnames %in% dataNames)) stop(paste("In 'fml', some variables are not in the data:\n", paste(linear.varnames[!linear.varnames%in%dataNames], collapse=', '), ".", sep=""))

		if(!missing(cluster) && length(cluster) != 0){
			# if dummies are provided, we make sure there is an
			# intercept so that factors can be handled properly
			linear.fml = update(linear.fml, ~.+1)
		}

		#
		# We construct the linear matrix
		#

		# we look at whether there are factor-like variables to be evaluated
		# if there is factors => model.matrix
		types = sapply(data[, dataNames %in% linear.varnames, FALSE], class)
		if(grepl("factor", deparse(linear.fml)) || any(types %in% c("character", "factor"))){
			useModel.matrix = TRUE
		} else {
			useModel.matrix = FALSE
		}

		if(useModel.matrix){
			linear.mat = stats::model.matrix(linear.fml, data)
		} else {
			linear.mat = prepare_matrix(linear.fml, data)
		}

		linear.params <- colnames(linear.mat)
		# N_linear <- length(linear.params)
		if(anyNA(linear.mat)){

			if(!na.rm){
				# Default behavior: no NA tolerance
				quiNA = apply(linear.mat, 2, anyNA)
				whoIsNA = linear.params[quiNA]
				text = show_vars_limited_width(whoIsNA)

				stop("Evaluation of the linear part returns NA. NAs are not supported, please remove them before running this function (or use na.rm). FYI the variables with NAs are:\n", text)
			} else {
				# If Na.rm => we keep track of the NAs
				isNA_L = is.na(rowSums(linear.mat))
				isNA_sample = isNA_sample || TRUE
			}

		}

		if(!is.numeric(linear.start)) stop("'linear.start' must be numeric!")

	} 	else {
		linear.params <- linear.start <- linear.varnames <- NULL
		useModel.matrix = FALSE
	}

	params <- c(nonlinear.params, linear.params)
	lparams <- length(params)
	varnames <- c(nonlinear.varnames, linear.varnames)

	# Attention les parametres non lineaires peuvent etre vides
	if(length(nonlinear.params)==0) isNL = FALSE
	else isNL = TRUE



	#
	# Handling Clusters ####
	#

	isDummy = FALSE
	if(!is.null(dots$clusterFromUpdate) && dots$clusterFromUpdate){
		# Cluster information coming from the update method

		# means that there is no modification of past clusters
		object = dots$object

		# we retrieve past information
		dum_all = object$id_dummies
		obs2remove = object$obsRemoved
		dummyOmises = object$clusterRemoved
		cluster = object$clusterNames
		nbCluster = object$clusterSize
		Q = length(nbCluster)

		# If obsRemoved => need to modify the data base
		if(length(obs2remove) > 0){
			data = data[-obs2remove, ]

			# We recreate the linear matrix
			if(isLinear) {
				if(useModel.matrix){
					# means there are factors
					linear.mat = stats::model.matrix(linear.fml, data)
				} else {
					linear.mat = linear.mat[-obs2remove, ]
				}
			}

			lhs = eval(fml[[2]], data)
		}

		# We still need to recreate some objects though
		isDummy = TRUE
		names(dum_all) = NULL
		dumMat_cpp = matrix(unlist(dum_all), ncol = Q) - 1

		dum_names = sum_y_all = obs_per_cluster_all = list()
		for(i in 1:Q){
			k = nbCluster[i]
			dum = dum_all[[i]]
			sum_y_all[[i]] = cpp_tapply_vsum(k, lhs, dum)
			obs_per_cluster_all[[i]] = cpp_table(k, dum)
			dum_names[[i]] = attr(dum_all[[i]], "clust_names")
		}

	} else if(!missing(cluster) && length(cluster)!=0){
		# The main cluster construction

		isDummy = TRUE

		isClusterFml = FALSE
		if(is.character(cluster) && any(!cluster %in% names(data))){
			var_problem = cluster[!cluster%in%names(data)]
			stop("The argument 'cluster' must be variable names! Cluster(s) not in the data: ", paste0(var_problem, collapse = ", "), ".")
		} else if(!is.character(cluster)){
			# means the cluster is a formula
			cluster_fml = cluster

			# we check that the cluster variables are indeed in the data
			cluster_vars = all.vars(cluster)
			if(!all(cluster_vars %in% names(data))){
				var_problem = cluster_vars[!cluster_vars %in% names(data)]
				stop("The following 'cluster' variable", ifelse(length(var_problem) == 1, " is", "s are"), " not in the data: ", paste0(var_problem, collapse = ", "), ".")
			}

			cluster_mat = model.frame(cluster_fml, data, na.action = NULL)
			# we change cluster to become a vector of characters
			cluster = names(cluster_mat)
			isClusterFml = TRUE
		} else {
			# the clusters are valid cluster names
			cluster_mat = data[, cluster, drop = FALSE]
		}

		# We change factors to character
		isFactor = sapply(cluster_mat, is.factor)
		if(any(isFactor)){
			for(i in which(isFactor)){
				cluster_mat[[i]] = as.character(cluster_mat[[i]])
			}
		}

		isNA_cluster = FALSE
		if(anyNA(cluster_mat)){
			if(!na.rm){
				# Default behavior, NA not allowed
				var_problem = cluster[sapply(cluster_mat, anyNA)]
				stop("The cluster variables contain NAs, this is not allowed (or use na.rm).\nFYI, the clusters with NA are: ", paste0(var_problem, collapse = ", "), ".")
			} else {
				# If Na.rm => we keep track of the NAs
				isNA_cluster = apply(cluster_mat, MARGIN = 1, anyNA)
				isNA_sample = isNA_sample || TRUE
			}
		}

		if(isNA_sample){
			# we remove NAs from the clusters only
			# after, we'll remove them from the data too
			isNA_full = isNA_y | isNA_L | isNA_NL | isNA_cluster
			nbNA = sum(isNA_full)
			nobs = nrow(data)
			if(nbNA == nobs){
				stop("All observations contain NAs. Estimation cannot be done.")
			}

			message_NA = paste0(nbNA, " observations removed because of NA values. (Breakup: LHS:", sum(isNA_y), ", RHS:", sum(isNA_L + isNA_NL), ", Cluster: ", sum(isNA_cluster), ")")

			# warning(nbNA, " observations removed because of NA values. (Breakup: LHS:", sum(isNA_y), ", RHS:", sum(isNA_L + isNA_NL), ", Cluster: ", sum(isNA_cluster), ")", immediate. = TRUE)

			# we drop the NAs from the cluster matrix
			cluster_mat = cluster_mat[!isNA_full, , drop = FALSE]
			obs2remove_NA = which(isNA_full)
			index_noNA = (1:nobs)[!isNA_full]

			# we change the LHS variable
			lhs = eval(fml[[2]], data[-obs2remove_NA, ])
		}

		Q = length(cluster)
		dum_all = dum_names = list()
		sum_y_all = obs_per_cluster_all = list()
		obs2remove = c()
		dummyOmises = list()
		for(i in 1:Q){

			dum_raw = cluster_mat[[cluster[i]]]

			# in order to avoid "unclassed" values > real nber of classes: we re-factor the cluster

			dum_names[[i]] = thisNames = getItems(dum_raw)
			dum = quickUnclassFactor(dum_raw)

			dum_all[[i]] = dum
			k = length(thisNames)

			# We delete "all zero" outcome
			sum_y_all[[i]] = sum_y_clust = cpp_tapply_vsum(k, lhs, dum)
			obs_per_cluster_all[[i]] = n_perClust = cpp_table(k, dum)

			if(family %in% c("poisson", "negbin")){
				qui = which(sum_y_clust==0)
			} else if(family == "logit"){
				qui = which(sum_y_clust==0 | sum_y_clust==n_perClust)
			} else if(family == "gaussian"){
				qui = NULL
			}

			if(length(qui>0)){
				# We first delete the data:
				dummyOmises[[i]] = thisNames[qui]
				obs2remove = unique(c(obs2remove, which(dum %in% qui)))
			} else {
				dummyOmises[[i]] = character(0)
			}
		}

		# We remove the problems
		if(length(obs2remove)>0){

			# update of the cluster matrix
			cluster_mat = cluster_mat[-obs2remove, , drop = FALSE]

			# update of the lhs
			lhs = lhs[-obs2remove]

			# Then we recreate the dummies
			for(i in 1:Q){

				dum_raw = cluster_mat[[cluster[i]]]

				dum_names[[i]] = getItems(dum_raw)
				dum = quickUnclassFactor(dum_raw)

				dum_all[[i]] = dum
				k = length(dum_names[[i]])

				# We also recreate these values
				sum_y_all[[i]] = cpp_tapply_vsum(k, lhs, dum)
				obs_per_cluster_all[[i]] = cpp_table(k, dum)

			}

			# Then the warning message
			nb_missing = sapply(dummyOmises, length)
			message_cluster = paste0(paste0(nb_missing, collapse = "/"), " cluster", ifelse(sum(nb_missing) == 1, "", "s"), " (", length(obs2remove), " observations) removed because of only ", ifelse(family=="logit", "zero/one", "zero"), " outcomes.")
			if(isNA_sample){
				warning(message_NA, "\n  ", message_cluster)
			} else {
				warning(message_cluster)
			}

			names(dummyOmises) = cluster
		}

		if(isNA_sample){
			# we update the value of obs2remove (will contain both NA and removed bc of outcomes)
			if(length(obs2remove) > 0){
				obs2remove_cluster = index_noNA[obs2remove]
			} else {
				obs2remove_cluster = c()
			}

			obs2remove = sort(c(obs2remove_NA, obs2remove_cluster))
		}

		# We compute two values that will be useful to compute the derivative wrt clusters
		dumMat_cpp = matrix(unlist(dum_all), ncol = Q) - 1
		nbCluster = sapply(dum_all, max)

	} else {
		# There is no cluster
		Q = 0

		# NA management is needed to create obs2remove
		if(isNA_sample){
			# after, we'll remove them from the data too
			isNA_full = isNA_y | isNA_L | isNA_NL
			nbNA = sum(isNA_full)
			nobs = nrow(data)
			if(nbNA == nobs){
				stop("All observations contain NAs. Estimation cannot be done.")
			}

			warning(nbNA, " observations removed because of NA values. (Breakup: LHS:", sum(isNA_y), ", RHS:", sum(isNA_L + isNA_NL), ")")

			# we drop the NAs from the cluster matrix
			obs2remove = which(isNA_full)
		} else {
			obs2remove = c()
		}
	}

	# NA & problem management
	if(length(obs2remove) > 0){
		# we kick out the problems (both NA related and cluster related)
		data = data[-obs2remove, ]

		# We recreate the linear matrix and the LHS
		if(isLinear) {
			if(useModel.matrix){
				# means there are factors
				linear.mat = stats::model.matrix(linear.fml, data)
			} else {
				linear.mat = linear.mat[-obs2remove, ]
			}
		}

		lhs = eval(fml[[2]], data)
	}

	# If presence of clusters => we exclude the intercept
	if(Q > 0){
		# If there is a linear intercept, we withdraw it
		# We drop the intercept:
		if(isLinear && "(Intercept)" == colnames(linear.mat)){
			var2remove = which(colnames(linear.mat) == "(Intercept)")
			if(ncol(linear.mat) == length(var2remove)){
				isLinear = FALSE
				linear.params = NULL
				params <- nonlinear.params
				lparams <- length(params)
				varnames <- nonlinear.varnames
			} else{
				linear.mat = linear.mat[, -var2remove, drop=FALSE]
				linear.params <- colnames(linear.mat)
				# N_linear <- length(linear.params)
				params <- c(nonlinear.params, linear.params)
				lparams <- length(params)
				varnames <- c(nonlinear.varnames, linear.varnames)
			}
		}
	}


	#
	# Checks for MONKEY TEST
	#

	if(lparams==0 & Q==0) stop("No parameter to be estimated.")
	if(!is.logical(useHessian)) stop("'useHessian' must be of type 'logical'!")

	#
	# Controls: The non linear part
	#

	if(isNL){
		if(missing(NL.start.init)){
			if(missing(NL.start)) stop("There must be starting values for NL parameters. Please use argument NL.start (or NL.start.init).")
			if(typeof(NL.start)!="list") stop("NL.start must be a list.")
			if(any(!nonlinear.params %in% names(NL.start))) stop(paste("Some NL parameters have no starting values:\n", paste(nonlinear.params[!nonlinear.params %in% names(NL.start)], collapse=", "), ".", sep=""))

			# we restrict NL.start to the nonlinear.params
			NL.start = NL.start[nonlinear.params]
		} else {
			if(length(NL.start.init)>1) stop("NL.start.init must be a scalar.")
			if(!is.numeric(NL.start.init)) stop("NL.start.init must be numeric!")
			if(!is.finite(NL.start.init)) stop("Infinites values as starting values, you must be kidding me...")

			if(missing(NL.start)){
				NL.start <- list()
				NL.start[nonlinear.params] <- NL.start.init
			} else {
				if(typeof(NL.start)!="list") stop("NL.start must be a list.")
				if(any(!names(NL.start) %in% params)) stop(paste("Some parameters in 'NL.start' are not in the formula:\n", paste(names(NL.start)[!names(NL.start) %in% params], collapse=", "), ".", sep=""))

				missing.params <- nonlinear.params[!nonlinear.params%in%names(NL.start)]
				NL.start[missing.params] <- NL.start.init
			}
		}
	} else {
		NL.start <- list()
	}

	#
	# The upper and lower limits
	#

	if(!missing(lower) && !is.null(lower)){

		if(typeof(lower)!="list"){
			stop("'lower' MUST be a list.")
		}

		lower[params[!params %in% names(lower)]] <- -Inf
		lower <- unlist(lower[params])

	}	else {
		lower <- rep(-Inf, lparams)
		names(lower) <- params
	}

	if(!missing(upper) && !is.null(upper)){

		if(typeof(upper)!="list"){
			stop("'upper' MUST be a list.")
		}

		upper[params[!params %in% names(upper)]] <- Inf
		upper <- unlist(upper[params])

	}	else {
		upper <- rep(Inf, lparams)
		names(upper) <- params
	}

	lower <- c(lower)
	upper <- c(upper)

	#
	# Controls: user defined gradient
	#

	if(!missing(nl.gradient)){
		isGradient = TRUE
		if(class(nl.gradient)!="formula" | length(nl.gradient)==3) stop("'nl.gradient' must be a formula like, for ex., ~f0(a1, x1, a2, x2). f0 giving the gradient.")
	} else {
		isGradient = FALSE
	}

	if(!is.null(d.hessian)){
		hessianArgs = list(d=d.hessian)
	} else hessianArgs = NULL
	assign("hessianArgs", hessianArgs, env)

	#
	# Offset
	#

	offset.value = 0
	if(!missing(offset) && !is.null(offset)){

		# control
		if(!class(offset) == "formula"){
			stop("Argument 'offset' must be a formula (e.g. ~ 1+x^2).")
		}

		if(length(offset) != 2){
			stop("Argument 'offset' must be a formula of the type (e.g.): ~ 1+x^2.")
		}

		offset.call = offset[[length(offset)]]
		vars.offset = all.vars(offset.call)

		if(any(!vars.offset %in% dataNames)){
			var_missing = vars.offset[!vars.offset %in% dataNames]
			stop("Some variable in the argument 'offset' are not in the data:\n", paste0(var_missing, sep=", "), ".")
		}

		offset.value = eval(offset.call, data)

		if(anyNA(offset.value)){
			stop("Evaluating the argument 'offset' lead to NA values.")
		}

	}

	assign("offset.value", offset.value, env)

	#
	# PRECISION
	#

	n = length(lhs)
	# The main precision
	if (!missing(precision.cluster) && !is.null(precision.cluster)){
		if(!length(precision.cluster)==1 || !is.numeric(precision.cluster) || precision.cluster <= 0 || precision.cluster >1){
			stop("If provided, argument 'precision.cluster' must be a strictly positive scalar lower than 1.")
		} else if(precision.cluster < 10000*.Machine$double.eps){
			stop("Argument 'precision.cluster' cannot be lower than ", signif(10000*.Machine$double.eps))
		}
		eps.cluster = precision.cluster
	} else {
		eps.cluster = 1e-5 # min(10**-(log10(n) + Q/3), 1e-5)
	}

	if(!is.numeric(itermax.cluster) || length(itermax.cluster) > 1 || itermax.cluster < 1){
		stop("Argument itermax.cluster must be an integer greater than 0.")
	}
	if(!is.numeric(itermax.deriv) || length(itermax.deriv) > 1 || itermax.deriv < 1){
		stop("Argument itermax.deriv must be an integer greater than 0.")
	}

	# other precisions
	eps.NR = ifelse(is.null(dots$eps.NR), eps.cluster/100, dots$eps.NR)
	eps.deriv = ifelse(is.null(dots$eps.deriv), 1e-4, dots$eps.deriv)


	# Initial checks are done
	nonlinear.params <- names(NL.start) #=> in the order the user wants

	# starting values of all parameters (initialized with the NL ones):
	start = NL.start

	# control for the linear start => we can provide coefficients
	# from past estimations. Coefficients that are not provided are set
	# to 0
	if(length(linear.start)>1){
		what = linear.start[linear.params]
		what[is.na(what)] = 0
		linear.start = what
	}
	start[linear.params] <- linear.start
	params <- names(start)
	start <- unlist(start)
	start <- c(start)
	lparams <- length(params)
	names(start) <- params

	# The right order of upper and lower
	upper = upper[params]
	lower = lower[params]

	#
	# The MODEL0 => to get the init of the theta for the negbin
	#

	# Bad location => rethink the design of the code
	assign(".famFuns", famFuns, env)
	assign(".family", family, env)
	assign("nobs", length(lhs), env)
	assign(".lhs", lhs, env)
	assign(".isMulticore", isMulticore, env)
	assign(".CORES", FENmlm_CORES, env)
	assign(".verbose", verbose, env)

	if(missing(theta.init)){
		theta.init = NULL
	} else {
		if(!is.null(theta.init) && (!is.numeric(theta.init) || length(theta.init)!=1 || theta.init<=0)){
			stop("the argument 'theta.init' must be a strictly positive scalar.")
		}
	}

	model0 <- get_model_null(env, theta.init)
	theta.init = model0$theta

	# For the negative binomial:
	if(family == "negbin"){
		params = c(params, ".theta")
		start = c(start, theta.init)
		names(start) = params
		upper = c(upper, 10000)
		lower = c(lower, 1e-3)
	}

	# On balance les donnees a utiliser dans un nouvel environnement
	for(i in varnames) assign(i, data[[i]], env)
	if(isLinear) assign("linear.mat", linear.mat, env)
	if(isGradient) assign(".call_gradient", nl.gradient[[2]], env)

	####
	#### Sending to the env ####
	####

	useExp_clusterCoef = family %in% c("poisson")

	# The dummies
	assign("isDummy", isDummy, env)
	if(isDummy){
		assign(".dum_all", dum_all, env)
		assign(".dumMat_cpp", dumMat_cpp, env)
		assign(".nbCluster", nbCluster, env)
		assign(".sum_y", sum_y_all, env)
		assign(".tableCluster", obs_per_cluster_all, env)

		# the saved dummies
		if(!is.null(dots$clusterStart)){
			# information on starting values coming from update method
			doExp = ifelse(useExp_clusterCoef, exp, I)

			if(dots$clusterFromUpdate || length(obs2remove) == 0){
				# Means it's the full cluster properly given
				assign(".savedDummy", doExp(dots$clusterStart), env)
			} else {
				assign(".savedDummy", doExp(dots$clusterStart[-obs2remove]), env)
			}

		} else if(useExp_clusterCoef){
			assign(".savedDummy", rep(1, length(lhs)), env)
		} else {
			assign(".savedDummy", rep(0, length(lhs)), env)
		}

		assign(".orderCluster", NULL, env)
		if(isMulticore || family %in% c("negbin", "logit")){
			# we also add this variable used in cpp
			orderCluster_all = list()
			for(i in 1:Q){
				orderCluster_all[[i]] = order(dum_all[[i]]) - 1
			}
			# orderCluster_mat = do.call("cbind", orderCluster_all)
			assign(".orderCluster", orderCluster_all, env)
		}
	}
	# other
	assign(".nl.call", nl.call, env)
	assign(".lhs", lhs, env)
	assign("isNL", isNL, env)
	assign("isLinear", isLinear, env)
	assign("isGradient", isGradient, env)
	assign("linear.params", linear.params, env)
	assign("nonlinear.params", nonlinear.params, env)
	assign("params", params, env)
	assign("nobs", length(lhs), env)
	assign(".verbose", verbose, env)
	assign("jacobian.method", jacobian.method, env)
	assign(".famFuns", famFuns, env)
	assign(".family", family, env)
	assign(".iter", 0, env)
	# Pour gerer les valeurs de mu:
	assign(".coefMu", list(), env)
	assign(".valueMu", list(), env)
	assign(".valueExpMu", list(), env)
	assign(".wasUsed", TRUE, env)
	# Pour les valeurs de la Jacobienne non lineaire
	assign(".JC_nbSave", 0, env)
	assign(".JC_nbMaxSave", 1, env)
	assign(".JC_savedCoef", list(), env)
	assign(".JC_savedValue", list(), env)
	# PRECISION
	assign(".eps.cluster", eps.cluster, env)
	assign(".eps.NR", eps.NR, env)
	assign(".eps.deriv", eps.deriv, env)
	# ITERATIONS
	assign(".itermax.cluster", itermax.cluster, env)
	assign(".itermax.deriv", itermax.deriv, env)
	# OTHER
	assign(".useAcc", useAcc, env)
	assign(".warn_0_Hessian", FALSE, env)
	assign(".warn_overfit_logit", FALSE, env)

	# To monitor how the clusters are computed (if the problem is difficult or not)
	assign(".firstIterCluster", 1e10, env) # the number of iterations in the first run
	assign(".firstRunCluster", TRUE, env) # flag for first enrty in get_dummies
	assign(".iterCluster", 1e10, env) # the previous number of cluster iterations
	assign(".evolutionLL", Inf, env) # diff b/w two successive LL
	assign(".pastLL", 0, env)
	assign(".iterLastPrecisionIncrease", 0, env) # the last iteration when precision was increased
	assign(".nbLowIncrease", 0, env) # number of successive evaluations with very low LL increments
	assign(".nbIterOne", 0, env) # nber of successive evaluations with only 1 iter to get the clusters
	assign(".difficultConvergence", FALSE, env)

	# Same for derivatives
	assign(".derivDifficultConvergence", FALSE, env)
	assign(".firstRunDeriv", TRUE, env) # flag for first entry in derivative
	assign(".accDeriv", TRUE, env) # Logical: flag for accelerating deriv
	assign(".iterDeriv", 1e10, env) # number of iterations in the derivatives step


	#
	# if there is only the intercept and cluster => we estimate only the clusters
	#

	if(!isLinear && !isNonLinear && Q>0){
		if(family == "negbin"){
			stop("To estimate the negative binomial model, you need at least one variable. (The estimation of the model with only the clusters is not implemented.)")
		}
		results = femlm_only_clusters(env, model0, cluster, dum_names)
		results$call = match.call()
		results$fml = fml
		if(length(obs2remove)>0){
			results$obsRemoved = obs2remove
			results$clusterRemoved = dummyOmises
		}
		results$onlyCluster = TRUE
		return(results)
	}


	# On teste les valeurs initiales pour informer l'utilisateur
	for(var in nonlinear.params) assign(var, start[var], env)

	if(isNL){
		mu = NULL
		try(mu <- eval(nl.call, envir= env), silent = FALSE)
		if(is.null(mu)){
			# the non linear part could not be evaluated - ad hoc message
			stop("The non-linear part (NL.fml) could not be evaluated. There may be a problem in 'NL.fml'.")
		}

		# DEPREC: the NL part should return stg the same length
		# # the special case of the constant
		# if(length(mu) == 1){
		# 	mu = rep(mu, nrow(data))
		# }

		# No numeric vectors
		if(!is.vector(mu) || !is.numeric(mu)){
			stop("Evaluation of NL.fml should return a numeric vector. (This is currently not the case)")
		}

		# Handling NL.fml errors
		if(length(mu) != nrow(data)){
			stop("Evaluation of NL.fml leads to ", length(mu), " observations while there are ", nrow(data), " observations in the data base. They should be the same lenght.")
		}

		if(anyNA(mu)){
			stop("Evaluating NL.fml leads to NA values (which are forbidden). Maybe it's a problem with the starting values, maybe it's another problem.")
		}

	} else {
		mu = eval(nl.call, envir = env)
	}


	# On sauvegarde les valeurs de la partie non lineaire
	assign(".nbMaxSave", NLsave, env) # nombre maximal de valeurs a sauvegarder
	assign(".nbSave", 1, env)  # nombre de valeurs totales sauvegardees
	assign(".savedCoef", list(start[nonlinear.params]), env)
	assign(".savedValue", list(mu), env)
	if(isLinear) {
		mu <- mu + c(linear.mat%*%unlist(start[linear.params]))
	}

	if(anyNA(mu)){
		stop("Evaluating the left hand side leads to NA values.")
	}

	# Check of the user-defined gradient, if given
	if(isGradient){
		for(nom in nonlinear.params) assign(nom, start[nom], env)
		test <- eval(nl.gradient[[2]], envir=env)
		if(!class(test)%in%c("list", "data.frame")) stop("The function called by 'nl.gradient' must return an object of type 'list' or 'data.frame'.")
		if(!all(nonlinear.params%in%names(test))) stop(paste("The gradient must return a value for each parameter. Some are missing:\n", paste(nonlinear.params[!nonlinear.params%in%names(test)], collapse=", "), ".", sep=""))
		if(!all(names(test)%in%nonlinear.params)) warning(paste("Some values given by 'nl.gradient' are not in the parameters:\n", paste(names(test)[!names(test)%in%nonlinear.params], collapse=", "), ".", sep=""))
		if(mean(sapply(test[nonlinear.params], length))!=length(lhs)) stop("Strange, the length of the vector returned by 'nl.gradient' does not match with the data.")
		#we save 1 gradient:
		jacob.mat = as.matrix(test[nonlinear.params])
		assign(".JC_nbSave", 1, env)
		assign(".JC_savedCoef", list(start[nonlinear.params]), env)
		assign(".JC_savedValue", list(jacob.mat), env)
	}

	# Mise en place du calcul du gradient
	gradient = femlm_gradient
	hessian <- NULL
	if(useHessian) hessian <- femlm_hessian

	# GIVE PARAMS
	if(!is.null(dots$give.params) && dots$give.params) return(list(coef=start, env=env))

	if(verbose >= 2) cat("Setup in ", (proc.time() - ptm)[3], "s\n", sep="")

	#
	# Maximizing the likelihood
	#

	opt <- NULL
	opt <- stats::nlminb(start=start, objective=femlm_ll, env=env, lower=lower, upper=upper, gradient=gradient, hessian=hessian, control=opt.control)

	if(is.null(opt)){
		stop("Could not achieve maximization.")
	}

	convStatus = TRUE
	warningMessage = ""
	if(!opt$message %in% c("X-convergence (3)", "relative convergence (4)", "both X-convergence and relative convergence (5)")){
		# warning("[femlm] The optimization algorithm did not converge, the results are not reliable. Use function diagnostic() to see what's wrong.", call. = FALSE)
		warningMessage = " The optimization algorithm did not converge, the results are not reliable."
		convStatus = FALSE
	}

	####
	#### After Maximization ####
	####

	coef <- opt$par

	# The Hessian
	hessian = femlm_hessian(coef, env=env)
	# we add the names of the non linear variables in the hessian
	if(isNonLinear || family == "negbin"){
		dimnames(hessian) = list(params, params)
	}

	# we create the Hessian without the bounded parameters
	hessian_noBounded = hessian

	# Handling the bounds
	if(!isNonLinear){
		NL.fml = NULL
		bounds = NULL
		isBounded = NULL
	} else {
		# we report the bounds & if the estimated parameters are bounded
		upper_bound = upper[nonlinear.params]
		lower_bound = lower[nonlinear.params]

		# 1: are the estimated parameters at their bounds?
		coef_NL = coef[nonlinear.params]
		isBounded = rep(FALSE, length(params))
		isBounded[1:length(coef_NL)] = (coef_NL == lower_bound) | (coef_NL == upper_bound)

		# 2: we save the bounds
		upper_bound_small = upper_bound[is.finite(upper_bound)]
		lower_bound_small = lower_bound[is.finite(lower_bound)]
		bounds = list()
		if(length(upper_bound_small) > 0) bounds$upper = upper_bound_small
		if(length(lower_bound_small) > 0) bounds$lower = lower_bound_small
		if(length(bounds) == 0){
			bounds = NULL
		}

		# 3: we update the Hessian (basically, we drop the bounded element)
		if(any(isBounded)){
			hessian_noBounded = hessian[-which(isBounded), -which(isBounded), drop = FALSE]

			boundText = ifelse(coef_NL == upper_bound, "Upper bounded", "Lower bounded")[isBounded]

			attr(isBounded, "type") = boundText
		}

	}

	# Variance

	var <- NULL
	try(var <- solve(hessian_noBounded), silent = TRUE)
	if(is.null(var)){
		warningMessage = paste(warningMessage, "The information matrix is singular (likely presence of collinearity).")

		var = hessian_noBounded*NA
		se = diag(var)
	} else {
		se = diag(var)
		se[se < 0] = NA
		se = sqrt(se)
	}

	# Warning message
	if(nchar(warningMessage) > 0){
		warning("[femlm]:", warningMessage, " Use function diagnostic() to see what's wrong.")
	}

	# To handle the bounded coefficient, we set its SE to NA
	if(any(isBounded)){
		se = se[params]
		names(se) = params
	}

	zvalue <- coef/se
	pvalue <- 2*pnorm(-abs(zvalue))

	# We add the information on the bound for the se & update the var to drop the bounded vars
	se_format = se
	if(any(isBounded)){
		se_format[!isBounded] = decimalFormat(se_format[!isBounded])
		se_format[isBounded] = boundText
	}

	coeftable <- data.frame("Estimate"=coef, "Std. Error"=se_format, "z value"=zvalue, "Pr(>|z|)"=pvalue, stringsAsFactors = FALSE)
	names(coeftable) <- c("Estimate", "Std. Error", "z value",  "Pr(>|z|)")
	row.names(coeftable) <- params

	attr(se, "type") = attr(coeftable, "type") = "Standard"

	mu_both = get_mu(coef, env, final = TRUE)
	mu = mu_both$mu
	exp_mu = mu_both$exp_mu

	# calcul pseudo r2
	loglik <- -opt$objective # moins car la fonction minimise
	ll_null <- model0$loglik
	# degres de liberte
	df_k = length(coef)
	if(isDummy) df_k = df_k + sum(sapply(dum_all, max) - 1) + 1
	# dummies are constrained, they don't have full dof (cause you need to take one value off for unicity)
	pseudo_r2 <- 1 - (loglik-df_k)/(ll_null-1)

	# Calcul residus
	expected.predictor = famFuns$expected.predictor(mu, exp_mu, env)
	residuals = lhs - expected.predictor

	# calcul squared corr
	if(sd(expected.predictor) == 0){
		sq.cor = NA
	} else {
		sq.cor = stats::cor(lhs, expected.predictor)**2
	}

	# The scores
	scores = femlm_scores(coef, env)
	if(isNonLinear){
		# we add the names of the non linear params in the score
		colnames(scores) = params
	}


	res <- list(coefficients=coef, coeftable=coeftable, loglik=loglik, iterations=opt$iterations, n=length(lhs), nparams=df_k, call=call, fml=fml, ll_null=ll_null, pseudo_r2=pseudo_r2, message=opt$message, convStatus=convStatus, sq.cor=sq.cor, fitted.values=expected.predictor, hessian=hessian, cov.unscaled=var, se=se, scores=scores, family=family, residuals=residuals)

	# Other optional elements
	if(!missing(offset)){
		res$offset = offset
	}

	if(!is.null(NL.fml)){
		res$NL.fml = NL.fml
		if(!is.null(bounds)){
			res$bounds = bounds
			res$isBounded = isBounded
		}
	}

	# Dummies
	if(isDummy){
		dummies = attr(mu, "mu_dummies")
		if(useExp_clusterCoef){
			dummies = rpar_log(dummies, env)
		}

		res$sumFE = dummies
		res$clusterNames = cluster

		id_dummies = list()
		for(i in 1:length(cluster)){
			dum = dum_all[[i]]
			attr(dum, "clust_names") = as.character(dum_names[[i]])
			id_dummies[[cluster[i]]] = dum
		}

		res$id_dummies = id_dummies
		clustSize = sapply(dum_all, max)
		names(clustSize) = cluster
		res$clusterSize = clustSize
		if(length(obs2remove)>0){
			res$obsRemoved = obs2remove
			res$clusterRemoved = dummyOmises
		}
	}

	if(family == "negbin"){
		theta = coef[".theta"]
		res$theta = theta

		if(theta > 1000){
			warning("Very high value of theta (", theta, "). There is no sign of overdisperion, you may consider a Poisson model.")
		}

	}

	class(res) <- "femlm"

	if(verbose > 0){
		cat("\n")
	}

	return(res)
}

femlm_only_clusters <- function(env, model0, cluster, dum_names){
	# Estimation with only the cluster coefficients

	#
	# 1st step => computing the dummies
	#

	nobs = get("nobs", env)
	family = get(".family", env)
	offset.value = get("offset.value", env)

	if(family == "negbin"){
		coef[[".theta"]] = model0$theta
	} else {
		coef = list()
	}

	# indicator of whether we compute the exp(mu)
	useExp = family %in% c("poisson", "logit", "negbin")
	useExp_clusterCoef = family %in% c("poisson")

	# mu, using the offset
	mu_noDum = offset.value
	if(length(mu_noDum) == 1) mu_noDum = rep(mu_noDum, nobs)

	# we create the exp of mu => used for later functions
	exp_mu_noDum = NULL
	if(useExp_clusterCoef){
		exp_mu_noDum = rpar_exp(mu_noDum, env)
	}

	dummies = getDummies(mu_noDum, exp_mu_noDum, env, coef)

	exp_mu = NULL
	if(useExp_clusterCoef){
		# despite being called mu, it is in fact exp(mu)!!!
		exp_mu = exp_mu_noDum*dummies
		mu = rpar_log(exp_mu, env)
	} else {
		mu = mu_noDum + dummies
		if(useExp){
			exp_mu = rpar_exp(mu, env)
		}
	}

	#
	# 2nd step => saving information
	#

	dum_all = get(".dum_all", env)
	famFuns = get(".famFuns", env)
	lhs = get(".lhs", env)

	# The log likelihoods
	loglik = famFuns$ll(lhs, mu, exp_mu, env, coef)
	ll_null = model0$loglik

	# degres de liberte
	df_k = sum(sapply(dum_all, max) - 1) + 1
	pseudo_r2 = 1 - loglik/ll_null # NON Adjusted

	# Calcul residus
	expected.predictor = famFuns$expected.predictor(mu, exp_mu, env)
	residuals = lhs - expected.predictor

	# calcul squared corr
	if(sd(expected.predictor) == 0){
		sq.cor = NA
	} else {
		sq.cor = stats::cor(lhs, expected.predictor)**2
	}

	# calcul r2 naif
	naive.r2 = 1 - sum(residuals**2) / sum((lhs - mean(lhs))**2)

	res = list(loglik=loglik, n=length(lhs), nparams=df_k, call=call, ll_null=ll_null, pseudo_r2=pseudo_r2, naive.r2=naive.r2, sq.cor=sq.cor, expected.predictor=expected.predictor, residuals=residuals, family=family)
	#
	# Information on the dummies

	if(useExp_clusterCoef){
		dummies = rpar_log(dummies, env)
	}

	res$sumFE = dummies
	res$clusterNames = cluster

	id_dummies = list()
	for(i in 1:length(cluster)){
		dum = dum_all[[i]]
		attr(dum, "clust_names") = as.character(dum_names[[i]])
		id_dummies[[cluster[i]]] = dum
	}

	res$id_dummies = id_dummies
	clustSize = sapply(dum_all, max)
	names(clustSize) = cluster
	res$clusterSize = clustSize

	if(family == "negbin"){
		theta = coef[".theta"]
		res$theta = theta
	}

	class(res) = "femlm"

	return(res)
}

femlm_hessian <- function(coef, env){
	# Computes the hessian

	verbose = get(".verbose", env)
	if(verbose >= 2) ptm = proc.time()
	params <- get("params", env)
	names(coef) <- params
	nonlinear.params <- get("nonlinear.params", env)
	k <- length(nonlinear.params)
	isNL <- get("isNL", env)
	hessianArgs = get("hessianArgs", env)
	famFuns = get(".famFuns", env)
	family = get(".family", env)
	y = get(".lhs", env)
	isDummy = get("isDummy", env)
	mu_both = get_savedMu(coef, env)
	mu = mu_both$mu
	exp_mu = mu_both$exp_mu

	jacob.mat = get_Jacobian(coef, env)

	ll_d2 = famFuns$ll_d2(y, mu, exp_mu, coef)
	if(isDummy){
		dxi_dbeta = deriv_xi(jacob.mat, ll_d2, env, coef)
		jacob.mat = jacob.mat + dxi_dbeta
	} else dxi_dbeta = 0

	hessVar = crossprod(jacob.mat, jacob.mat * ll_d2)

	if(isNL){
		# we get the 2nd derivatives
		z = numDeriv::genD(evalNLpart, coef[nonlinear.params], env=env, method.args = hessianArgs)$D[, -(1:k), drop=FALSE]
		ll_dl = famFuns$ll_dl(y, mu, exp_mu, coef=coef, env=env)
		id_r = rep(1:k, 1:k)
		id_c = c(sapply(1:k, function(x) 1:x), recursive=TRUE)
		H = matrix(0, nrow=k, ncol=k)
		H[cbind(id_r, id_c)] = H[cbind(id_r, id_c)] = colSums(z*ll_dl)
	} else H = 0

	# on ajoute la partie manquante
	if(isNL) hessVar[1:k, 1:k] = hessVar[1:k, 1:k] + H

	if(family=="negbin"){
		theta = coef[".theta"]
		ll_dx_dother = famFuns$ll_dx_dother(y, mu, exp_mu, coef, env)

		if(isDummy){
			dxi_dother = deriv_xi_other(ll_dx_dother, ll_d2, env, coef)
		} else {
			dxi_dother = 0
		}

		# calcul des derivees secondes vav de theta
		h.theta.L = famFuns$hess.thetaL(theta, jacob.mat, y, dxi_dbeta, dxi_dother, ll_d2, ll_dx_dother)
		hessVar = cbind(hessVar, h.theta.L)
		h.theta = famFuns$hess_theta_part(theta, y, mu, exp_mu, dxi_dother, ll_dx_dother, ll_d2, env)
		hessVar = rbind(hessVar, c(h.theta.L, h.theta))

	} else if(family=="tobit"){
		sigma = coef[".sigma"]
		h.sigma.L = famFuns$hess.sigmaL(sigma, jacob.mat, y, mu, dxi_dbeta, ll_d2)
		hessVar = cbind(hessVar, h.sigma.L)
		h.sigma = famFuns$hess.sigma(sigma, y, mu)
		hessVar = rbind(hessVar, c(h.sigma.L, h.sigma))
	}

	if(anyNA(hessVar)){
		stop("NaN in the Hessian, can be due to a possible overfitting problem.\nIf so, to have an idea of what's going on, you can reduce the value of the argument 'rel.tol' of the nlminb algorithm using the argument 'opt.control = list(rel.tol=?)' with ? the new value.")
	}

	# warn_0_Hessian = get(".warn_0_Hessian", env)
	# if(!warn_0_Hessian && any(diag(hessVar) == 0)){
	# 	# We apply the warning only once
	# 	var_problem = params[diag(hessVar) == 0]
	# 	warning("Some elements of the diagonal of the hessian are equal to 0: likely presence of collinearity. FYI the problematic variables are: ", paste0(var_problem, collapse = ", "), ".", immediate. = TRUE)
	# 	assign(".warn_0_Hessian", TRUE, env)
	# }

	# if(verbose >= 2) cat("Hessian: ", (proc.time()-ptm)[3], "s\n", sep="")
	- hessVar
}

femlm_gradient <- function(coef, env){
	# cat("gradient:\n") ; print(as.vector(coef))

	params = get("params", env)
	names(coef) = params
	nonlinear.params = get("nonlinear.params", env)
	linear.params = get("linear.params", env)
	famFuns = get(".famFuns", env)
	family = get(".family", env)
	y = get(".lhs", env)
	mu_both = get_savedMu(coef, env)
	mu = mu_both$mu
	exp_mu = mu_both$exp_mu

	# calcul de la jacobienne
	res <- list() #stocks the results

	# cat("\tgetting jacobian")
	# ptm = proc.time()
	jacob.mat = get_Jacobian(coef, env)
	# cat("in", (proc.time()-ptm)[3], "s.\n")

	# cat("\tComputing gradient ")
	# ptm = proc.time()
	# res = famFuns$grad(jacob.mat, y, mu, env, coef)
	res = getGradient(jacob.mat, y, mu, exp_mu, env, coef)
	# cat("in", (proc.time()-ptm)[3], "s.\n")
	names(res) = c(nonlinear.params, linear.params)

	if(family=="negbin"){
		theta = coef[".theta"]
		res[".theta"] = famFuns$grad.theta(theta, y, mu, exp_mu, env)
	}

	return(-unlist(res[params]))
}

femlm_scores <- function(coef, env){
	# Computes the scores (Jacobian)
	params = get("params", env)
	names(coef) <- params
	famFuns = get(".famFuns", env)
	family = get(".family", env)
	y = get(".lhs", env)
	mu_both = get_savedMu(coef, env)
	mu = mu_both$mu
	exp_mu = mu_both$exp_mu

	jacob.mat = get_Jacobian(coef, env)
	scores = getScores(jacob.mat, y, mu, exp_mu, env, coef)

	if(family=="negbin"){
		theta = coef[".theta"]
		score.theta = famFuns$scores.theta(theta, y, mu, exp_mu, env)
		scores = cbind(scores, score.theta)
	}

	return(scores)
}

femlm_ll <- function(coef, env){
	# Log likelihood

	# misc funs
	iter = get(".iter", env) + 1
	assign(".iter", iter, env)
	pastLL = get(".pastLL", env)
	verbose = get(".verbose", env)
	ptm = proc.time()
	if(verbose >= 1){
		coef_names = sapply(names(coef), charShorten, width = 10)
		coef_line = paste0(coef_names, ": ", signif(coef), collapse = " -- ")
		cat("\nIter", iter, "- Coefficients:", coef_line, "\n")
	}

	# computing the LL
	famFuns = get(".famFuns", env)
	family = get(".family", env)
	y <- get(".lhs", env)

	if(any(is.na(coef))) stop("Divergence... (some coefs are NA)\nTry option verbose=2 to figure out the problem.")

	mu_both = get_mu(coef, env)
	mu = mu_both$mu
	exp_mu = mu_both$exp_mu

	# for the NEGBIN, we add the coef
	ll = famFuns$ll(y, mu, exp_mu, env, coef)

	evolutionLL = ll - pastLL
	assign(".evolutionLL", evolutionLL, env)
	assign(".pastLL", ll, env)

	if(iter == 1) evolutionLL = "--"

	if(verbose >= 1) cat("LL = ", ll, " (", (proc.time()-ptm)[3], "s)\tevol = ", evolutionLL, "\n", sep = "")
	if(ll==(-Inf)) return(1e308)
	return(-ll) # je retourne -ll car la fonction d'optimisation minimise
}

evalNLpart = function(coef, env){
	# cat("Enter evalNLpart : ", as.vector(coef), "\n")
	# fonction qui evalue la partie NL
	isNL = get("isNL", env)
	if(!isNL) return(0)

	nonlinear.params <- get("nonlinear.params", env)
	nl.call <- get(".nl.call", env)
	nbSave = get(".nbSave", env)
	nbMaxSave = get(".nbMaxSave", env)
	savedCoef = get(".savedCoef", env)
	savedValue = get(".savedValue", env)
	if(!is.null(names(coef))){
		coef = coef[nonlinear.params]
	} else if (length(coef) != length(nonlinear.params)){
		stop("Problem with the length of the NL coefficients.")
	}

	if(nbMaxSave == 0){
		for(var in nonlinear.params) assign(var, coef[var], env)
		y_nl <- eval(nl.call, envir= env)

		# we check problems
		if(anyNA(y_nl)){
			stop("Evaluation of non-linear part returns NAs. The coefficients were: ", paste0(nonlinear.params, " = ", signif(coef[nonlinear.params], 3)), ".")
		}

		return(y_nl)
	}

	for(i in nbSave:1){
		#les valeurs les + recentes sont en derniere position
		if(all(coef == savedCoef[[i]])){
			return(savedValue[[i]])
		}
	}

	#Si la valeur n'existe pas, on la sauvegarde
	#on met les valeurs les plus recentes en derniere position
	for(var in nonlinear.params) assign(var, coef[var], env)
	y_nl = eval(nl.call, envir = env)

	# we check problems
	if(anyNA(y_nl)){
		stop("Evaluation of non-linear part returns NAs. The coefficients were: ", paste0(nonlinear.params, " = ", signif(coef[nonlinear.params], 3)), ".")
	}

	if(nbSave<nbMaxSave){
		savedCoef[[nbSave+1]] = coef
		savedValue[[nbSave+1]] = y_nl
		assign(".nbSave", nbSave+1, env)
	} else if(nbMaxSave>1){
		tmp = list()
		tmp[[nbSave]] = coef
		tmp[1:(nbSave-1)] = savedCoef[2:nbSave]
		savedCoef = tmp

		tmp = list()
		tmp[[nbSave]] = y_nl
		tmp[1:(nbSave-1)] = savedValue[2:nbSave]
		savedValue = tmp
	} else{
		savedCoef = list(coef)
		savedValue = list(y_nl)
	}

	# cat("computed NL part:", as.vector(coef), "\n")

	assign(".savedCoef", savedCoef, env)
	assign(".savedValue", savedValue, env)
	return(y_nl)
}

get_mu = function(coef, env, final = FALSE){
	# This function computes the RHS of the equation
	# mu_L => to save one matrix multiplication
	isNL = get("isNL", env)
	isLinear = get("isLinear", env)
	isDummy = get("isDummy", env)
	nobs = get("nobs", env)
	params = get("params", env)
	family = get(".family", env)
	offset.value = get("offset.value", env)
	names(coef) = params

	# UseExp: indicator if the family needs to use exp(mu) in the likelihoods:
	#     this is useful because we have to compute it only once (save computing time)
	# useExp_clusterCoef: indicator if we use the exponential of mu to obtain the cluster coefficients
	#     if it is TRUE, it will mean that the dummy will be equal
	#     to exp(mu_dummies) despite being named mu_dummies
	useExp = family %in% c("poisson", "logit", "negbin")
	useExp_clusterCoef = family %in% c("poisson")

	# For managing mu:
	coefMu = get(".coefMu", env)
	valueMu = get(".valueMu", env)
	valueExpMu = get(".valueExpMu", env)
	wasUsed = get(".wasUsed", env)
	if(wasUsed){
		coefMu = valueMu = valueExpMu = list()
		assign(".wasUsed", FALSE, env)
	}

	if(length(coefMu)>0){
		for(i in 1:length(coefMu)){
			if(all(coef==coefMu[[i]])){
				return(list(mu = valueMu[[i]], exp_mu = valueExpMu[[i]]))
			}
		}
	}

	if(isNL){
		muNL = evalNLpart(coef, env)
	} else muNL = 0

	if(isLinear){
		linear.params = get("linear.params", env)
		linear.mat = get("linear.mat", env)
		mu_L = c(linear.mat %*% coef[linear.params])
	} else mu_L = 0

	mu_noDum = muNL + mu_L + offset.value

	# Detection of overfitting issues with the logit model:
	if(family == "logit"){
		warn_overfit_logit = get(".warn_overfit_logit", env)

		if(!warn_overfit_logit && max(abs(mu_noDum)) >= 300){
			# overfitting => now finding the precise cause
			# browser()
			if(!isNL || (isLinear && max(abs(mu_L)) >= 100)){
				# we create the matrix with the coefficients to find out the guy
				mat_L_coef = linear.mat * matrix(coef[linear.params], nrow(linear.mat), 2, byrow = TRUE)
				max_var = apply(abs(mat_L_coef), 2, max)
				best_suspect = linear.params[which.max(max_var)]
				warning("in femlm(): Likely presence of an overfitting problem. One suspect variable is: ", best_suspect, ".", immediate. = TRUE, call. = FALSE)
			} else {
				warning("in femlm(): Likely presence of an overfitting problem due to the non-linear part.", immediate. = TRUE, call. = FALSE)
			}

			assign(".warn_overfit_logit", TRUE, env)
		}
	}


	# we create the exp of mu => used for later functions
	exp_mu_noDum = NULL
	if(useExp_clusterCoef){
		exp_mu_noDum = rpar_exp(mu_noDum, env)
	}

	if(isDummy){
		# we get back the last dummy
		mu_dummies = getDummies(mu_noDum, exp_mu_noDum, env, coef, final)
	} else {
		if(useExp_clusterCoef){
			mu_dummies = 1
		} else {
			mu_dummies = 0
		}
	}

	# We add the value of the dummy to mu and we compute the exp if necessary
	exp_mu = NULL
	if(useExp_clusterCoef){
		# despite being called mu_dummies, it is in fact exp(mu_dummies)!!!
		exp_mu = exp_mu_noDum*mu_dummies
		mu = rpar_log(exp_mu, env)
	} else {
		mu = mu_noDum + mu_dummies
		if(useExp){
			exp_mu = rpar_exp(mu, env)
		}
	}

	if(isDummy){
		# BEWARE, if useExp_clusterCoef, it is equal to exp(mu_dummies)
		attr(mu, "mu_dummies") = mu_dummies
	}

	if(length(mu)==0) mu = rep(mu, nobs)

	# we save the value of mu:
	coefMu = append(coefMu, list(coef))
	valueMu = append(valueMu, list(mu))
	valueExpMu = append(valueExpMu, list(exp_mu))
	assign(".coefMu", coefMu, env)
	assign(".valueMu", valueMu, env)
	assign(".valueExpMu", valueExpMu, env)

	return(list(mu = mu, exp_mu = exp_mu))
}

get_savedMu = function(coef, env){
	# This function gets the mu without computation
	# It follows a LL evaluation
	coefMu = get(".coefMu", env)
	valueMu = get(".valueMu", env)
	valueExpMu = get(".valueExpMu", env)
	assign(".wasUsed", TRUE, env)

	if(length(coefMu)>0) for(i in 1:length(coefMu)) if(all(coef==coefMu[[i]])){
		# cat("coef nb:", i, "\n")
		return(list(mu = valueMu[[i]], exp_mu = valueExpMu[[i]]))
	}

	stop("Problem in \"get_savedMu\":\n gradient did not follow LL evaluation.")
}

get_Jacobian = function(coef, env){
	# retrieves the Jacobian of the "rhs"
	params <- get("params", env)
	names(coef) <- params
	isNL <- get("isNL", env)
	isLinear <- get("isLinear", env)
	isGradient = get("isGradient", env)

	if(isNL){
		nonlinear.params = get("nonlinear.params", env)
		jacob.mat = get_NL_Jacobian(coef[nonlinear.params], env)
	} else jacob.mat = c()

	if(isLinear){
		linear.mat = get("linear.mat", env)
		if(is.null(dim(jacob.mat))){
			jacob.mat = linear.mat
		} else {
			jacob.mat = cbind(jacob.mat, linear.mat)
		}
	}

	return(jacob.mat)
}

get_NL_Jacobian = function(coef, env){
	# retrieves the Jacobian of the non linear part
	#cat("In NL JAC:\n")
	#print(coef)
	nbSave = get(".JC_nbSave", env)
	nbMaxSave = get(".JC_nbMaxSave", env)
	savedCoef = get(".JC_savedCoef", env)
	savedValue = get(".JC_savedValue", env)

	nonlinear.params <- get("nonlinear.params", env)
	coef = coef[nonlinear.params]

	if(nbSave>0) for(i in nbSave:1){
		#les valeurs les + recentes sont en derniere position
		if(all(coef == savedCoef[[i]])){
			# cat("Saved value:", as.vector(coef), "\n")
			return(savedValue[[i]])
		}
	}

	#Si la valeur n'existe pas, on la sauvegarde
	#on met les valeurs les plus recentes en derniere position
	isGradient <- get("isGradient", env)
	if(isGradient){
		call_gradient <- get(".call_gradient", env)
		#we send the coef in the environment
		for(var in nonlinear.params) assign(var, coef[var], env)
		jacob.mat <- eval(call_gradient, envir=env)
		jacob.mat <- as.matrix(as.data.frame(jacob.mat[nonlinear.params]))
	} else {
		jacobian.method <- get("jacobian.method", env)
		jacob.mat <- numDeriv::jacobian(evalNLpart, coef, env=env, method=jacobian.method)
	}

	#Controls:
	if(anyNA(jacob.mat)){
		qui <- which(apply(jacob.mat, 2, function(x) anyNA(x)))
		variables <- nonlinear.params[qui]
		stop("ERROR: The Jacobian of the nonlinear part has NA!\nThis concerns the following variables:\n", paste(variables, sep=" ; "))
	}

	#Sauvegarde
	if(nbSave<nbMaxSave){
		savedCoef[[nbSave+1]] = coef
		savedValue[[nbSave+1]] = jacob.mat
		assign(".JC_nbSave", nbSave+1, env)
	} else if(nbMaxSave>1){
		tmp = list()
		tmp[[nbSave]] = coef
		tmp[1:(nbSave-1)] = savedCoef[2:nbSave]
		savedCoef = tmp

		tmp = list()
		tmp[[nbSave]] = jacob.mat
		tmp[1:(nbSave-1)] = savedValue[2:nbSave]
		savedValue = tmp
	} else{
		savedCoef = list(coef)
		savedValue = list(jacob.mat)
	}

	# print(colSums(jacob.mat))

	# cat("computed NL Jacobian:", as.vector(coef), "\n")
	# print(savedCoef)

	assign(".JC_savedCoef", savedCoef, env)
	assign(".JC_savedValue", savedValue, env)
	return(jacob.mat)
}

get_model_null <- function(env, theta.init){
	# I have the closed form of the ll0
	famFuns = get(".famFuns", env)
	family = get(".family", env)
	N = get("nobs", env)
	y = get(".lhs", env)
	verbose = get(".verbose", env)
	ptm = proc.time()

	# one of the elements to be returned
	theta = NULL

	if(family == "poisson"){
		# There is a closed form

		if(".lfactorial" %in% names(env)){
			lfact = get(".lfactorial", env)
		} else {
			# lfactorial(x) == lgamma(x+1)
			# lfact = sum(lfactorial(y))
			lfact = sum(rpar_lgamma(y + 1, env))
			assign(".lfactorial", lfact, env)
		}

		sy = sum(y)
		constant = log(sy / N)
		# loglik =  sy*log(sy) - sy*log(N) - sy - sum(lfactorial(y))
		loglik =  sy*log(sy) - sy*log(N) - sy - lfact
	} else if(family == "gaussian"){
		# there is a closed form
		constant = mean(y)
		ss = sum( (y - constant)**2 )
		sigma = sqrt( ss / N )
		loglik = -1/2/sigma^2*ss - N*log(sigma) - N*log(2*pi)/2
	} else if(family == "logit"){
		# there is a closed form
		sy = sum(y)
		constant = log(sy) - log(N - sy)
		loglik = sy*log(sy) - sy*log(N-sy) - N*log(N) + N*log(N-sy)
	} else if(family=="negbin"){

		if(".lgamma" %in% names(env)){
			lgamm = get(".lgamma", env)
		} else {
			# lgamm = sum(lgamma(y + 1))
			lgamm = sum(rpar_lgamma(y + 1, env))
			assign(".lgamma", lgamm, env)
		}

		sy = sum(y)
		constant = log(sy / N)

		mean_y = mean(y)
		invariant = sum(y*constant) - lgamm

		if(is.null(theta.init)){
			theta.guess = max(mean_y**2 / max((var(y) - mean_y), 1e-4), 0.05)
		} else {
			theta.guess = theta.init
		}

		# I set up a limit of 0.05, because when it is too close to 0, convergence isnt great

		opt <- nlminb(start=theta.guess, objective=famFuns$ll0_theta, y=y, gradient=famFuns$grad0_theta, lower=1e-3, mean_y=mean_y, invariant=invariant, hessian = famFuns$hess0_theta, env=env)

		loglik = -opt$objective
		theta = opt$par
	}

	if(verbose >= 2) cat("Null model in ", (proc.time()-ptm)[3], "s. ", sep ="")

	return(list(loglik=loglik, constant=constant, theta = theta))
}

getGradient = function(jacob.mat, y, mu, exp_mu, env, coef, ...){
	famFuns = get(".famFuns", env)
	ll_dl = famFuns$ll_dl(y, mu, exp_mu, coef=coef, env=env)
	c(crossprod(jacob.mat, ll_dl))
}

getScores = function(jacob.mat, y, mu, exp_mu, env, coef, ...){
	famFuns = get(".famFuns", env)
	isDummy = get("isDummy", env)

	ll_dl = famFuns$ll_dl(y, mu, exp_mu, coef=coef, env=env)
	scores = jacob.mat* ll_dl

	if(isDummy){
		ll_d2 = famFuns$ll_d2(y, mu, exp_mu, coef=coef, env=env)
		dxi_dbeta = deriv_xi(jacob.mat, ll_d2, env, coef)
		scores = scores + dxi_dbeta * ll_dl
	}

	return(as.matrix(scores))
}

getDummies = function(mu, exp_mu, env, coef, final = FALSE){
	# function built to get all the dummy variables
	# We retrieve past dummies (that are likely to be good
	# starting values)
	mu_dummies = get(".savedDummy", env)
	family = get(".family", env)
	eps.cluster = get(".eps.cluster", env)
	verbose = get(".verbose", env)
	if(verbose >= 2) ptm = proc.time()

	#
	# Dynamic precision
	#

	iterCluster = get(".iterCluster", env)
	evolutionLL = get(".evolutionLL", env)
	nobs = get("nobs", env)
	iter = get(".iter", env)
	iterLastPrecisionIncrease = get(".iterLastPrecisionIncrease", env)

	nbIterOne = get(".nbIterOne", env)
	if(iterCluster == 1){
		nbIterOne = nbIterOne + 1
	} else { # we reinitialise
		nbIterOne = 0
	}
	assign(".nbIterOne", nbIterOne, env)

	# nber of times LL almost didn't increase
	nbLowIncrease = get(".nbLowIncrease", env)
	if(evolutionLL/nobs < 1e-8){
		nbLowIncrease = nbLowIncrease + 1
	} else { # we reinitialise
		nbLowIncrease = 0
	}
	assign(".nbLowIncrease", nbLowIncrease, env)

	if(!final && eps.cluster > .Machine$double.eps*10000 && iterCluster==1 && nbIterOne >= 2 && nbLowIncrease >= 2 && (iter - iterLastPrecisionIncrease) >= 3){
		eps.cluster = eps.cluster/10
		if(verbose >= 2) cat("Precision increased to", eps.cluster, "\n")
		assign(".eps.cluster", eps.cluster, env)
		assign(".iterLastPrecisionIncrease", iter, env)

		# If the precision increases, we must also increase the precision of the dummies!
		if(family %in% c("negbin", "logit")){
			assign(".eps.NR", eps.cluster / 100, env)
		}

		# we also set acceleration to on
		assign(".useAcc", TRUE, env)
	} else if(final){
		# we don't need ultra precision for these last dummies
		eps.cluster = eps.cluster * 10**(iterLastPrecisionIncrease != 0)
		if(family %in% c("negbin", "logit")){
			assign(".eps.NR", eps.cluster / 100, env)
		}
	}

	iterMax = get(".itermax.cluster", env)
	nbCluster = get(".nbCluster", env)
	Q = length(nbCluster)

	# whether we use the eponentiation of mu
	useExp_clusterCoef = family %in% c("poisson")
	if(useExp_clusterCoef){
		mu_in = exp_mu * mu_dummies
	} else {
		mu_in = mu + mu_dummies
	}

	#
	# Computing the optimal mu
	#

	index = order(nbCluster, decreasing = TRUE)

	useAcc = get(".useAcc", env)
	carryOn = FALSE

	# Finding the complexity of the problem

	firstRunCluster = get(".firstRunCluster", env)
	if(firstRunCluster && Q >= 3){
		# First iteration: we check if the problem is VERY difficult (for Q = 3+)
		useAcc = TRUE
		assign(".useAcc", TRUE, env)
		res = convergence(coef, mu_in, env, index, iterMax = 15)
		if(res$iter == 15){
			assign(".difficultConvergence", TRUE, env)
			carryOn = TRUE
		}
	} else if(useAcc){
		res = convergence(coef, mu_in, env, index, iterMax)
		if(res$iter <= 2){
			# if almost no iteration => no acceleration next time
			assign(".useAcc", FALSE, env)
		}
	} else {
		res = convergence(coef, mu_in, env, index, iterMax = 15)
		if(res$iter == 15){
			carryOn = TRUE
		}
	}

	if(carryOn){
		# the problem is difficult => acceleration on
		useAcc = TRUE
		assign(".useAcc", TRUE, env)

		res = convergence(coef, res$mu_new, env, index, iterMax)
	}

	mu_new = res$mu_new
	iter = res$iter

	#
	# Retrieving the value of the dummies
	#

	if(useExp_clusterCoef){
		mu_dummies = mu_new / exp_mu
	} else {
		mu_dummies = mu_new - mu
	}

	# Warning messages if necessary:
	if(iter==iterMax) warning("[Getting cluster coefficients] iteration limit reached (", iterMax, ").", call. = FALSE, immediate. = TRUE)

	assign(".iterCluster", iter, env)

	# we save the dummy:
	assign(".savedDummy", mu_dummies, env)

	if(verbose >= 2){
		acc_info = ifelse(useAcc, "+Acc. ", "-Acc. ")
		cat("Cluster Coef.: ", (proc.time()-ptm)[3], "s (", acc_info, "iter:", iter, ")\t", sep = "")
	}

	# we update the flag
	assign(".firstRunCluster", FALSE, env)

	mu_dummies
}

irons_tuck_iteration = function(X, GX, GGX){
	# We compute one iteration

	delta_X = GX - X
	delta_GX = GGX - GX
	delta2_X = delta_GX - delta_X

	GGX - c(delta_GX %*% delta2_X) / c(crossprod(delta2_X)) * delta_GX
}

computeDummies = function(dum, mu, env, coef, sum_y, orderCluster, tableCluster){
	family = get(".family", env)
	famFuns = get(".famFuns", env)
	y = get(".lhs", env)
	eps.NR = get(".eps.NR", env)

	# if family is poisson or gaussian, there is a closed form
	if(family%in%c("poisson", "gaussian")){
		return(famFuns$closedFormDummies(dum, y, mu, env, sum_y, orderCluster, tableCluster))
	}

	#
	# For non closed-form dummies:
	#

	x1 = dichoNR_Cpp(dum, mu, env, coef, sum_y, orderCluster, tableCluster)

	x1
}

dichoNR_Cpp = function(dum, mu, env, coef, sum_y, orderCluster, tableCluster){
	# function that combines dichotomy
	# We are sure that there is a zero, so dichotomy will work 100%

	nb_cases = length(tableCluster)
	y = get(".lhs", env)
	N = length(y)

	# 1st we get the boundaries
	famFuns = get(".famFuns", env)
	eps.NR = get(".eps.NR", env)

	# Get the init, the inf and sup boundaries
	minMax = cpp_conditional_minMax(nb_cases, N, mu, dum, tableCluster)

	# mu est la valeur estimee en les parametres (hors les dummies)
	# plus les mu sont eleves, plus la dummy doit etre petite pour compenser
	borne_inf = famFuns$guessDummy(sum_y, tableCluster, minMax[, 2]) # Value > 0, (on fait moins muMax)
	borne_sup = famFuns$guessDummy(sum_y, tableCluster, minMax[, 1]) # Value < 0

	family = get(".family", env)
	family_nb = switch(family, poisson=1, negbin=2, logit=3, gaussian=4, probit=5)

	theta = coef[".theta"]
	# in the case there are no params to estimate
	if(!is.numeric(theta)) theta = as.numeric(NA)


	#
	# We apply multi cluster if necessary
	#

	isMulticore = get(".isMulticore", env)
	FENmlm_CORES = get(".CORES", env)

	if(!isMulticore){
		res = cpp_DichotomyNR(N = N, K = nb_cases, family = family_nb, theta,
									 epsDicho=eps.NR, lhs=y, mu, borne_inf,
									 borne_sup, orderCluster, tableCluster)
	} else {
		res = cpppar_DichotomyNR(FENmlm_CORES, K = nb_cases, family = family_nb, theta,
										 epsDicho=eps.NR, lhs=y, mu, borne_inf,
										 borne_sup, orderCluster, tableCluster)
	}

	res
}

deriv_xi = function(jacob.mat, ll_d2, env, coef){
	# Derivative of the cluster coefficients

	# data:
	iterMax = get(".itermax.deriv", env)
	nbCluster = get(".nbCluster", env)
	Q = length(nbCluster)

	verbose = get(".verbose", env)
	if(verbose >= 2) ptm = proc.time()

	#
	# initialisation of dxi_dbeta
	#

	if(Q >= 2){
		# We set the initial values for the first run
		if(!".sum_deriv" %in% names(env)){
			# init of the sum of the derivatives => 0
			dxi_dbeta = matrix(0, nrow(jacob.mat), ncol(jacob.mat))
		} else {
			dxi_dbeta = get(".sum_deriv", env)
		}
	} else {
		# no need if only 1, direct solution
		dxi_dbeta = NULL
	}

	#
	# Computing the optimal dxi_dbeta
	#

	index = order(nbCluster, decreasing = TRUE)

	accDeriv = get(".accDeriv", env)
	carryOn = FALSE

	# Finding the complexity of the problem

	firstRunDeriv = get(".firstRunDeriv", env)
	if(firstRunDeriv){
		# set accDeriv: we use information on cluster deriv
		iterCluster = get(".iterCluster", env)
		diffConv = get(".difficultConvergence", env)
		if(iterCluster < 20 & !diffConv){
			accDeriv = FALSE
			assign(".accDeriv", FALSE, env)
		}
	}

	if(firstRunDeriv && accDeriv && Q >= 3){
		# First iteration: we check if the problem is VERY difficult (for Q = 3+)
		assign(".accDeriv", TRUE, env)
		res = dconvergence(dxi_dbeta, jacob.mat, ll_d2, env, index, iterMax = 15)
		if(res$iter == 15){
			assign(".derivDifficultConvergence", TRUE, env)
			carryOn = TRUE
		}
	} else if(accDeriv){
		res = dconvergence(dxi_dbeta, jacob.mat, ll_d2, env, index, iterMax)
		if(res$iter <= 10){
			# if almost no iteration => no acceleration next time
			assign(".accDeriv", FALSE, env)
		}
	} else {
		res = dconvergence(dxi_dbeta, jacob.mat, ll_d2, env, index, iterMax = 50)
		if(res$iter == 50){
			carryOn = TRUE
		}
	}

	if(carryOn){
		# the problem is difficult => acceleration on
		accDeriv = TRUE
		assign(".accDeriv", TRUE, env)
		res = dconvergence(res$dxi_dbeta, jacob.mat, ll_d2, env, index, iterMax)
	}

	dxi_dbeta = res$dxi_dbeta
	iter = res$iter

	if(iter == iterMax) warning("[Getting cluster derivatives] Maximum iterations reached (", iterMax, ").")

	assign(".firstRunDeriv", FALSE, env)
	assign(".sum_deriv", dxi_dbeta, env)

	if(verbose >= 2){
		acc_info = ifelse(accDeriv, "+Acc. ", "-Acc. ")
		cat("  Derivatives: ", (proc.time()-ptm)[3], "s (", acc_info, "iter:", iter, ")\n", sep = "")
	}

	return(dxi_dbeta)
}

deriv_xi_other = function(ll_dx_dother, ll_d2, env, coef){
	# derivative of the dummies wrt an other parameter
	dumMat_cpp = get(".dumMat_cpp", env)
	nbCluster = get(".nbCluster", env)
	dum_all = get(".dum_all", env)
	eps.deriv = get(".eps.deriv", env)
	tableCluster_all = get(".tableCluster", env)
	orderCluster_all = get(".orderCluster", env)
	Q = length(dum_all)
	iterMax = 5000

	if(Q==1){
		dum = dum_all[[1]]
		k = max(dum)
		S_Jmu = rpar_tapply_vsum(k, ll_dx_dother, dum, orderCluster_all[[1]], tableCluster_all[[1]], env)
		S_mu = rpar_tapply_vsum(k, ll_d2, dum, orderCluster_all[[1]], tableCluster_all[[1]], env)
		dxi_dother = - S_Jmu[dum] / S_mu[dum]
	} else {
		# The cpp way:

		N = length(ll_d2)

		# We set the initial values for the first run
		if(!".sum_deriv_other" %in% names(env)){
			init = rep(0, N)
		} else {
			init = get(".sum_deriv_other", env)
		}

		dxi_dother <- RcppPartialDerivative_other(iterMax, Q, N, epsDeriv = eps.deriv, ll_d2, ll_dx_dother, init, dumMat_cpp, nbCluster)

		# we save the values
		assign(".sum_deriv_other", dxi_dother, env)

	}

	as.matrix(dxi_dother)
}

show_vars_limited_width = function(charVect, nbChars = 60){
	# There are different cases

	n = length(charVect)

	if(n==1){
		text = paste0(charVect, ".")
		return(text)
	}

	nb_char_vect = nchar(charVect)
	sumChar = cumsum(nb_char_vect) + (0:(n-1))*2 + 3 + 1

	if(max(sumChar) < nbChars){
		text = paste0(paste0(charVect[-n], collapse = ", "), " and ", charVect[n])
		return(text)
	}

	qui = max(which.max(sumChar > nbChars - 8) - 1, 1)

	nb_left = n - qui

	if(nb_left == 1){
		text = paste0(paste0(charVect[1:qui], collapse = ", "), " and ", nb_left, " other.")
	} else {
		text = paste0(paste0(charVect[1:qui], collapse = ", "), " and ", nb_left, " others.")
	}

	return(text)
}

####
#### Convergence FE ####
####

# In this seciton, we add all the functions used to compute the cluster coefficients

convergence = function(coef, mu_in, env, index, iterMax){
	# computes the new mu wrt the cluster coefficients
	# only for the clusters indicated by index
	# index: an interger vector

	Q = length(index)
	useAcc = get(".useAcc", env)
	diffConv = get(".difficultConvergence", env)

	if(useAcc && diffConv && Q > 2){
		# in case of complex cases: it's more efficient
		# to initialize the first two clusters

		res = convergence(mu_in, env, index[1:2], iterMax)
		mu_in = res$mu_new
	}

	if(Q == 1){
		mu_new = conv_single(coef, mu_in, env, index)
		iter = 1
	} else if(Q >= 2){
		# Dynamic setting of acceleration

		if(!useAcc){
			res = conv_seq(coef, mu_in, env, index, iterMax = iterMax)
		} else if(useAcc){
			res = conv_acc(coef, mu_in, env, index, iterMax = iterMax)
		}

		mu_new = res$mu_new
		iter = res$iter
	}

	# we return a list with: new mu and iterations
	list(mu_new = mu_new, iter = iter)
}

conv_single = function(coef, mu_in, env, index){
	# convergence for a single cluster, the one of index
	# it returns: the new mu (NOT mu_dummies)

	stopifnot(length(index) == 1)
	family = get(".family", env)
	dum_all = get(".dum_all", env)
	sum_y_all = get(".sum_y", env)
	tableCluster_all = get(".tableCluster", env)
	orderCluster_all = get(".orderCluster", env)
	nbCluster = get(".nbCluster", env)
	Q = length(dum_all)


	if(family == "poisson"){
		# For this family, mu_in is in exponential form
		exp_cluster_coef = sum_y_all[[index]] / rpar_tapply_vsum(nbCluster[index], mu_in, dum_all[[index]], orderCluster_all[[index]], tableCluster_all[[index]], env)
		# multiplication since exponential form
		mu_new = mu_in * exp_cluster_coef[dum_all[[index]]]
	} else {
		cluster_coef = computeDummies(dum_all[[index]], mu_in, env, coef, sum_y_all[[index]], orderCluster_all[[index]], tableCluster_all[[index]])
		mu_new = mu_in + cluster_coef[dum_all[[index]]]
	}

	return(mu_new)
}

conv_seq = function(coef, mu_in, env, index, iterMax){
	# convergence of cluster coef without acceleration

	family = get(".family", env)
	Q = length(index)

	if(family == "poisson"){
		if(Q == 2){
			res = conv_seq_pois_2(mu_in, env, index, iterMax)
		} else {
			res = conv_seq_pois_gen(mu_in, env, index, iterMax)
		}
	} else if(family == "gaussian" && Q == 2){
		res = conv_seq_gaus_2(mu_in, env, index, iterMax)
	} else {
		res = conv_seq_genl_gen(coef, mu_in, env, index, iterMax)
	}

	return(res)
}

conv_seq_pois_2 = function(exp_mu_in, env, index, iterMax){
	# computes the new mu for twoclusters in the poisson case
	# this is in exponential form

	# data:
	sum_y_all = get(".sum_y", env)
	dum_all = get(".dum_all", env)
	nbCluster = get(".nbCluster", env)
	eps.cluster = get(".eps.cluster", env)

	# setting up
	setup_poisson_fixedcost(env)
	info = get(".indexOrdered", env)

	i = index[1]
	j = index[2]

	Ab = sum_double_index(info$n[i], info$n[j], info$index[[i]], info$index[[j]], exp_mu_in[info$order])

	# the main loop
	alpha = rep(1, nbCluster[i])
	ca = as.numeric(sum_y_all[[i]])
	cb = as.numeric(sum_y_all[[j]])

	for(iter in 1:iterMax){
		alpha_old = alpha
		alpha = ca / (Ab %m% (cb / (Ab %tm% alpha)))

		diff = max(abs(log(range(alpha_old / alpha))))
		if(diff < eps.cluster) break
	}

	beta = cb / (Ab %tm% alpha)
	alpha = ca / (Ab %m% beta) # we update the biggest vector last
	# final update of mu, with the equilibrium cluster coefficients
	exp_mu_new = exp_mu_in * alpha[dum_all[[i]]] * beta[dum_all[[j]]]

	# return
	list(mu_new = exp_mu_new, iter = iter)
}

conv_seq_pois_gen = function(exp_mu_in, env, index, iterMax){
	# computes the new mu for clusters > 2 in the poisson case
	# this is in exponential form

	# Last update should be the one with the highest nber of cases:
	index = rev(index) # since index is decreasingly sorted

	# data:
	sum_y_all = get(".sum_y", env)
	dum_all = get(".dum_all", env)
	nbCluster = get(".nbCluster", env)
	tableCluster_all = get(".tableCluster", env)
	orderCluster_all = get(".orderCluster", env)
	eps.cluster = get(".eps.cluster", env)

	for(iter in 1:iterMax){
		# cat("----\n")
		exp_mu0 = exp_mu_in

		for(q in index){
			dum = dum_all[[q]]

			# get the dummies
			exp_mu_dum = sum_y_all[[q]] / rpar_tapply_vsum(nbCluster[q], exp_mu_in, dum, orderCluster_all[[q]], tableCluster_all[[q]], env)
			# add them to the stack
			exp_mu_in = exp_mu_in * exp_mu_dum[dum]
		}

		diff <- max(abs(log(range(exp_mu0/exp_mu_in)))) # max error
		if(diff < eps.cluster) break
	}

	# return
	list(mu_new = exp_mu_in, iter = iter)
}

conv_seq_gaus_2 = function(mu_in, env, index, iterMax){

	# data:
	dum_all = get(".dum_all", env)
	eps.cluster = get(".eps.cluster", env)

	# setup
	setup_gaussian_fixedcost(env)

	mat_X = get(".mat_X", env)
	mat_Xx = get(".mat_Xx", env)
	lhs = get(".lhs", env)

	i = index[1]
	j = index[2]

	A = mat_X[[i]]
	B = mat_X[[j]]

	Ab = mat_Xx[[paste0(i, j)]]
	Ba = mat_Xx[[paste0(j, i)]]

	resid = lhs - mu_in

	const_a = A %m% resid
	const_b = B %m% resid
	a_tilde = const_a - (Ab %m% const_b)

	alpha = a_tilde

	for(iter in 1:iterMax){
		alpha_old = alpha
		alpha = as.vector(a_tilde + (Ab %m% (Ba %m% alpha)))

		diff = max(abs(range(alpha_old - alpha)))
		if(diff < eps.cluster) break
	}

	beta = const_b - (Ba %m% alpha)
	alpha = const_a - (Ab %m% beta) # we update the biggest vector last
	# final update of mu, with the equilibrium cluster coefficients
	mu_new = mu_in + alpha[dum_all[[i]]] + beta[dum_all[[j]]]

	# return
	list(mu_new = mu_new, iter = iter)
}

conv_seq_genl_gen = function(coef, mu_in, env, index, iterMax){
	# Sequential case, general

	# Last update should be the one with the highest nber of cases:
	index = rev(index) # since index is decreasingly sorted

	# data:
	sum_y_all = get(".sum_y", env)
	dum_all = get(".dum_all", env)
	nbCluster = get(".nbCluster", env)
	tableCluster_all = get(".tableCluster", env)
	orderCluster_all = get(".orderCluster", env)
	eps.cluster = get(".eps.cluster", env)

	for(iter in 1:iterMax){

		mu0 = mu_in

		for(q in index){
			dum = dum_all[[q]]
			sum_y = sum_y_all[[q]]
			tableCluster = tableCluster_all[[q]]
			orderCluster = orderCluster_all[[q]]
			# get the dummies
			mu_dum = computeDummies(dum, mu_in, env, coef, sum_y, orderCluster, tableCluster)
			# add them to the stack
			mu_in = mu_in + mu_dum[dum]
		}

		if(anyNA(mu_in)) stop("Dummies could not be computed. Unknown error.")

		diff <- max(abs(mu0-mu_in)) # max error
		# cat("iter =", iter, " Diff =", diff, " sum =", sum(abs(mu0-mu_in)), "\n")
		if(diff < eps.cluster) break
	}

	list(mu_new = mu_in, iter = iter)
}

conv_acc = function(coef, mu_in, env, index, iterMax){
	# Convergence with acceleration

	family = get(".family", env)

	if(family == "poisson"){
		res = conv_acc_pois(mu_in, env, index, iterMax)
		if(res$anyNegative){
			# means the IT algorithm does not work
			res = conv_seq(res$mu_new, env, index, iterMax)
		}
	} else {
		res = conv_acc_genl(coef, mu_in, env, index, iterMax)
	}

	return(res)
}

conv_acc_pois = function(exp_mu_in, env, index, iterMax){
	# Convergence + acceleration, poisson case

	# data:
	sum_y_all = get(".sum_y", env)
	dum_all = get(".dum_all", env)
	nbCluster = get(".nbCluster", env)
	tableCluster_all = get(".tableCluster", env)
	orderCluster_all = get(".orderCluster", env)
	eps.cluster = get(".eps.cluster", env)
	Q = length(index)

	# additional data:
	nbCluster_new = nbCluster[index]
	start_cluster = 1 + c(0, cumsum(nbCluster_new))
	end_cluster = cumsum(nbCluster_new)

	if(Q == 2){

		# setup
		setup_poisson_fixedcost(env)
		info = get(".indexOrdered", env)

		i = index[1]
		j = index[2]

		Ab = sum_double_index(info$n[i], info$n[j], info$index[[i]], info$index[[j]], exp_mu_in[info$order])

		# the main loop
		alpha = rep(1, nbCluster[i])
		ca = as.numeric(sum_y_all[[i]])
		cb = as.numeric(sum_y_all[[j]])

		G = function(X){
			X_new = as.vector(ca / (Ab %m% (cb / (Ab %tm% X))))
			X_new
		}

	} else {

		G = function(X){
			# We compute the first item from the other items

			# Note that the order of the data in X_list is different from the real order of the clusters
			# we need to take special care of this specificity

			X_list = list()
			for(i in 1:(Q-1)){
				X_list[[i]] = X[start_cluster[i]:end_cluster[i]]
			}

			# We start from Q because at the moment there is only Q-1 values in X_list
			# starting with it will create its value

			for(i in Q:1){
				q = index[i]

				# getting the value of mu
				exp_mu_current = exp_mu_in
				for(i_other in (1:Q)[-i]){
					# We sum the cluster coefficients of the other clusters
					q_other = index[i_other]
					exp_mu_current = exp_mu_current * X_list[[i_other]][dum_all[[q_other]]]
				}

				# We update the value of X_list
				X_list[[i]] = sum_y_all[[q]] / rpar_tapply_vsum(nbCluster[q], exp_mu_current, dum_all[[q]], orderCluster_all[[q]], tableCluster_all[[q]], env)
			}

			X_new = unlist(X_list[1:(Q-1)])

			# we return them
			X_new
		}

	}


	#
	# The main loop
	#

	anyNegative = FALSE

	X = rep(1, sum(nbCluster[index[1:(Q-1)]]))
	GX = G(X)

	diff = max(abs(X - GX))
	if(diff < eps.cluster){
		iter = 1
	} else {
		for(iter in 1:iterMax){
			GGX = G(GX)
			X_new = irons_tuck_iteration(X, GX, GGX)

			# the values must not be negative!
			if(any(X_new < 0)) {
				X_new[X_new < 0 ] = 1e-12
				anyNegative = TRUE
				break
			}

			GX = G(X_new)

			# Stopping criterion:
			diff = max(abs(X_new - GX))

			# update
			X = X_new

			if(diff < eps.cluster) break
		}
	}

	# since we computed it, we use it:
	X = GX

	#
	# Last iteration -- getting the value of **all** cluster coefficients
	#

	X_list = list()
	for(i in 1:(Q-1)){
		X_list[[i]] = X[start_cluster[i]:end_cluster[i]]
	}

	for(i in Q:1){
		# the real ordered q:
		q = index[i]

		# getting the value of mu
		exp_mu_current = exp_mu_in
		for(i_other in (1:Q)[-i]){
			# We sum the cluster coefficients of the other clusters
			q_other = index[i_other]
			exp_mu_current = exp_mu_current * X_list[[i_other]][dum_all[[q_other]]]
		}

		# We update the value of X_list
		X_list[[i]] = sum_y_all[[q]] / rpar_tapply_vsum(nbCluster[q], exp_mu_current, dum_all[[q]], orderCluster_all[[q]], tableCluster_all[[q]], env)
	}

	# last update of mu based on last update of cluster coef
	q = index[1]
	exp_mu_current = exp_mu_current * X_list[[1]][dum_all[[q]]]

	# the return
	list(mu_new = exp_mu_current, iter = iter, anyNegative = anyNegative)
}

conv_acc_genl = function(coef, mu_in, env, index, iterMax){

	# data:
	sum_y_all = get(".sum_y", env)
	dum_all = get(".dum_all", env)
	nbCluster = get(".nbCluster", env)
	tableCluster_all = get(".tableCluster", env)
	orderCluster_all = get(".orderCluster", env)
	eps.cluster = get(".eps.cluster", env)
	family = get(".family", env)
	Q = length(index)

	# additional data:
	nbCluster_new = nbCluster[index]
	start_cluster = 1 + c(0, cumsum(nbCluster_new))
	end_cluster = cumsum(nbCluster_new)


	# Function defined with some variables local to getDummies function
	if(family == "gaussian" && Q == 2){

		# setup
		setup_gaussian_fixedcost(env)

		mat_X = get(".mat_X", env)
		mat_Xx = get(".mat_Xx", env)
		lhs = get(".lhs", env)

		i = index[1]
		j = index[2]

		A = mat_X[[i]]
		B = mat_X[[j]]

		Ab = mat_Xx[[paste0(i, j)]]
		Ba = mat_Xx[[paste0(j, i)]]

		resid = lhs - mu_in

		const_a = A %m% resid
		const_b = B %m% resid
		a_tilde = const_a - (Ab %m% const_b)

		G = function(X){
			X_new = as.vector(a_tilde + (Ab %m% (Ba %m% X)))
			X_new
		}

	} else {

		G = function(X){
			# We compute the first item from the other items

			# X_list: the list of the cluster coefficients
			# Note that the order of the data in X_list is different from the real order of the clusters
			# we need to take special care of this specificity

			X_list = list()
			for(i in 1:(Q-1)){
				X_list[[i]] = X[start_cluster[i]:end_cluster[i]]
			}

			# We start from Q because at the moment there is only Q-1 values in X_list
			# starting with it will create its value

			for(i in Q:1){
				# the real ordered q:
				q = index[i]

				# getting the value of mu
				mu_current = mu_in
				for(i_other in (1:Q)[-i]){
					# We sum the cluster coefficients of the other clusters
					q_other = index[i_other]
					mu_current = mu_current + X_list[[i_other]][dum_all[[q_other]]]
				}

				# We update the value of X_list
				X_list[[i]] = computeDummies(dum_all[[q]], mu_current, env, coef, sum_y_all[[q]], orderCluster_all[[q]], tableCluster_all[[q]])
			}

			X_new = unlist(X_list[1:(Q-1)])

			# we return them
			X_new
		}
	}

	#
	# The main loop
	#

	X = rep(0, sum(nbCluster[index[1:(Q-1)]]))
	GX = G(X)

	diff = max(abs(X - GX))
	if(diff < eps.cluster){
		iter = 1 # => there is convergence
	} else {
		for(iter in 1:iterMax){
			GGX = G(GX)
			X_new = irons_tuck_iteration(X, GX, GGX)
			GX = G(X_new)

			# Stopping criterion:
			diff = max(abs(X_new - GX))

			# update
			X = X_new

			if(diff < eps.cluster) break
		}
	}

	# since we've made one more iteration, we use it:
	X = GX

	#
	# Last iteration: getting the value of **all** cluster coef.
	#

	X_list = list()
	for(i in 1:(Q-1)){
		X_list[[i]] = X[start_cluster[i]:end_cluster[i]]
	}

	for(i in Q:1){
		q = index[i]

		# getting the value of mu
		mu_current = mu_in
		for(i_other in (1:Q)[-i]){
			# We sum the cluster coefficients of the other clusters
			q_other = index[i_other]
			mu_current = mu_current + X_list[[i_other]][dum_all[[q_other]]]
		}

		# We update the value of X_list
		X_list[[i]] = computeDummies(dum_all[[q]], mu_current, env, coef, sum_y_all[[q]], orderCluster_all[[q]], tableCluster_all[[q]])
	}

	# we add the last item
	q_other = index[1]
	mu_current = mu_current + X_list[[1]][dum_all[[q_other]]]

	# the return
	list(mu_new = mu_current, iter = iter)
}

####
#### Convergence Deriv ####
####

dconvergence = function(dxi_dbeta, jacob.mat, ll_d2, env, index, iterMax){
	# core function to compute the derivatives of the cluster coefficients

	accDeriv = get(".accDeriv", env)
	derivDiffConv = get(".derivDifficultConvergence", env)
	Q = length(index)

	if(accDeriv && derivDiffConv && Q > 2){
		# in case of complex cases: it's more efficient
		# to initialize the first two clusters

		res = dconvergence(dxi_dbeta, jacob.mat, ll_d2, env, index[1:2], iterMax)
		dxi_dbeta = res$dxi_dbeta
	}


	if(Q == 1){

		dxi_dbeta = dconv_single(jacob.mat, ll_d2, env, index)
		iter = 1

	} else {

		# The convergence algorithms
		if(accDeriv){
			res = dconv_acc(dxi_dbeta, jacob.mat, ll_d2, env, index, iterMax)
			dxi_dbeta = res$dxi_dbeta
			iter = res$iter
		} else {
			res = dconv_seq(dxi_dbeta, jacob.mat, ll_d2, env, index, iterMax)
			dxi_dbeta = res$dxi_dbeta
			iter = res$iter
		}
	}

	return(list(dxi_dbeta = dxi_dbeta, iter = iter))
}

dconv_single = function(jacob.mat, ll_d2, env, index){

	# data:
	dum_all = get(".dum_all", env)
	nbCluster = get(".nbCluster", env)
	tableCluster_all = get(".tableCluster", env)
	orderCluster_all = get(".orderCluster", env)

	dum = dum_all[[index]]
	k = nbCluster[index]
	S_Jmu = cpp_tapply_sum(k, jacob.mat*ll_d2, dum)
	S_mu = rpar_tapply_vsum(k, ll_d2, dum, orderCluster_all[[index]], tableCluster_all[[index]], env)
	dxi_dbeta = - S_Jmu[dum, ] / S_mu[dum]

	return(dxi_dbeta)
}

dconv_seq = function(deriv_init, jacob.mat, ll_d2, env, index, iterMax){
	# deriv_init: past value of dxi_dbeta

	# data:
	Q = length(index)
	dumMat_cpp = get(".dumMat_cpp", env)
	nbCluster = get(".nbCluster", env)
	family = get(".family", env)
	eps.deriv = get(".eps.deriv", env)

	N = nrow(jacob.mat)
	K = ncol(jacob.mat)

	# we reorder the index: largest groups at the end
	index = rev(index)
	nbCluster = nbCluster[index]
	dumMat_cpp = dumMat_cpp[, index, drop = FALSE]

	if(family == "gaussian"){

		res = RcppPartialDerivative_gaussian_new(iterMax, Q, N, K, epsDeriv = eps.deriv, jacob.mat, deriv_init, dumMat_cpp, nbCluster)

	} else {

		res = RcppPartialDerivative_new(iterMax, Q, N, K, epsDeriv = eps.deriv, ll_d2, jacob.mat, deriv_init, dumMat_cpp, nbCluster)

	}

	return(list(dxi_dbeta = res$dxi_dbeta, iter = res$iter))
}

dconv_acc = function(deriv_init, jacob.mat, ll_d2, env, index, iterMax){

	# data:
	dum_all = get(".dum_all", env)
	nbCluster = get(".nbCluster", env)
	tableCluster_all = get(".tableCluster", env)
	orderCluster_all = get(".orderCluster", env)
	eps.deriv = get(".eps.deriv", env)
	Q = length(index)

	nobs = length(ll_d2)

	# The weights:
	weight_all = list()

	for(q in 1:Q){
		sum_lld2 = -1/rpar_tapply_vsum(nbCluster[q], ll_d2, dum_all[[q]], orderCluster_all[[q]], tableCluster_all[[q]], env)
		weight_all[[q]] = ll_d2 * sum_lld2[dum_all[[q]]]
	}

	# additional data:
	nbCluster_new = nbCluster[index]
	start_cluster = 1 + c(0, cumsum(nbCluster_new))
	end_cluster = cumsum(nbCluster_new)

	newDerivFlag = FALSE

	if(newDerivFlag && Q == 2){

		# setup
		setup_deriv_fixedcost(env)
		order_both = get(".order_ij_ji", env)

		i = index[1]
		j = index[2]

		dum_1 = dum_all[[i]]
		dum_2 = dum_all[[j]]

		n1 = nbCluster[i]
		n2 = nbCluster[j]

		# les matrices types normalisees
		coefmat_i = weight_all[[i]]
		A = list(n_i = n1, n_j = nobs, index_i = dum_1 - 1L, index_j = (1:nobs) - 1L, coefmat = coefmat_i)
		coefmat_j = weight_all[[j]]
		B = list(n_i = n2, n_j = nobs, index_i = dum_2 - 1L, index_j = (1:nobs) - 1L, coefmat = coefmat_j)

		# Other matrices
		order_ij = order_both$order_ij
		Ab = sum_double_index(n1, n2, dum_1[order_ij], dum_2[order_ij], coefmat_i[order_ij])

		order_ji = order_both$order_ji
		Ba = sum_double_index(n2, n1, dum_2[order_ji], dum_1[order_ji], coefmat_j[order_ji])

		G = function(X){
			X_new = const_a + (Ab %m% (Ba %m% X))
			X_new
		}

	} else {
		G = function(X){
			# Jk doit etre egal a la valeur de "Jk" + "sum deriv init"

			X_list = list()
			for(i in 1:(Q-1)){
				X_list[[i]] = X[start_cluster[i]:end_cluster[i]]
			}

			# We start from Q because at the moment there is only Q-1 values in X_list
			# starting with it will create its value

			for(i in Q:1){
				q = index[i]

				# getting the value of the sum of derivatives
				Jk_sum_deriv = Jk
				for(i_other in (1:Q)[-i]){
					# We sum the cluster derivatives of the other clusters
					q_other = index[i_other]
					Jk_sum_deriv = Jk_sum_deriv + X_list[[i_other]][dum_all[[q_other]]]
				}

				# We update the value of X_list
				X_list[[i]] = rpar_tapply_vsum(nbCluster[q], weight_all[[q]]*Jk_sum_deriv, dum_all[[q]], orderCluster_all[[q]], tableCluster_all[[q]], env)
			}

			X_new = unlist(X_list[1:(Q-1)])

			# we return them
			X_new
		}
	}

	# saving the number of iterations
	max_iter = 1

	dxi_dbeta_list = list()
	K = ncol(jacob.mat)
	for(k in 1:K){
		Jk = jacob.mat[, k] + deriv_init[, k]

		if(newDerivFlag && Q == 2){
			# set up the constants:
			a = A %m% Jk
			b = B %m% Jk
			const_a = a + (Ab %m% b)
		}

		X = rep(0, sum(nbCluster[index[1:(Q-1)]]))
		GX = G(X)

		diff <- max(abs(X-GX))
		if(diff < eps.deriv){
			iter = 1 # convergence
		} else {
			for(iter in 1:iterMax){
				GGX = G(GX)
				X_new = irons_tuck_iteration(X, GX, GGX)
				GX = G(X_new)

				diff <- max(abs(X_new-GX))

				# update
				X = X_new

				if(diff < eps.deriv) break
			}
		}

		# we use the last iteration
		X = GX

		if(iter > max_iter){
			max_iter = iter
		}

		#
		# Last iteration, we end with the largest cluster
		#


		X_list = list()
		for(i in 1:(Q-1)){
			X_list[[i]] = X[start_cluster[i]:end_cluster[i]]
		}

		for(i in Q:1){
			q = index[i]

			# getting the value of the sum of derivatives
			Jk_sum_deriv = Jk
			for(i_other in (1:Q)[-i]){
				# We sum the cluster derivatives of the other clusters
				q_other = index[i_other]
				Jk_sum_deriv = Jk_sum_deriv + X_list[[i_other]][dum_all[[q_other]]]
			}

			# We update the value of X_list
			X_list[[i]] = rpar_tapply_vsum(nbCluster[q], weight_all[[q]]*Jk_sum_deriv, dum_all[[q]], orderCluster_all[[q]], tableCluster_all[[q]], env)
		}

		q_other = index[1]
		dxi_dbeta_list[[k]] = Jk_sum_deriv + X_list[[1]][dum_all[[q_other]]] - jacob.mat[, k]

	}

	dxi_dbeta = do.call("cbind", dxi_dbeta_list)

	# we save the values
	return(list(dxi_dbeta = dxi_dbeta, iter = max_iter))
}

####
#### Misc FE ####
####

setup_poisson_fixedcost = function(env){

	# We set up only one
	if(".indexOrdered" %in% names(env)){
		return(NULL)
	}

	ptm = proc.time()

	dum_all = get(".dum_all",env)

	dum_A = as.integer(dum_all[[1]])
	dum_B = as.integer(dum_all[[2]])

	myOrder = order(dum_A, dum_B)
	index_i = dum_A[myOrder]
	index_j = dum_B[myOrder]

	res = list(n = c(max(dum_A), max(dum_B)), index = list(index_i, index_j), order = myOrder)

	assign(".indexOrdered", res, env)

	verbose = get(".verbose", env)
	if(verbose >= 2) cat("Poisson fixed-cost setup: ", (proc.time()-ptm)[3], "s\n", sep = "")
}


setup_gaussian_fixedcost = function(env){

	# We set up only one
	if(".mat_X" %in% names(env)){
		return(NULL)
	}

	ptm = proc.time()

	lhs = get(".lhs", env)
	tableCluster_all = get(".tableCluster", env)
	dum_all = get(".dum_all", env)

	n = length(lhs)
	mat_X = list()
	mat_Xx = list()

	# NEW

	dum_1 = as.integer(dum_all[[1]])
	dum_2 = as.integer(dum_all[[2]])

	n1 = max(dum_1)
	n2 = max(dum_2)

	# les matrices normalisees: A, B
	coefmat_i = 1/(tableCluster_all[[1]][dum_1])
	mat_X[[1]] = list(n_i = n1, n_j = n, index_i = dum_1 - 1L, index_j = (1:n) - 1L, coefmat = coefmat_i)

	coefmat_j = 1/(tableCluster_all[[2]][dum_2])
	mat_X[[2]] = list(n_i = n2, n_j = n, index_i = dum_2 - 1L, index_j = (1:n) - 1L, coefmat = coefmat_j)

	# les matrices: Ab, Ba
	order_ij = order(dum_1, dum_2)
	mat_Xx[["12"]] = sum_double_index(n1, n2, dum_1[order_ij], dum_2[order_ij], coefmat_i[order_ij])

	order_ji = order(dum_2, dum_1)
	mat_Xx[["21"]] = sum_double_index(n2, n1, dum_2[order_ji], dum_1[order_ji], coefmat_j[order_ji])

	assign(".mat_X", mat_X, env)
	assign(".mat_Xx", mat_Xx, env)
	assign(".order_ij_ji", list(order_ij = order_ij, order_ji = order_ji), env)

	verbose = get(".verbose", env)
	if(verbose >= 2) cat("Gaussian fixed-cost setup: ", (proc.time()-ptm)[3], "s\n", sep = "")
}

setup_deriv_fixedcost = function(env){
	# We set up only one
	if(".order_ij_ji" %in% names(env)){
		return(NULL)
	}

	ptm = proc.time()

	dum_all = get(".dum_all", env)

	dum_1 = dum_all[[1]]
	dum_2 = dum_all[[2]]

	# Matrices Ab et Ba
	order_ij = order(dum_1, dum_2)
	order_ji = order(dum_2, dum_1)

	assign(".order_ij_ji", list(order_ij = order_ij, order_ji = order_ji), env)

	verbose = get(".verbose", env)
	if(verbose >= 2) cat("Deriv. fixed-cost setup: ", (proc.time()-ptm)[3], "s\n", sep = "")
}

# matrix multiply
"%m%" = function(mat, x){
	mmult(mat$n_i, mat$index_i, mat$index_j, mat$coefmat, x)
}

"%tm%" = function(mat, x){
	mmult(mat$n_j, mat$index_j, mat$index_i, mat$coefmat, x)
}

####
#### Parallel Functions ####
####

# In this section, we create all the functions that will be parallelized

rpar_exp = function(x, env){
	# fast exponentiation
	isMulticore = get(".isMulticore", env)
	FENmlm_CORES = get(".CORES", env)

	if(!isMulticore){
		# simple exponentiation
		return(exp(x))
	} else {
		# parallelized one
		return(cpppar_exp(x, FENmlm_CORES))
	}

}

rpar_log = function(x, env){
	# fast log
	isMulticore = get(".isMulticore", env)
	FENmlm_CORES = get(".CORES", env)

	if(!isMulticore){
		# simple log
		return(log(x))
	} else {
		# parallelized one
		return(cpppar_log(x, FENmlm_CORES))
	}

}

rpar_lgamma = function(x, env){
	# fast lgamma
	isMulticore = get(".isMulticore", env)
	FENmlm_CORES = get(".CORES", env)

	if(!isMulticore){
		# lgamma via cpp is faster
		return(cpp_lgamma(x))
	} else {
		# parallelized one
		return(cpppar_lgamma(x, FENmlm_CORES))
	}

}

rpar_digamma = function(x, env){
	isMulticore = get(".isMulticore", env)
	FENmlm_CORES = get(".CORES", env)

	if(!isMulticore){
		# digamma via cpp is as fast => no need
		return(digamma(x))
	} else {
		# parallelized one
		return(cpppar_digamma(x, FENmlm_CORES))
	}

}

rpar_trigamma = function(x, env){
	isMulticore = get(".isMulticore", env)
	FENmlm_CORES = get(".CORES", env)

	if(!isMulticore){
		# trigamma via cpp is as fast => no need
		return(trigamma(x))
	} else {
		# parallelized one
		return(cpppar_trigamma(x, FENmlm_CORES))
	}

}

rpar_log_a_exp = function(a, mu, exp_mu, env){
	# compute log_a_exp in a fast way
	isMulticore = get(".isMulticore", env)
	FENmlm_CORES = get(".CORES", env)

	if(!isMulticore){
		# cpp is faster
		return(cpp_log_a_exp(a, mu, exp_mu))
	} else {
		# parallelized one
		return(cpppar_log_a_exp(FENmlm_CORES, a, mu, exp_mu))
	}
}

rpar_tapply_vsum = function(K, x, dum, obsCluster, tableCluster, env){
	isMulticore = get(".isMulticore", env)
	FENmlm_CORES = get(".CORES", env)

	if(isMulticore){
		res <- cpppar_tapply_vsum(FENmlm_CORES, K, x, obsCluster, tableCluster)
	} else {
		res <- cpp_tapply_vsum(K, x, dum)
	}

	res
}

