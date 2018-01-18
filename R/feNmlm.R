# Commands to genereate the help files:
# file.sources = list.files("R", pattern="*.R$", full.names=TRUE, ignore.case=TRUE)
# sapply(file.sources, source, .GlobalEnv)
# roxygen2::roxygenise(roclets = "rd")

# TODO:
# Add only cluster with negbin (needs to create the functions to maximize it)
# update res2table and res2tex for models with onlyCluster
# change name: id_dummies and dummies to id_clusters and sumClusterCoefs (or a better name)


#' Fixed effects maximum likelihood models
#'
#' This function estimates maximum likelihood models (e.g., Poisson or Logit) and is efficient to handle any number of fixed effects (i.e. cluster variables). It further allows for nonlinear in parameters right hand sides.
#'
#' @param fml A formula. This formula gives the linear formula to be estimated (it is similar to a \code{lm} formula), for example: \code{fml = z~x+y}. To include cluster variables, you can 1) either insert them in this formula using a pipe (e.g. \code{fml = z~x+y|cluster1+cluster2}), or 2) either use the argment \code{cluster}. You can add a non-linear element in this formula by using the argment \code{NL.fml}. If you want to estimate only a non-linear formula without even the intercept, you can use \code{fml = z~0} in combination with \code{NL.fml}.
#' @param NL.fml A formula. If provided, this formula represents the non-linear part of the right hand side (RHS). Note that contrary to the \code{fml} argument, the coefficients must explicitely appear in this formula. For instance, it can be \code{~a*log(b*x + c*x^3)}, where \code{a}, \code{b}, and \code{c} are the coefficients to be estimated. Note that only the RHS of the formula is to be provided, and NOT the left hand side.
#' @param data A data.frame containing the necessary variables to run the model. The variables of the non-linear right hand side of the formula are identified with this \code{data.frame} names. Note that no \code{NA} is allowed in the variables to be used in the estimation.
#' @param family Character scalar. It should provide the family. The possible values are "poisson" (Poisson model with log-link, the default), "negbin" (Negative Binomial model with log-link), "logit" (LOGIT model with log-link), "gaussian" (Gaussian model).
#' @param cluster Character vector. The name/s of a/some variable/s within the dataset to be used as clusters. These variables should contain the identifier of each observation (e.g., think of it as a panel identifier).
#' @param start A list. Starting values for the non-linear parameters. ALL the parameters are to be named and given a staring value. Example: \code{start=list(a=1,b=5,c=0)}. Though, there is an exception: if all parameters are to be given the same starting value, you can use the argument \code{start.init}.
#' @param lower A list. The lower bound for each of the non-linear parameters that requires one. Example: \code{lower=list(b=0,c=0)}. Beware, if the estimated parameter is at his lower bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.
#' @param upper A list. The upper bound for each of the non-linear parameters that requires one. Example: \code{upper=list(a=10,c=50)}. Beware, if the estimated parameter is at his upper bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.
#' @param env An environment. You can provide an environement in which the non-linear part will be evaluated. (May be useful for some particular non-linear functions.)
#' @param start.init Numeric scalar. If the argument \code{start} is not provided, or only partially filled (i.e. there remain non-linear parameters with no starting value), then the starting value of all remaining non-linear parameters is set to \code{start.init}.
#' @param offset A formula. An offset can be added to the estimation. It should be a formula of the form (for example) ~0.5*x**2. This offset is linearily added to the elements of the main formula 'fml'. Note that when using the argument 'NL.fml', you can directly add the offset there.
#' @param nl.gradient A formula. The user can prodide a function that computes the gradient of the non-linear part. The formula should be of the form \code{~f0(a1,x1,a2,a2)}. The important point is that it should be able to be evaluated by: \code{eval(nl.gradient[[2]], env)} where \code{env} is the working environment of the algorithm (which contains all variables and parameters). The function should return a list or a data.frame whose names are the non-linear parameters.
#' @param linear.start Numeric named vector. The starting values of the linear part. Note that you can
#' @param jacobian.method Character scalar. Provides the method used to numerically compute the jacobian of the non-linear part. Can be either \code{"simple"} or \code{"Richardson"}. Default is \code{"simple"}. See the help of \code{\link[numDeriv]{jacobian}} for more information.
#' @param useHessian Logical. Should the Hessian be computed in the optimization stage? Default is \code{TRUE}.
#' @param opt.control List of elements to be passed to the optimization method (\code{\link[stats]{nlminb}}.
#' @param debug Logical. If \code{TRUE} then the log-likelihood as well as all parameters are printed at each iteration. Default is \code{FALSE}.
#' @param theta.init Positive numeric scalar. The starting value of the dispersion parameter if \code{family="negbin"}. By default, the algorithm uses as a starting value the theta obtained from the model with only the intercept.
#' @param noWarning Logical, default is \code{FALSE}. Should the warnings be displayed?
#' @param ... Not currently used.
#'
#' @return
#' An \code{femlm} object.
#' \item{coef}{The coefficients.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
#' \item{loglik}{The loglikelihood.}
#' \item{iterations}{Number of iterations of the algorithm.}
#' \item{n}{The number of observations.}
#' \item{k}{The number of parameters of the model.}
#' \item{call}{The call.}
#' \item{nonlinear.fml}{The nonlinear formula of the call. It also contains the dependent variable.}
#' \item{linear.formula}{The linear formula of the call.}
#' \item{ll_null}{Log-likelihood of the null model (i.e. with the intercept only).}
#' \item{pseudo_r2}{The adjusted pseudo R2.}
#' \item{naive.r2}{The R2 as if the expected predictor was the linear predictor in OLS.}
#' \item{message}{The convergence message from the optimization procedures.}
#' \item{sq.cor}{Squared correlation between the dependent variable and its expected value as given by the optimization.}
#' \item{expected.predictor}{The expected predictor is the expected value of the dependent variable.}
#' \item{cov.unscaled}{The variance covariance matrix of the parameters.}
#' \item{se}{The standard-error of the parameters.}
#'
#' @seealso
#' See also \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
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
#' # 2) Log-Log Gaussian estimation
#' est_gaus = femlm(log(Euros+1) ~ log(dist_km)|Origin+Destination+Product, trade, family="gaussian")
#'
#' # 3) Negative Binomial estimation
#' est_nb = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade, family="negbin")
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
#' # Non-linear examples
#' #
#'
#' # Generating data for a simple example
#' n = 100
#' x = rnorm(n, 1, 5)**2
#' y = rnorm(n, -1, 5)**2
#' z = rpois(n, x*y) + rpois(n, 2)
#' base = data.frame(x, y, z)
#'
#' # Comparing the results of a 'linear' function using a 'non-linear' call
#' est0L = femlm(z~log(x)+log(y), base)
#' est0NL = femlm(z~1, base, NL.fml = ~a*log(x)+b*log(y), start = list(a=0, b=0))
#' # we compare the estimates with the function res2table
#' res2table(est0L, est0NL)
#'
#' # Generating a non-linear relation
#' z2 = rpois(n, x + y) + rpois(n, 1)
#' base$z2 = z2
#'
#' # Using a non-linear form
#' est1NL = femlm(z2~0, base, NL.fml = ~log(a*x + b*y), start = list(a=1, b=2), lower = list(a=0, b=0))
#' # we can't estimate this relation linearily
#' # => closest we can do:
#' est1L = femlm(z2~log(x)+log(y), base)
#'
#' res2table(est1L, est1NL)
#'
#' # Using a custom Jacobian for the function log(a*x + b*y)
#' myGrad = function(a,x,b,y){
#' 	# Custom Jacobian
#' 	s = a*x+b*y
#' 	data.frame(a = x/s, b = y/s)
#' }
#'
#' est1NL_grad = femlm(z2~0, base, NL.fml = ~log(a*x + b*y), start = list(a=1,b=2),
#'                      nl.gradient = ~myGrad(a,x,b,y))
#'
#'
femlm <- function(fml, data, family=c("poisson", "negbin", "logit", "gaussian"), NL.fml, cluster, start, lower, upper, env, start.init, offset, nl.gradient, linear.start=0, jacobian.method=c("simple", "Richardson"), useHessian=TRUE, opt.control=list(), debug=FALSE, theta.init, noWarning = FALSE, ...){

	# use of the conjugate gradient in the gaussian case to get
	# the cluster coefficients
	useCG = FALSE

	# opt_method = c("nlminb", "optim")
	# opt_method <- match.arg(opt_method)
	opt_method = "nlminb"
	optim.method="BFGS" # deprec

	jacobian.method <- match.arg(jacobian.method)
	family = match.arg(family)

	# Some settings (too complicated to be tweaked by the user)
	# Nber of evaluations of the NL part to be kept in memory
	# Default keeps the two last evaluations
	NLsave = 2

	# PRECISION
	dots = list(...)

	# I initially called the cluster dummies... I keep it for compatibility
	if(missing(cluster) && "dummy" %in% names(dots)) cluster = dots$dummy
	if("linear.fml" %in% names(dots)) stop("Argument 'linear.fml' is deprecated, now use 'fml' in combination with 'NL.fml'.")

	#
	# The clusters => controls + setup
	fml_char = as.character(fml)[3]
	n_pipes = length(strsplit(fml_char, "|", fixed=TRUE)[[1]]) - 1
	if(n_pipes >= 2) stop("The argument 'fml' cannot contain more than one '|'.")
	extract_fml = extractCluster(fml)
	if(!is.null(extract_fml$cluster)){
		if(missing(cluster) || length(cluster) == 0){
			cluster = extract_fml$cluster
			fml = extract_fml$fml
		} else {
			stop("To add cluster variables: either include them in argument 'fml' using a pipe ('|'), either use the argument 'cluster'. You cannot use both!")
		}
	}

	# Other custom parameters
	eps.cluster = ifelse(is.null(dots$eps.cluster), 1e-5, dots$eps.cluster)
	eps.NR = ifelse(is.null(dots$eps.NR), eps.cluster/100, dots$eps.NR)
	eps.deriv = ifelse(is.null(dots$eps.deriv), 1e-4, dots$eps.deriv)
	# if it is null => OK
	d.hessian = dots$d.hessian

	famFuns = switch(family,
						 poisson = ml_poisson(),
						 negbin = ml_negbin(),
						 logit = ml_logit(),
						 tobit = ml_tobit(),
						 probit = ml_probit(),
						 gaussian = ml_gaussian())

	stopifnot(class(fml)=="formula")
	if(length(fml)!=3) stop("The formula must be two sided.\nEG: a~exp(b/x), or a~0 if there is no nonlinear part.")
	call = match.call()
	dataNames = names(data)

	# The LHS must contain only values in the DF
	namesLHS = all.vars(fml[[2]])
	if(!all(namesLHS%in%dataNames)) stop("Some elements on the LHS of the formula are not in the dataset:\n", paste0(namesLHS[!namesLHS%in%dataNames], collapse=", "))

	# Now the nonlinear part:

	if(!missing(NL.fml)){
		isNonLinear = TRUE
		nl.call = NL.fml[[length(NL.fml)]]
		# allnames = all.vars(fml[[3]])
		allnames = all.vars(nl.call)
		nonlinear.params = allnames[!allnames %in% dataNames]
		nonlinear.varnames = allnames[allnames %in% dataNames]

		if(length(nonlinear.params) == 0){
			warning("As there is no parameter to estimate in argument 'NL.fml', this argument is ignored.\nIf you want to add an offset, use argument 'offset'.")
		}

	} else {
		isNonLinear = FALSE
		nl.call = 0
		allnames = nonlinear.params = nonlinear.varnames = character(0)
	}

	# The conversion of the data (due to data.table)
	if("data.table" %in% class(data)){
		class(data) = "data.frame"
	}

	# The dependent variable: lhs==left_hand_side
	lhs = as.vector(eval(fml[[2]], data))

	# creation de l'environnement
	if(missing(env)) env <- new.env()
	else stopifnot(class(env)=="environment")

	#
	# First check
	#

	# NA are not allowed !!!
	if(anyNA(lhs)) stop("The left hand side of the fomula has NA values. Please provide data without NA.")
	if(family%in%c("poisson", "negbin") & any(lhs<0)) stop("Negative values of the dependant variable \nare not allowed for the \"", family, "\" family.", sep="")
	if(family%in%c("logit") & !all(lhs==0 | lhs==1)) stop("The dependant variable has values different from 0 or 1.\nThis is not allowed with the \"logit\" family.")

	# Add checks for the clusters

	#
	# Controls and setting of the linear part:
	#

	isLinear = FALSE
	linear.varnames = all.vars(fml[[3]])

	if(length(linear.varnames)>0 || attr(terms(fml), "intercept")==1){
		isLinear = TRUE
		linear.fml = fml
	}

	# if(!missing(linear.fml)){
	# 	isLinear <- TRUE
	if(isLinear){
		# if(class(linear.fml)!="formula" || length(linear.fml)!=2) stop("'linear.fml' must be a formula like, for ex., ~x1+x2-1")
		# linear.varnames <- all.vars(linear.fml)
		if(!all(linear.varnames%in%dataNames)) stop(paste("In 'fml', some variables are not in the data:\n", paste(linear.varnames[!linear.varnames%in%dataNames], collapse=', '), ".", sep=""))
		if(!missing(cluster) && length(cluster)!=0){
			#if dummies are provided, we make sure there is an
			#intercept so that factors can be handled properly
			linear.fml = update(linear.fml, ~.+1)
		}
		linear.mat = stats::model.matrix(linear.fml, data)
		linear.params <- colnames(linear.mat)
		N_linear <- length(linear.params)
		if(anyNA(linear.mat)){
			quiNA = apply(linear.mat, 2, anyNA)
			whoIsNA = linear.params[quiNA]

			text = show_vars_limited_width(whoIsNA)

			stop("Evaluation of the linear part returns NA. NAs are not supported, please remove them before running this function. FYI the variables with NAs are:\n", text)
		}
		if(!is.numeric(linear.start)) stop("'linear.start' must be numeric!")
	} 	else linear.params <- linear.start <- linear.varnames <- NULL

	params <- c(nonlinear.params, linear.params)
	lparams <- length(params)
	varnames <- c(nonlinear.varnames, linear.varnames)

	# Attention les parametres non lineaires peuvent etre vides
	if(length(nonlinear.params)==0) isNL = FALSE
	else isNL = TRUE

	# Control for NAs
	if(anyNA(data[, varnames])){
		varWithNA = varnames[which(apply(data[, varnames, FALSE], 2, anyNA))]
		text = show_vars_limited_width(varWithNA)
		stop("Some variables used for this estimation contain NA. NAs are not supported, please remove them first.\nFYI, the variables are:\n", text, call. = FALSE)
	}

	#
	# Handling Clusters ####
	#

	isDummy = FALSE
	if(!missing(cluster) && length(cluster)!=0){
		# TODO:
		# - mettre un meilleur controle des classes pures (mettre un while, ici ne suffit pas)
		# - ajouter la possibilite de mettre une "reference" (ie oter une dummy)

		isDummy = TRUE
		if(class(cluster)!="character" | any(!cluster%in%names(data))){
			var_problem = cluster[!cluster%in%names(data)]
			stop("The argument 'cluster' must be a variable name!\nCluster(s) not in the data: ", paste0(var_problem, collapse = ", "), ".")
		}
 		#qui = which(sapply(cluster, function(x) is.factor(class(data[[x]]))))
 		#if(length(qui)>0) stop("The variable(s) defining the cluster(s) must be a factor!\nIt concerns:", paste(cluster[qui], collapse=", "))

		Q = length(cluster)
		dum_all = dum_names = list()
		sum_y_all = obs_per_cluster_all = list()
		obs2remove = c()
		for(i in 1:Q){
			dum_raw = data[[cluster[i]]]
			# in order to avoid "unclassed" values > real nber of classes: we re-factor the cluster

			# DEPREC
			# dum_names[[i]] = thisNames = sort(unique(dum_raw))
			# if(class(dum_raw) == "factor"){
			# 	dum = unclass(as.factor(unclass(dum_raw)))
			# } else {
			# 	dum = unclass(as.factor(dum_raw))
			# }
			# END DEPREC

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
				dummyOmises = thisNames[qui] # not currently used
				obs2remove = unique(c(obs2remove, which(dum %in% qui)))
			}
		}

		# We remove the problems
		if(length(obs2remove)>0){
			data = data[-obs2remove, ]

			# We recreate the linear matrix and the LHS
			if(isLinear) linear.mat = stats::model.matrix(linear.fml, data)
			lhs = eval(fml[[2]], data)

			# Then we recreate the dummies
			for(i in 1:Q){
				dum_raw = data[[cluster[i]]]

				# DEPREC
				# dum_names[[i]] = sort(unique(dum_raw))
				# if(class(dum_raw) == "factor"){
				# 	dum = unclass(as.factor(unclass(dum_raw)))
				# } else {
				# 	dum = unclass(as.factor(dum_raw))
				# }
				# END DEPREC

				dum_names[[i]] = getItems(dum_raw)
				dum = quickUnclassFactor(dum_raw)

				dum_all[[i]] = dum
				k = length(dum_names[[i]])

				# We also recreate these values
				sum_y_all[[i]] = cpp_tapply_vsum(k, lhs, dum)
				obs_per_cluster_all[[i]] = cpp_table(k, dum)

			}

			# Then the warning message
			if(!noWarning) warning(length(dummyOmises), " clusters (", length(obs2remove), " observations) removed because of only ", ifelse(family=="logit", "zero/one", "zero"), " outcomes.", call. = FALSE, immediate. = TRUE)
		}

		# We compute two values that will be useful to compute the derivative wrt clusters
		# Trick to go faster
		dumMat_cpp = matrix(unlist(dum_all), ncol = Q) - 1

		nbCluster = sapply(dum_all, function(x) length(unique(x)))

		# If there is a linear intercept, we withdraw it
		# We drop the intercept:
		if("(Intercept)" %in% colnames(linear.mat)){
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
				N_linear <- length(linear.params)
				params <- c(nonlinear.params, linear.params)
				lparams <- length(params)
				varnames <- c(nonlinear.varnames, linear.varnames)
			}
		}
	} else {
		# There is no cluster
		Q = 0
	}

	# browser()

	#
	# Checks for MONKEY TEST
	#

	if(lparams==0) stop("No parameter to be estimated.")
	if(!is.logical(useHessian)) stop("'useHessian' must be of type 'logical'!")

	#
	# Controls: The non linear part
	#

	if(isNL){
		if(missing(start.init)){
			if(missing(start)) stop("There must be starting values.")
			if(typeof(start)!="list") stop("start must be a list.")
			if(any(!names(start) %in% params)) stop(paste("Some parameters in 'start' are not in the formula:\n", paste(names(start)[!names(start) %in% params], collapse=", "), ".", sep=""))
			if(any(!nonlinear.params %in% names(start))) stop(paste("Hey, some parameters have no starting values:\n", paste(nonlinear.params[!nonlinear.params%in%names(start)], collapse=", "), ".", sep=""))
		}
		else{
			if(length(start.init)>1) stop("start.init musn't be a vector.")
			if(class(start.init)!="numeric") stop("start.init must be numeric!")
			if(!is.finite(start.init)) stop("Infinites values as starting values, you must be kidding me...")

			if(missing(start)){
				start <- list()
				start[nonlinear.params] <- start.init
			}
			else{
				if(typeof(start)!="list") stop("start must be a list.")
				if(any(!names(start) %in% params)) stop(paste("Some parameters in 'start' are not in the formula:\n", paste(names(start)[!names(start) %in% params], collapse=", "), ".", sep=""))

				missing.params <- nonlinear.params[!nonlinear.params%in%names(start)]
				start[missing.params] <- start.init
			}
		}
	} else start <- list()

	#
	# Controls: The upper and lower limits
	#

	if(!missing(lower)){
		if(typeof(lower)!="list") stop("'lower' MUST be a list.")
		if(any(!names(lower)%in%params)){
			text <- paste("Hey, some parameters in 'lower' are not in the formula:\n", paste(names(lower)[!names(lower)%in%params], collapse=", "), ".", sep="")
			stop(text)
		}
	}
	if(!missing(upper)){
		if(typeof(upper)!="list") stop("'upper' MUST be a list.")
		if(any(!names(upper)%in%params)){
			text <- paste("Hey, some parameters in 'upper' are not in the formula:\n", paste(names(upper)[!names(upper)%in%params], collapse=", "), ".", sep="")
			stop(text)
		}
	}

	# Now setting upper and lower

	if(!missing(lower)){
		lower[params[!params%in%names(lower)]] <- -Inf
		lower <- unlist(lower[params])
	}	else {
		#TODO B
		lower <- rep(-Inf, lparams)
		names(lower) <- params
	}
	if(!missing(upper)){
		upper[params[!params%in%names(upper)]] <- Inf
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
		isGradient=TRUE
		if(class(nl.gradient)!="formula" | length(nl.gradient)==3) stop("'nl.gradient' must be a formula like, for ex., ~f0(a1, x1, a2, x2). f0 giving the gradient.")
	} else {
		isGradient=FALSE
	}

	if(!is.null(d.hessian)){
		hessianArgs=list(d=d.hessian)
	} else hessianArgs = NULL
	assign("hessianArgs", hessianArgs, env)

	#
	# Offset
	#

	offset.value = 0
	if(!missing(offset) && !class(offset) == "formula"){
		stop("Argument 'offset' must be a formula (e.g. ~ 1+x^2).")
	} else if(!missing(offset)){

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

	# Initial checks are done
	nonlinear.params <- names(start) #=> in the order the user wants

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

	if(missing(theta.init)){
		theta.init = NULL
	} else {
		if(!is.numeric(theta.init) || length(theta.init)!=1 || theta.init<=0) stop("the argument 'theta.init' must be a strictly positive scalar.")
	}

	model0 <- get_model_null(env, theta.init)
	theta.init = model0$theta

	# For the negative binomial:
	if(family=="negbin"){
		params = c(params, ".theta")
		start = c(start, theta.init)
		names(start) = params
		upper = c(upper, Inf)
		lower = c(lower, 1e-3)
		# DEPREC:
		# 		theta = 1 # no matter the init, its handled in the function getting theta
		# 		assign(".theta", theta, env)
	} else if(family=="tobit"){
		params = c(params, ".sigma")
		start = c(start, 1)
		names(start) = params
		upper = c(upper, Inf)
		lower = c(lower, 1e-3)
	}

	# On balance les donnees a utiliser dans un nouvel environnement
	for(i in varnames) assign(i, data[[i]], env)
	if(isLinear) assign("linear.mat", linear.mat, env)
	if(isGradient) assign(".call_gradient", nl.gradient[[2]], env)

	####
	#### Dummy Special cases ####
	####

	if( (family == "poisson" && Q==2) || (useCG && family == "gaussian" && Q>=2) ){
		# Creation of the matrices

		n = length(lhs)
		mat_all = list()

		for(q in 1:Q){
			mat = Matrix(0, nbCluster[q], n, sparse = TRUE)
			mat[cbind(dum_all[[q]], 1:n)] = 1
			mat_all[[q]] = mat
		}

		assign(".mat_all", mat_all, env)
	}

	# This is specific to the Gaussian case (for multi clusters)
	if(useCG && family == "gaussian" && Q>=2){

		n = length(lhs)

		# 2) Creation of the main (big) matrix A

		A1 = Matrix(0, sum(nbCluster), sum(nbCluster), sparse = TRUE)

		adj_nb =  c(0, cumsum(nbCluster))

		for(q1 in 1:(Q-1)){
			for(q2 in (q1+1):Q){

				mat1 = mat_all[[q1]]
				mat2 = mat_all[[q2]]

				nb = tcrossprod(mat1, mat2)
				qui = which(nb>0, arr.ind = TRUE)
				A1[cbind(adj_nb[q1] + qui[,1], adj_nb[q2] + qui[, 2])] = nb[qui]
			}
		}

		A2 = (A1 + t(A1))
		A = Diagonal(sum(nbCluster), x=unlist(obs_per_cluster_all)) + A2

		# Save
		assign(".A", A, env)
	}


	####
	#### Sending to the env ####
	####

	useExp = family %in% c("poisson", "logit", "negbin")

	# The dummies
	assign("isDummy", isDummy, env)
	if(isDummy){
		assign(".dummy", dum_all, env)
		assign(".dumMat_cpp", dumMat_cpp, env)
		assign(".nbCluster", nbCluster, env)
		assign(".sum_y", sum_y_all, env)
		assign(".tableCluster", obs_per_cluster_all, env)

		# the saved dummies
		if(useExp){
			assign(".savedDummy", rep(1, length(lhs)), env)
		} else {
			assign(".savedDummy", rep(0, length(lhs)), env)
		}


		assign(".orderCluster", NULL, env)
		if(family %in% c("negbin", "logit")){
			# we also add this variable used in cpp
			orderCluster_all = list()
			for(i in 1:Q){
				orderCluster_all[[i]] = order(dum_all[[i]]) - 1
			}
			assign(".orderCluster", orderCluster_all, env)
		}

		assign(".useCG", useCG, env)
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
	assign("debug", debug, env)
	assign("jacobian.method", jacobian.method, env)
	assign(".famFuns", famFuns, env)
	assign(".family", family, env)
	assign("iter", 0, env)
	# Pour gerer les valeurs de mu:
	assign(".coefMu", list(), env)
	assign(".valueMu", list(), env)
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

	#
	# if there is only the intercept and cluster => we estimate only the clusters
	#

	if(!isLinear && !isNonLinear && Q>0){
		if(family == "negbin"){
			stop("To estimate the negative binomial model, you need at least one variable. (The estimation of the model with only the clusters is to be implemented.)")
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
	mu = eval(nl.call, envir= env)

	# On sauvegarde les valeurs de la partie non lineaire
	assign(".nbMaxSave", NLsave, env) # nombre maximal de valeurs a sauvegarder
	assign(".nbSave", 1, env)  # nombre de valeurs totales sauvegardees
	assign(".savedCoef", list(start[nonlinear.params]), env)
	assign(".savedValue", list(mu), env)
	if(isLinear) mu <- mu + c(linear.mat%*%unlist(start[linear.params]))

	if(length(mu)!=nrow(data)) stop("Wow, must be a big problem... length(lhs)!=length(eval(fml))")
	if(anyNA(mu)) stop("Hum, must be a problem, evaluating formula returns NA.\nMaybe only these starting values can't be computed, or maybe there's another BIGGER problem.")

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

	#Mise en place du calcul du gradient
	gradient = NULL
	if(opt_method=="nlminb" | optim.method%in%c("BFGS", "CG", "L-BFGS-B")) gradient = ll_glm_gradient
	hessian <- NULL
	if(useHessian) hessian <- ll_glm_hessian

	# GIVE PARAMS
	if(!is.null(dots$give.params) && dots$give.params) return(list(coef=start, env=env))

	#
	# Maximizing the likelihood
	#

	opt <- NULL
	if(opt_method=="nlminb"){
# 		try(opt <- nlminb(start=start, objective=ll_glm, env=env, lower=lower, upper=upper, gradient=gradient, hessian=hessian), silent=FALSE)
		opt <- stats::nlminb(start=start, objective=ll_glm, env=env, lower=lower, upper=upper, gradient=gradient, hessian=hessian, control=opt.control)
	} else if(opt_method=="optim"){
		try(opt <- stats::optim(par=start, fn=ll_glm, gr=gradient, env=env, hessian=FALSE, control=opt.control), silent=FALSE)
		opt$objective <- opt$value
		opt$message = opt$convergence
	}

	if(is.null(opt)){
		stop("Could not achieve maximization.")
	}

	convStatus = TRUE
	if(opt_method=="nlminb" && !opt$message %in% c("relative convergence (4)", "both X-convergence and relative convergence (5)")){
		warning("The result is not reliable, the optimization did not converge.", call. = FALSE)
		convStatus = FALSE
	} else if(opt_method=="optim"){
		if(opt$convergence!=0){
			warning("The result is not reliable, the optimization did not converge.\nConvergence code", opt$convergence, ". See optim help page for more info.", call. = FALSE)
			opt$message = paste0("No convergence. Status ", opt$convergence)
			convStatus = FALSE
		} else opt$message = "Convergence"
	}

	####
	#### After Maximization ####
	####

	coef <- opt$par

	# The Hessian
	hessian = ll_glm_hessian(coef, env=env)
	# we add the names of the non linear variables in the hessian
	if(isNonLinear){
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
			hessian_noBounded = hessian[-which(isBounded), -which(isBounded)]

			boundText = ifelse(coef_NL == upper_bound, "Upper bounded", "Lower bounded")[isBounded]

			attr(isBounded, "type") = boundText
		}

	}

	# Variance

	var <- NULL
	try(var <- solve(hessian_noBounded), silent = TRUE)
	if(is.null(var)){
		#var <- MASS::ginv(hessian_noBounded)
		if(!isDummy){
			Qr = qr(hessian_noBounded)
			collvar = params[Qr$pivot[-(1:Qr$rank)]]
			var <- MASS::ginv(hessian_noBounded)
			warning("Could not achieve to get the covariance matrix. The information matrix was singular.\nTry to manually suppress some collinear variables. FYI the suspects are:\n", paste(collvar, collapse=", "), call. = FALSE)
		} else {
			var <- MASS::ginv(hessian_noBounded)
			collvar = params[diag(var)==0]
			if(length(collvar)>0){
				warning("Could not achieve to get the covariance matrix. The information matrix was singular.\nTry to manually suppress some collinear variables. FYI the suspects (collinear with the clusters) are:\n", paste(collvar, collapse=", "), call. = FALSE)
			} else {
				warning("Could not achieve to get the covariance matrix. The information matrix was singular.\nTry to manually suppress some collinear variables. They may be collinear with the clusters.", call. = FALSE)
			}
		}

		# WARNING: we DO NOT give the covariance matrix => too dangerous ! and it is not interpretable anyway
		# => we add a message in "summary" and an option to compute the generalized inverse of the hessian

		var = var*NA
		se = diag(var)
	} else {
		se = diag(var)
		se[se<0] = NA
		if(anyNA(se)) warning("CAUTION: Variance needs to be 'eigenfixed'.", call. = FALSE)
		se = sqrt(se)
	}

	# To handle the bounded coefficient, we set its SE to NA
	if(any(isBounded)){
		se = se[params]
		names(se) = params
	}

	zvalue <- coef/se
	pvalue <- 2*pnorm(-abs(zvalue))

	# We add the information on the bound for the se & update the var to drop the bounded vars
	if(any(isBounded)){
		se[!isBounded] = decimalFormat(se[!isBounded])
		se[isBounded] = boundText
	}

	# coeftable <- cbind(coef, se, zvalue, pvalue)
	# colnames(coeftable) <- c("Estimate", "Std. Error", "z value",  "Pr(>|z|)")
	# rownames(coeftable) <- params
	# class(coeftable) <- "coeftest"

	coeftable <- data.frame("Estimate"=coef, "Std. Error"=se, "z value"=zvalue, "Pr(>|z|)"=pvalue, stringsAsFactors = FALSE)
	names(coeftable) <- c("Estimate", "Std. Error", "z value",  "Pr(>|z|)")
	row.names(coeftable) <- params

	attr(se, "type") = attr(coeftable, "type") = "Standard"

	mu = get_mu(coef, env)

	# calcul pseudo r2
	loglik <- -opt$objective # moins car la fonction minimise
	ll_null <- model0$loglik
	# degres de liberte
	df_k = length(coef)
	if(isDummy) df_k = df_k + sum(sapply(dum_all, max) - 1) + 1
	# dummies are constrained, they don't have full dof (cause you need to take one value off for unicity)
	pseudo_r2 <- 1 - (loglik-df_k)/(ll_null-1)
	# pseudo_r2 <- 1 - loglik/ll_null # NON Adjusted => NO

	# deprecated resids
	# null.resids <- lhs-model0$constant
	# new null resids
	# null.resids <- lhs-mean(lhs)

	# Calcul residus
	expected.predictor = famFuns$expected.predictor(mu, env)
	resids = lhs - expected.predictor

	# calcul squared corr
	if(sd(expected.predictor) == 0){
		sq.cor = NA
	} else {
		sq.cor = stats::cor(lhs, expected.predictor)**2
	}

	# calcul r2 naif
	naive.r2 = 1 - sum(resids**2)/sum((lhs-mean(lhs))**2)

	# The scores
	scores = ll_glm_scores(coef, env)
	if(isNonLinear){
		# we add the names of the non linear params in the score
		colnames(scores) = params
	}


	res <- list(coef=coef, coeftable=coeftable, loglik=loglik, iterations=opt$iterations, n=length(lhs), k=df_k, call=call, fml=fml, NL.fml=NL.fml, ll_null=ll_null, pseudo_r2=pseudo_r2, naive.r2=naive.r2, message=opt$message, convStatus=convStatus, sq.cor=sq.cor, expected.predictor=expected.predictor, hessian=hessian, cov.unscaled=var, bounds=bounds, isBounded=isBounded, se=se, scores=scores, family=family, resids=resids)

	# Dummies
	if(isDummy){
		# NOT YET IMPLEMENTED FOR VARIOUS CLUSTERS
		dummies = attr(mu, "mu_dummies")
		if(useExp){
			dummies = log(dummies)
		}

		if(FALSE & all(dum_names == "(Intercept)")){
			# get the variance of the intercept
			# and include it in the coefficients
			ptm = proc.time()
			res$Intercept = dummies
			jacob.mat = get_Jacobian(coef, env)
			se = famFuns$dumVar(jacob.mat, mu, env, coef)
			zvalue = dummies/se
			pvalue = 2*pnorm(-abs(zvalue))
			line_intercept = matrix(c(dummies, se, zvalue, pvalue), 1, 4)
			rownames(line_intercept) = "(Intercept)"
			coeftable = rbind(line_intercept, res$coeftable)
			print(coeftable)
			cat("SD of the intercept in:", (proc.time()-ptm)[3], "\n")
		} else {
			res$dummies = dummies
			res$clusterNames = cluster

			id_dummies = list()
			for(i in 1:length(cluster)){
				# id_dummies[[cluster[i]]] = factor(dum_all[[i]], labels=dum_names[[i]])
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
	}

	if(family == "negbin"){
		theta = coef[".theta"]
		res$theta = theta
	}

	class(res) <- "femlm"

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

	# mu, using the offset
	mu_noDum = offset.value
	if(length(mu_noDum) == 1) mu_noDum = rep(mu_noDum, nobs)

	# we create the exp of mu => used for later functions
	exp_mu_noDum = NULL
	if(useExp){
		exp_mu_noDum = exp(mu_noDum)
	}

	dummies = getDummies(mu_noDum, exp_mu_noDum, env, coef)

	if(useExp){
		# despite being called mu, it is in fact exp(mu)!!!
		mu = exp_mu_noDum*dummies
	} else {
		mu = mu_noDum + dummies
	}

	#
	# 2nd step => saving information
	#

	dum_all = get(".dummy", env)
	famFuns = get(".famFuns", env)
	lhs = get(".lhs", env)

	# The log likelihoods
	loglik = famFuns$ll(lhs, mu, env, coef)
	ll_null = model0$loglik

	# degres de liberte
	df_k = sum(sapply(dum_all, max) - 1) + 1
	pseudo_r2 = 1 - loglik/ll_null # NON Adjusted

	# Calcul residus
	expected.predictor = famFuns$expected.predictor(mu, env)
	resids = lhs - expected.predictor

	# calcul squared corr
	if(sd(expected.predictor) == 0){
		sq.cor = NA
	} else {
		sq.cor = stats::cor(lhs, expected.predictor)**2
	}

	# calcul r2 naif
	naive.r2 = 1 - sum(resids**2) / sum((lhs - mean(lhs))**2)

	res = list(loglik=loglik, n=length(lhs), k=df_k, call=call, ll_null=ll_null, pseudo_r2=pseudo_r2, naive.r2=naive.r2, sq.cor=sq.cor, expected.predictor=expected.predictor, family=family)
	#
	# Information on the dummies

	if(useExp){
		dummies = log(dummies)
	}

	res$dummies = dummies
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

ll_glm_hessian <- function(coef, env){
	# Computes the hessian
	# cat("in Hessian:", as.vector(coef), "\n")
	debug = get("debug", env)
	if(debug) ptm = proc.time()
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
	mu = get_savedMu(coef, env)

	jacob.mat = get_Jacobian(coef, env)

	ll_d2 = famFuns$ll_d2(y, mu, coef)
	if(isDummy){
		dxi_dbeta = deriv_xi(jacob.mat, ll_d2, env, coef)
		jacob.mat = jacob.mat + dxi_dbeta
	} else dxi_dbeta = 0

	hessVar = crossprod(jacob.mat, jacob.mat * ll_d2)

	if(isNL){
		#we get the 2nd derivatives
		z = numDeriv::genD(evalNLpart, coef[nonlinear.params], env=env, method.args = hessianArgs)$D[, -(1:k), drop=FALSE]
		ll_dl = famFuns$ll_dl(y, mu, coef=coef, env=env)
		id_r = rep(1:k, 1:k)
		id_c = c(sapply(1:k, function(x) 1:x), recursive=TRUE)
		H = matrix(0, nrow=k, ncol=k)
		H[cbind(id_r, id_c)] = H[cbind(id_r, id_c)] = colSums(z*ll_dl)
	} else H = 0

	# on ajoute la partie manquante
	if(isNL) hessVar[1:k, 1:k] = hessVar[1:k, 1:k] + H

	if(family=="negbin"){
		theta = coef[".theta"]
		ll_dx_dother = famFuns$ll_dx_dother(y, mu, coef, env)

		if(isDummy){
			dxi_dother = deriv_xi_other(ll_dx_dother, ll_d2, env, coef)
		} else {
			dxi_dother = 0
		}

		# calcul des derivees secondes vav de theta
		h.theta.L = famFuns$hess.thetaL(theta, jacob.mat, y, dxi_dbeta, dxi_dother, ll_d2, ll_dx_dother)
		hessVar = cbind(hessVar, h.theta.L)
		h.theta = famFuns$hess_theta_part(theta, y, mu, dxi_dother, ll_dx_dother, ll_d2)
		hessVar = rbind(hessVar, c(h.theta.L, h.theta))

# 		theta = attr(mu, ".theta")
# 		print(theta)
		#Jacob.mat est la derivee totale de mu vav de beta (ie avec les dummies)
# 		hessVar = hessVar + famFuns$hess_theta_part(theta, jacob.mat, y, mu)
	} else if(family=="tobit"){
		sigma = coef[".sigma"]
		h.sigma.L = famFuns$hess.sigmaL(sigma, jacob.mat, y, mu, dxi_dbeta, ll_d2)
		hessVar = cbind(hessVar, h.sigma.L)
		h.sigma = famFuns$hess.sigma(sigma, y, mu)
		hessVar = rbind(hessVar, c(h.sigma.L, h.sigma))
		print(hessVar)
	}

	if(debug) cat("\nHessian: ", (proc.time()-ptm)[3], "s")
# print(hessVar) ; print(class(hessVar))
	- hessVar
}

ll_glm_gradient <- function(coef, env){
	# cat("gradient:\n") ; print(as.vector(coef))

	params = get("params", env)
	names(coef) = params
	nonlinear.params = get("nonlinear.params", env)
	linear.params = get("linear.params", env)
	famFuns = get(".famFuns", env)
	family = get(".family", env)
	y = get(".lhs", env)
	mu = get_savedMu(coef, env)

	# calcul de la jacobienne
	res <- list() #stocks the results

	# cat("\tgetting jacobian")
	# ptm = proc.time()
	jacob.mat = get_Jacobian(coef, env)
	# cat("in", (proc.time()-ptm)[3], "s.\n")

	# cat("\tComputing gradient ")
	# ptm = proc.time()
	# res = famFuns$grad(jacob.mat, y, mu, env, coef)
	res = getGradient(jacob.mat, y, mu, env, coef)
	# cat("in", (proc.time()-ptm)[3], "s.\n")
	names(res) = c(nonlinear.params, linear.params)

	if(family=="negbin"){
		theta = coef[".theta"]
		res[".theta"] = famFuns$grad.theta(theta, y, mu)
	}
# print(res)
	return(-unlist(res[params]))
}

ll_glm_scores <- function(coef, env){
	# Computes the scores (Jacobian)
	params = get("params", env)
	names(coef) <- params
	famFuns = get(".famFuns", env)
	family = get(".family", env)
	y = get(".lhs", env)
	mu = get_savedMu(coef, env)

	jacob.mat = get_Jacobian(coef, env)
	scores = getScores(jacob.mat, y, mu, env, coef)

	if(family=="negbin"){
		theta = coef[".theta"]
		score.theta = famFuns$scores.theta(theta, y, mu)
		scores = cbind(scores, score.theta)
		# DEPREC (theta-conditionned)
# 		isDummy = get("isDummy", env)
# 		if(isDummy){
# 			ll_d2 = famFuns$ll_d2(y, mu, coef)
# 			dxi_dbeta = deriv_xi(jacob.mat, ll_d2, env, coef)
# 			jacob.mat = jacob.mat + dxi_dbeta
# 		}
# 		theta = attr(mu, ".theta")
# 		scores = scores + famFuns$scores_theta_part(theta, jacob.mat, y, mu)
	}

	return(scores)
}

ll_glm <- function(coef, env){
	# Log likelihood
	# cat("LL:\n") ; print(coef)
	# misc funs
	iter = get("iter", env) + 1
	assign("iter", iter, env)
	debug = get("debug", env)
	ptm = proc.time()
	if(debug) cat("\nIter", iter, "- Evaluation LL:", sprintf("%1.1e", as.vector(coef)), "\n")

	# computing the LL
	famFuns = get(".famFuns", env)
	family = get(".family", env)
	y <- get(".lhs", env)

	if(any(is.na(coef))) stop("Divergence... (some coefs are NA)\nTry option debug=TRUE to see the problem.")

	mu = get_mu(coef, env)

	# for the NEGBIN, we add the coef
	ll = famFuns$ll(y, mu, env, coef)

	if(debug) cat("LL =", ll, " (", (proc.time()-ptm)[3], " s)", sep = "")
	if(ll==(-Inf)) return(1e308)
	return(-ll) # je retourne -ll car la fonction d'optimisation minimise
}

ll_glm_TEST_score <- function(coef, env){
	# Used to compute the scores numerically
	# Not user oriented

	debug <- get("debug", env)
	if(debug) print(coef)

	# computing the LL
	famFuns = get(".famFuns", env)
	y <- get(".lhs", env)

	if(any(is.na(coef))) stop("Divergence... (some coefs are NA)\nTry option debug=TRUE to see the problem.")

	mu = get_mu(coef, env)

	ll = famFuns$ll_TEST_score(y, mu, env, coef)

	return(ll)
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
	if(!is.null(names(coef))) coef = coef[nonlinear.params]
	else if (length(coef)!=length(nonlinear.params)) stop("Problem with the length of the NL coefficients.")

	if(nbMaxSave == 0){
		for(var in nonlinear.params) assign(var, coef[var], env)
		y_nl <- eval(nl.call, envir= env)
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

get_mu = function(coef, env){
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

	# indicator of whether we compute the exp(mu)
	# if usExp==TRUE, alors mu_dummies est aussi l'exponentielle des dummies
	# BEWARE, the variable mu will be called mu although it can be equal to exp(mu)!!!!!!!
	# this is a notational abuse in order to simplify the coding
	useExp = family %in% c("poisson", "logit", "negbin")

	# For managing mu:
	coefMu = get(".coefMu", env)
	valueMu = get(".valueMu", env)
	wasUsed = get(".wasUsed", env)
	if(wasUsed){
		coefMu = valueMu = list()
		assign(".wasUsed", FALSE, env)
	}

	if(length(coefMu)>0) for(i in 1:length(coefMu)) if(all(coef==coefMu[[i]])) return(valueMu[[i]])

	if(isNL){
		muNL = evalNLpart(coef, env)
	} else muNL = 0

	if(isLinear){
		linear.params = get("linear.params", env)
		linear.mat = get("linear.mat", env)
		mu_L = c(linear.mat%*%coef[linear.params])
	} else mu_L = 0

	mu_noDum = muNL + mu_L + offset.value

	# we create the exp of mu => used for later functions
	exp_mu_noDum = NULL
	if(useExp){
		exp_mu_noDum = exp(mu_noDum)
	}

	if(isDummy){
		# we get back the last dummy
		mu_dummies = getDummies(mu_noDum, exp_mu_noDum, env, coef)
	} else {
		if(useExp){
			mu_dummies = 1
		} else {
			mu_dummies = 0
		}
	}

	# We add the value of the dummy to mu
	if(useExp){
		# despite being called mu, it is in fact exp(mu)!!!
		mu = exp_mu_noDum*mu_dummies
	} else {
		mu = mu_noDum + mu_dummies
	}

	if(isDummy){
		# BEWARE, if useExp, it is equal to exp(dummies)
		attr(mu, "mu_dummies") = mu_dummies
	}

	if(length(mu)==0) mu = rep(mu, nobs)

	# we save the value of mu:
	coefMu = append(coefMu, list(coef))
	valueMu = append(valueMu, list(mu))
	assign(".coefMu", coefMu, env)
	assign(".valueMu", valueMu, env)

	return(mu)
}

get_savedMu = function(coef, env){
	# This function gets the mu without computation
	# It follows a LL evaluation
	coefMu = get(".coefMu", env)
	valueMu = get(".valueMu", env)
	assign(".wasUsed", TRUE, env)

	if(length(coefMu)>0) for(i in 1:length(coefMu)) if(all(coef==coefMu[[i]])){
		# cat("coef nb:", i, "\n")
		return(valueMu[[i]])
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
			jacob.mat = cbind(jacob.mat, linear.mat)
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

ll0_nlglm <- function(lhs, env){
	#I have the closed form of the ll0
	famFuns = get(".famFuns", env)
	family = get(".family", env)
	y = get(".lhs", env)

	if(family=="negbin"){
		start = c(0, 1)
		lower = c(-Inf, 1e-4)
	} else {
		start = 0
		lower = NULL
	}

	opt <- nlminb(start=start, objective=famFuns$ll0, y=y, gradient=famFuns$grad0, lower=lower)

	return(list(loglik=-opt$objective, constant=opt$par[1]))
}

get_model_null <- function(env, theta.init){
	# I have the closed form of the ll0
	famFuns = get(".famFuns", env)
	family = get(".family", env)
	N = get("nobs", env)
	y = get(".lhs", env)

	# one of the elements to be returned
	theta = NULL

	if(family == "poisson"){
		# There is a closed form

		if(".lfactorial" %in% names(env)){
			lfact = get(".lfactorial", env)
		} else {
			lfact = sum(lfactorial(y))
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
			lgamm = sum(lgamma(y + 1))
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

		opt <- nlminb(start=theta.guess, objective=famFuns$ll0_theta, y=y, gradient=famFuns$grad0_theta, lower=1e-3, mean_y=mean_y, invariant=invariant, hessian = famFuns$hess0_theta)

		loglik = -opt$objective
		theta = opt$par
	}

	return(list(loglik=loglik, constant=constant, theta = theta))
}

getGradient = function(jacob.mat, y, mu, env, coef, ...){
	famFuns = get(".famFuns", env)
	ll_dl = famFuns$ll_dl(y, mu, coef=coef, env=env)
	c(crossprod(jacob.mat, ll_dl))
}

getScores = function(jacob.mat, y, mu, env, coef, ...){
	famFuns = get(".famFuns", env)
	isDummy = get("isDummy", env)

	ll_dl = famFuns$ll_dl(y, mu, coef=coef, env=env)
	scores = jacob.mat* ll_dl

	if(isDummy){
		ll_d2 = famFuns$ll_d2(y, mu, coef=coef, env=env)
		dxi_dbeta = deriv_xi(jacob.mat, ll_d2, env, coef)
		scores = scores + dxi_dbeta * ll_dl
	}

	return(as.matrix(scores))
}

getDummies = function(mu, exp_mu, env, coef){
	# function built to get all the dummy variables
	# We retrieve past dummies (that are likely to be good
	# starting values)
	mu_dummies = get(".savedDummy", env)
	family = get(".family", env)
	eps.cluster = get(".eps.cluster", env)
	debug = get("debug", env)
	if(debug) ptm = proc.time()

	iterMax = 10000
	dum_all = get(".dummy", env)
	sum_y_all = get(".sum_y", env)
	tableCluster_all = get(".tableCluster", env)
	orderCluster_all = get(".orderCluster", env)
	nbCluster = get(".nbCluster", env)
	Q = length(dum_all)

	# whether we use the eponentiation of mu
	useExp = family %in% c("poisson", "logit", "negbin")

	if(useExp){
		exp_mu_in = exp_mu * mu_dummies
	} else {
		mu_in = mu + mu_dummies
	}
# browser()
	useCG = get(".useCG", env)

	if(useCG && family == "gaussian" && Q>=2){
		# in this case, we can use the algorithm of the conjugate gradient
		# This algorithm is efficient only for complicated cases
		# otherwise, there is no much value added

		A = get(".A", env)
		# we get the value of the constant: we need to compute the sum of the x's for each cluster
		sum_x_all = list()
		for(q in 1:Q){
			dum = dum_all[[q]]
			k  = nbCluster[q]
			sum_x_all[[q]] = cpp_tapply_vsum(k, mu_in, dum)
		}

		b = (unlist(sum_y_all) - unlist(sum_x_all))

		all_clusters_coef = simple_CG(A, b, eps.cluster)

		# update of the dummies
		adj_nb =  c(0, cumsum(nbCluster))
		for(q in 1:Q){
			dum = dum_all[[q]]
			mu_dummies = mu_dummies + all_clusters_coef[(adj_nb[q]+1):adj_nb[q+1]][dum]
		}

	} else if(family == "poisson" && Q==2){
		# This method is usually faster, however, in simple cases it might be not efficient because of the
		# setup of the matrices
		# for complicated cases, it is much more efficient

		# Specific version that loops only over the dummy with the lowest number of cases

		# We select the cluster with the lowest nber of cases
		if(nbCluster[1] < nbCluster[2]){
			i = 1
			j = 2
		} else {
			i = 2
			j = 1
		}

		# Creation of the necessary preliminary data

		mat_all = get(".mat_all", env)
		M1 = mat_all[[i]]
		M2 = mat_all[[j]]

		At = M1%*%(t(M2)*exp_mu_in)
		# Bt = M2%*%(t(M1)*exp_mu_in)
		Bt = t(At)

		# the main loop
		alpha = rep(1, nbCluster[i])
		ca = as.numeric(sum_y_all[[i]])
		cb = as.numeric(sum_y_all[[j]])
		for(iter in 1:iterMax){
			alpha_old = alpha
			alpha = ca / (At%*%(cb / (Bt%*%alpha)))

			diff = max(abs(log(range(alpha_old / alpha))))
			if(diff < eps.cluster) break
		}

		# DEPREC: now we use directly the exponentiation
		# lbeta = log(cb) - log(Bt%*%alpha)
		# mu_dummies = mu_dummies + log(alpha)[dum_all[[i]]] + lbeta[dum_all[[j]]]

		beta = cb / Bt%*%alpha
		mu_dummies = mu_dummies * alpha[dum_all[[i]]] * beta[dum_all[[j]]]

	} else if (family == "poisson"){
		# Quicker than before, 100% sure

		# We loop directly without using logs

		for(iter in 1:iterMax){

			exp_mu0 = exp_mu_in

			for(q in 1:Q){
				dum = dum_all[[q]]

				# get the dummies
				exp_mu_dum = sum_y_all[[q]] / cpp_tapply_vsum(nbCluster[q], exp_mu_in, dum)
				# add them to the stack
				exp_mu_in = exp_mu_in*exp_mu_dum[dum]
			}

			if(Q == 1) break

			diff <- max(abs(log(range(exp_mu0/exp_mu_in)))) # max error
			if(diff < eps.cluster) break
		}

		# print(iter)
		if(iter==iterMax) warning("[getting dummies] iteration limit reached (max diff=", signif(diff), ").", call. = FALSE, immediate. = TRUE)

		# mu_dummies = log(exp_mu_in) - mu
		mu_dummies = exp_mu_in / exp_mu

	} else if (TRUE && useCG && family %in% c("negbin","logit")){
		# The new versions of the logit and negbin
		# no much more efficient because the NR is inefficient when we try to find the exp() of the dummy
		# => so we compute the NR on the log of the dummy and then exponentiate it, not very efficient...

		for(iter in 1:iterMax){

			exp_mu0 = exp_mu_in

			for(q in 1:Q){
				dum = dum_all[[q]]
				sum_y = sum_y_all[[q]]
				tableCluster = tableCluster_all[[q]]
				orderCluster = orderCluster_all[[q]]
				# get the dummies
				mu_dum = dichoNR_Cpp_exp(dum, exp_mu_in, env, coef, sum_y, orderCluster, tableCluster)
				# add them to the stack
				exp_mu_in = exp_mu_in * exp(mu_dum)[dum]
			}

			if(Q==1) break

			diff <- max(abs(log(range(exp_mu0/exp_mu_in)))) # max error
			# print(diff)

			if(anyNA(exp_mu_in)) stop("NA in cluster coefficients, unknown origin.")

			if(diff < eps.cluster) break

		}

		# print(iter)
		if(iter==iterMax) warning("[getting dummies] iteration limit reached (max diff=", signif(diff), ").", call. = FALSE, immediate. = TRUE)

		# mu_dummies = log(exp_mu_in) - mu
		mu_dummies = exp_mu_in / exp_mu

	} else {
		# in the case we are in the useExp case, we need to make some changes for consistency
		if(useExp) mu_in = log(exp_mu_in)

		for(iter in 1:iterMax){

			mu0 = mu_in

			for(q in 1:Q){
				dum = dum_all[[q]]
				sum_y = sum_y_all[[q]]
				tableCluster = tableCluster_all[[q]]
				orderCluster = orderCluster_all[[q]]
				# get the dummies
				mu_dum = computeDummies(dum, mu_in, env, coef, sum_y, orderCluster, tableCluster)
				# add them to the stack
				mu_in = mu_in + mu_dum[dum]
			}

			if(anyNA(mu_in)) stop("Dummies could not be computed.")

			if(Q==1) break

			diff <- max(abs(mu0-mu_in)) # max error
			# cat("iter =", iter, " Diff =", diff, " sum =", sum(abs(mu0-mu_in)), "\n")
			if(diff < eps.cluster) break
		}

		# print(iter)
		if(iter==iterMax) warning("[getting dummies] iteration limit reached (max diff=", signif(diff), ").", call. = FALSE, immediate. = TRUE)

		mu_dummies = mu_in - mu

		# Idem, if useExp, we need to make some changes
		if(useExp) mu_dummies = exp(mu_dummies)

	}

	# we save the dummy:
	assign(".savedDummy", mu_dummies, env)

	if(debug) cat("\nDummies: ", (proc.time()-ptm)[3], "s\t")

	mu_dummies
}

computeDummies = function(dum, mu, env, coef, sum_y, orderCluster, tableCluster){
	family = get(".family", env)
	famFuns = get(".famFuns", env)
	y = get(".lhs", env)
	eps.NR = get(".eps.NR", env)

	# if family is poisson or gaussian, there is a closed form
	if(family%in%c("poisson", "gaussian")){
		return(famFuns$closedFormDummies(dum, y, mu, env, sum_y, tableCluster))
	}

	#
	# For non closed-form dummies:
	#

	x1 = dichoNR_Cpp(dum, mu, env, coef, sum_y, orderCluster, tableCluster)

	x1
}

get_mu_noDum = function(coef, env){
	# TO REMOVE
	# This function computes the RHS of the equation
	# mu_L => to save one matrix multiplication
	isNL <- get("isNL", env)
	isLinear <- get("isLinear", env)
	params <- get("params", env)
	names(coef) <- params

	if(isNL){
		muNL = evalNLpart(coef, env)
	} else muNL = 0

	if(isLinear){
		linear.params <- get("linear.params", env)
		linear.mat <- get("linear.mat", env)
		mu_L <- c(linear.mat%*%coef[linear.params])
	} else mu_L = 0

	mu_noDum = muNL + mu_L

	return(mu_noDum)
}

dichoNR_Cpp_exp = function(dum, exp_mu, env, coef, sum_y, orderCluster, tableCluster){
	# function that combines dichotomy
	# We are sure that there is a zero, so dichotomy will work 100%
	# browser()
	nb_cases = length(tableCluster)
	y = get(".lhs", env)
	N = length(y)

	# 1st we get the boundaries
	famFuns = get(".famFuns", env)
	eps.NR = get(".eps.NR", env)

	# Get the init, the inf and sup boundaries
	minMax = cpp_conditional_minMax(nb_cases, N, exp_mu, dum, tableCluster)

	# Init => 0, car a la fin, cela doit converger vers 0
	init = rep(0, nb_cases)

	# mu est la valeur estimee en les parametres (hors les dummies)
	# plus les mu sont eleves, plus la dummy doit etre petite pour compenser
	borne_inf = famFuns$guessExpDummy(sum_y, tableCluster, minMax[, 2]) # Value > 0, (on fait moins muMax)
	borne_sup = famFuns$guessExpDummy(sum_y, tableCluster, minMax[, 1]) # Value < 0

	family = get(".family", env)
	family_nb = switch(family, poisson=1, negbin=2, logit=3, gaussian=4, probit=5)

	theta = coef[".theta"]
	# ptm <- proc.time()
	res = new_RcppDichotomyNR(N = N, K = nb_cases, family = family_nb, theta,
								 epsDicho=eps.NR, lhs=y, exp_mu, log(borne_inf),
								 log(borne_sup), orderCluster, tableCluster)

	# print(proc.time() - ptm)
	res
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

	# Init => 0, car a la fin, cela doit converger vers 0
	init = rep(0, nb_cases)

	# mu est la valeur estimee en les parametres (hors les dummies)
	# plus les mu sont eleves, plus la dummy doit etre petite pour compenser
	borne_inf = famFuns$guessDummy(sum_y, tableCluster, minMax[, 2]) # Value > 0, (on fait moins muMax)
	borne_sup = famFuns$guessDummy(sum_y, tableCluster, minMax[, 1]) # Value < 0

	family = get(".family", env)
	family_nb = switch(family, poisson=1, negbin=2, logit=3, gaussian=4, probit=5)

	theta = coef[".theta"]
	# ptm <- proc.time()
	res = RcppDichotomyNR(N = N, K = nb_cases, family = family_nb, theta,
							  epsDicho=eps.NR, lhs=y, mu, borne_inf,
							  borne_sup, orderCluster, tableCluster)

	# print(proc.time() - ptm)
	res
}

deriv_xi = function(jacob.mat, ll_d2, env, coef){
	#Derivee des dummies
	dum_all = get(".dummy", env)
	dumMat_cpp = get(".dumMat_cpp", env)
	nbCluster = get(".nbCluster", env)
	family = get(".family", env)
	eps.deriv = get(".eps.deriv", env)
	Q = length(dum_all)

	# Controls (for clusters >= 2)
	nobs = get("nobs", env)
	iterMax = 100

	if(Q == 1){
		dum = dum_all[[1]]
		k = max(dum)
		S_Jmu = cpp_tapply_sum(k, jacob.mat*ll_d2, dum)
		S_mu = cpp_tapply_vsum(k, ll_d2, dum)
		dxi_dbeta = - S_Jmu[dum, ] / S_mu[dum]
	} else {
		# The cpp way:

		N = length(ll_d2)
		K = length(coef) - (family=="negbin") # to avoid theta

		# We set the initial values for the first run
		if(!".sum_deriv" %in% names(env)){
			# init of the sum of the derivatives => 0
			init = matrix(0, N, K)
		} else {
			init = get(".sum_deriv", env)
		}

		if(family == "gaussian"){
			dxi_dbeta <- RcppPartialDerivative_gaussian(Q, N, K, epsDeriv = eps.deriv, jacob.mat, init, dumMat_cpp, nbCluster)
		} else {
			dxi_dbeta <- RcppPartialDerivative(Q, N, K, epsDeriv = eps.deriv, ll_d2, jacob.mat, init, dumMat_cpp, nbCluster)
		}

		# we save the values
		assign(".sum_deriv", dxi_dbeta, env)

	}

	as.matrix(dxi_dbeta)
}

deriv_xi_other = function(ll_dx_dother, ll_d2, env, coef){
	# derivative of the dummies wrt an other parameter
	dumMat_cpp = get(".dumMat_cpp", env)
	nbCluster = get(".nbCluster", env)
	dum_all = get(".dummy", env)
	eps.deriv = get(".eps.deriv", env)
	Q = length(dum_all)

	if(Q==1){
		dum = dum_all[[1]]
		k = max(dum)
		S_Jmu = cpp_tapply_vsum(k, ll_dx_dother, dum)
		S_mu = cpp_tapply_vsum(k, ll_d2, dum)
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

		dxi_dother <- RcppPartialDerivative_other(Q, N, epsDeriv = eps.deriv, ll_d2, ll_dx_dother, init, dumMat_cpp, nbCluster)

		# we save the values
		assign(".sum_deriv_other", dxi_dother, env)

	}

	as.matrix(dxi_dother)
}

simple_CG = function(A, b, eps.CG = 1e-6){
	# Simple algorithm to solve a linear system of equation (conjugate gradient)
	# The matrix A is always symmetric & positive definite by definition

	iterMax = 2000

	# Initialisations
	x_old = rep(0, ncol(A))
	r_old = as.numeric(b - A%*%x_old)
	p_old = r_old

	cp_r_old = crossprod(r_old)

	for(i in 1:iterMax){
		# some values
		g = as.numeric(A%*%p_old)

		alpha = as.numeric(cp_r_old / t(p_old)%*%g)
		x = x_old + alpha*p_old
		r = r_old - alpha*g

		cp_r_new = crossprod(r)
		beta = cp_r_new / cp_r_old
		p = r + beta*p_old

		# exit condition
		if(max(abs(r)) < eps.CG) break
		# print(max(abs(r)))

		# MAJ
		r_old = r
		p_old = p
		x_old = x

		cp_r_old = cp_r_new
	}

	if(i == iterMax){
		warning("[Conjugate gradient] Max. iterations reached (", iterMax, ") maxDiff = ", max(abs(r)))
	}

	x
}


gt = function(ptm) (proc.time() - ptm)[[3]]


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




