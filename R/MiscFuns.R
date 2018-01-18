

# Checks: getFE (verifier que tout marche bien) // summary => do.unclass (idem)

# URGENT: add possibility to compute or not the SE for lower/upper bounded coef

# ADD other diagnostics (like odd-ratio,  APE,  PAE for probit/logit + other R2)
# ADD model0 with dummies

#' A print facility for \code{femlm} objects. It can compute different types of standard errors.
#'
#' This function is very similar to usual \code{summary} functions as it provides the table of coefficients along with other information on the fit of the estimation.
#'
#' @method print femlm
#'
#' @param x A femlm object. Obtained using \code{\link[FENmlm]{femlm}}.
#' @param ... Other arguments to be passed to \code{\link[FENmlm]{summary.femlm}}.
#'
#' @seealso
#' See also the main estimation functions \code{\link[FENmlm]{femlm}} and \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # The data
#' n = 100
#' x = rnorm(n, 1, 5)**2
#' y = rnorm(n, -1, 5)**2
#' z = rpois(n, x*y) + rpois(n, 10)
#' base = data.frame(x, y, z)
#'
#' # Results of the Poisson..
#' est_poisson = femlm(z~log(x)+log(y), base, family="poisson")
#' # .. and of the Negative Binomial
#' est_negbin = femlm(z~log(x)+log(y), base, family="negbin")
#'
#' # Displaying the results
#' print(est_poisson)
#' print(est_negbin)
#'
#' # Changing the way the standard errors are computed:
#' summary(est_poisson, se="white") # similar to print(est_poisson, se="white")
#' summary(est_negbin, se="white")
#'
#' #
#' # Now with fixed-effects
#' #
#'
#' # Bilateral network
#' nb = 20
#' n = nb**2
#' k = nb
#' id1 = factor(rep(1:k, each=n/k))
#' id2 = factor(rep(1:(n/k), times=k))
#' x = rnorm(n, 1, 5)**2
#' y = rnorm(n, -1, 5)**2
#' z = rpois(n, x*y+rnorm(n, sd = 3)**2)
#' base = data.frame(x, y, z, id1, id2)
#'
#' # We want to use the ID's of each observation as a variable: we use the option cluster
#' est_poisson = femlm(z~log(x)+log(y), base, family="poisson", cluster=c("id1","id2"))
#' # Displaying the results
#' est_poisson
#' # now with twoway clustered santard-errors
#' summary(est_poisson, "twoway")
#'
print.femlm <- print.feNmlm <- function(x, ...){
	#Ajouter le sandwich estimator des variances dans le print

	x = summary(x, fromPrint = TRUE, ...)

	# if(!is.null(x$clusterNames)) cat("Standard errors clustered by \"", x$clusterNames[1], "\"\n", sep="")

	coeftable = x$coeftable
	# The type of SE
	se.type = attr(coeftable, "type")
	family_format = c(poisson="Poisson", negbin="Negative Binomial", logit="Logit", gaussian="Gaussian")

	msg = ifelse(is.null(x$call$NL.fml), "", "Non-linear ")
	cat(msg, "ML esimation, family = ", family_format[x$family], "\n", sep="")
	cat("Observations:", addCommas(x$n), "\n")
	if(!is.null(x$clusterSize)) cat("Cluster sizes: ", paste0(x$clusterNames, ": ", x$clusterSize, collapse=",  "), "\n", sep="")

	if(is.null(x$onlyCluster)){

		cat("Standard-errors type:", se.type, "\n")

		# The matrix of coefficients
		if(x$family=="negbin"){
			if(nrow(coeftable) == 2){
				# new_table = matrix(coeftable[1, ], 1, 4)
				# colnames(new_table) = colnames(coeftable)
				# rownames(new_table) = rownames(coeftable)[1]
				# class(new_table) = "coeftest"
				new_table = coeftable[1, , FALSE]
			} else {
				new_table = coeftable[-nrow(coeftable), ]
			}

			myPrintCoefTable(new_table)
			# stats::printCoefmat(new_table)
			cat("\nDispersion parameter: theta =", coeftable[".theta", 1], "\n")
		} else if(x$family=="tobit"){
			if(nrow(coeftable) == 2){
				# new_table = matrix(coeftable[1, ], 1, 4)
				# colnames(new_table) = colnames(coeftable)
				# rownames(new_table) = rownames(coeftable)[1]
				# class(new_table) = "coeftest"
				new_table = coeftable[1, , FALSE]
			} else {
				new_table = coeftable[-nrow(coeftable), ]
			}

			# stats::printCoefmat(new_table)
			myPrintCoefTable(new_table)
			cat("\nSigma =", coeftable[".sigma", 1], "\n")
		} else {
			# stats::printCoefmat(coeftable)
			myPrintCoefTable(coeftable)
		}

		cat("\n# Evaluations:", x$iterations, "\n")

	}

	bic <- -2*x$loglik+x$k*log(x$n)
	# aic <- -2*x$loglik+2*x$k
	# cat("Log-likelihood:", x$loglik, "\nBIC:", bic, "\nAIC:", aic, "\n")
	bic_format = ifelse(log10(abs(bic))<=9, addCommas(bic), bic)
	LL_format = ifelse(log10(abs(x$loglik))<=9, addCommas(x$loglik), x$loglik)

	cat("Log-likelihood:", LL_format, "\nBIC:", bic_format, "\n")
	cat("Pseudo-R2:", x$pseudo_r2, "\n")
	cat("Squared Cor.:", x$sq.cor, "\n")
	if(is.null(x$onlyCluster)) cat("Convergence state:", x$message, "\n")
}

##

#' Summary of a \code{femlm} object. Computes different types of standard errors.
#'
#' This function is similar to \code{print.femlm}. It provides the table of coefficients along with other information on the fit of the estimation. It can compute different types of standard errors. The new variance covariance matrix is an object returned.
#'
#' @method summary femlm
#'
#' @param se Character scalar. Which kind of standard error should be prompted: \dQuote{standard} (default), \dQuote{White}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}?
#' @param cluster A list of vectors. Used only if \code{se="cluster"}, \dQuote{se=twoway}, \dQuote{se=threeway} or \dQuote{se=fourway}. The vectors should give the cluster of each observation. Note that if the estimation was run using \code{cluster}, the standard error is automatically clustered along the cluster given in \code{\link[FENmlm]{femlm}}.
#' @param object A femlm object. Obtained using \code{\link[FENmlm]{femlm}}.
#' @param dof_correction Logical, default is \code{TRUE}. Should there be a degree of freedom correction to the standard errors of the coefficients?
#' @param forceCovariance Logical, default is \code{FALSE}. In the peculiar case where the obtained Hessian is not invertible (usually because of collinearity of some variables), use this option force the covariance matrix, by using a generalized inverse of the Hessian. This can be useful to spot where possible problems come from.
#' @param keepBounded Logical, default is \code{FALSE}. If \code{TRUE}, then the bounded coefficients (if any) are treated as unrestricted coefficients and their S.E. is computed (otherwise it is not).
#' @param ... Not currently used.
#'
#' @return
#' It returns a \code{femlm} object with:
#' \item{cov.scaled}{The new variance-covariance matrix (computed according to the argument \code{se}).}
#' \item{coeftable}{The table of coefficients with the new standard errors.}
#'
#' @seealso
#' See also the main estimation function \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # The data
#' n = 100
#' x = rnorm(n,1,5)**2
#' y = rnorm(n,-1,5)**2
#' z = rpois(n,x*y) + rpois(n, 2)
#' base = data.frame(x,y,z)
#'
#' # Comparing the results of a 'linear' function
#' est0L = femlm(z~log(x)+log(y), base, family="poisson")
#'
#' # Displaying the summary
#' summary(est0L, se="white")
#' myWhiteVcov = summary(est0L, se="white")$cov.scaled
#'
summary.femlm <- summary.feNmlm <- function(object, se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), cluster, dof_correction=FALSE, forceCovariance = FALSE, keepBounded = FALSE, ...){
	# computes the clustered SD and returns the modified vcov and coeftable
	# computes the BIC/AIC
	x = object

	if(!is.null(x$onlyCluster)){
		# means that the estimation is done without variables
		return(x)
	}

	dots = list(...)

	sd.val = match.arg(se)

	# If cov.scaled exists => means that it has already been computed
	if(!is.null(x$cov.scaled) && "fromPrint" %in% names(dots)) return(x)

	# => method to compute the clustered SE:
	quick = TRUE

	isBounded = object$isBounded

	if(is.null(isBounded)){
		isBounded = rep(FALSE, length(x$coef))
	}

	# We handle the bounded parameters:

	# keepBounded = dots$keepBounded
	# if(is.null(keepBounded)) keepBounded = FALSE

	if(any(isBounded)){
		if(keepBounded){
			# we treat the bounded parameters as regular variables
			myScore = x$score
			x$cov.unscaled = solve(x$hessian)
		} else {
			myScore = x$score[, -which(isBounded)]
		}
	} else {
		myScore = x$score
	}


	n = x$n
	k = x$k

	correction.dof = n / (n-k*dof_correction)

	x$bic <- -2*x$loglik + k*log(n)
	x$aic <- -2*x$loglik + 2*k

	if(anyNA(x$cov.unscaled)){
		if(!forceCovariance){
			warning("Standard errors set to NA because of collinearity. You can use option 'forceCovariance' to bypass this message.", call. = FALSE)
			return(x)
		} else {
			x$cov.unscaled = MASS::ginv(x$hessian)
			if(anyNA(x$cov.unscaled)) stop("The covariance matrix could not be 'forced'.")
			return(summary(x, se=sd.val, cluster=cluster, dof_correction=dof_correction))
		}
	} else if(sd.val == "standard"){
		vcov = x$cov.unscaled * correction.dof
	} else if(sd.val == "white"){
		vcov = crossprod(myScore%*%x$cov.unscaled) * correction.dof
	} else {
		# Clustered SD!
		nway = switch(sd.val, cluster=1, twoway=2, threeway=3, fourway=4)

		#
		# Controls
		#

		# Controlling the clusters
		do.unclass = TRUE
		if(missing(cluster)){
			if(is.null(x$id_dummies)){
				stop("To display clustered standard errors, you must provide the argument 'cluster'.")
			} else if(length(x$id_dummies) < nway) {
				stop("Since the result was not clustered with ", nway, " clusters, you need to provide the argument 'cluster' with ", nway, "clusters.")
			} else {
				cluster = x$id_dummies[1:nway]
				if(length(x$clusterNames)>1) warning("Standard-errors clustered w.r.t. ", paste0(x$clusterNames[1:nway], collapse = " & "), call. = FALSE, immediate. = TRUE)
				# in that specific case, there is no need of doing unclass.factor because already done
				do.unclass = FALSE
			}
		} else {

			isS = ifelse(nway>1, "s, each", "")
			if(! (is.list(cluster) && length(cluster)==nway) ) stop("The 'cluster' must be a list containing ", nway, " element", isS, " being the vector of IDs of each observation.")
			cluster = as.list(cluster)
		}

		# now we check the lengths:
		n_per_cluster = sapply(cluster, length)
		if(!all(diff(n_per_cluster) == 0)) stop("The vectors of the argument 'cluster' must be of the same length.")
		# Either they are of the same length of the data
		if(n_per_cluster[1] != x$n){
			# Then two cases: either the user introduces the original data and it is OK
			if(n_per_cluster[1] == (x$n + length(x$obsRemoved))){
				# We modify the clusters
				for(i in 1:nway) cluster[[i]] = cluster[[i]][-x$obsRemoved]
			} else {
				# If this is not the case: there is a problem
				stop("The length of the cluster does not match the original data.")
			}
		}

		#
		# Calculus
		#

		# initialisation
		vcov = x$cov.unscaled * 0

		if(!quick){
			for(i in 1:nway){

				myComb = combn(nway, i)

				for(j in 1:ncol(myComb)){
					if(i > 1){
						myDots = cluster[myComb[, j]]
						myDots$sep = "_"
						index = do.call(paste, myDots)
					} else {
						# for i==1, no need of the conversion to character
						index = cluster[[myComb[, j]]]
					}

					if(i > 1){
						do.unclass = TRUE
					}

					vcov = vcov + (-1)**(i+1) * vcovClust(index, x$cov.unscaled, myScore, dof_correction, do.unclass)
					# cat(sprintf("i=%i, j=%i\n", i, j))
					# print(vcovClust(index, x$cov.unscaled, myScore, dof_correction, do.unclass))
				}

			}
		} else {

			if(do.unclass){
				for(i in 1:nway){
					cluster[[i]] = quickUnclassFactor(cluster[[i]])
				}
			}

			for(i in 1:nway){

				myComb = combn(nway, i)

				power = floor(1 + log10(sapply(cluster, max)))

				for(j in 1:ncol(myComb)){

					if(i == 1){
						index = cluster[[myComb[, j]]]
					} else if(i > 1){

						vars = myComb[, j]

						if(sum(power[vars]) > 14){
							myDots = cluster[vars]
							myDots$sep = "_"
							index = do.call(paste, myDots)
						} else {
							# quicker, but limited by the precision of integers
							index = cluster[[vars[1]]]
							for(k in 2:length(vars)){
								index = index + cluster[[vars[k]]]*10**sum(power[vars[1:(k-1)]])
							}
						}

						index = quickUnclassFactor(index)

					}

					vcov = vcov + (-1)**(i+1) * vcovClust(index, x$cov.unscaled, myScore, dof_correction, do.unclass=FALSE)

				}
			}
		}
	}

	# Eigenfix if the vcov is not semi-positive definite
	if(any(diag(vcov)<0)){
		e = eigen(vcov)
		val = e$values
		if(is.complex(val)){
			# val = Re(val)
			warning("Variance is not positive definite. Some S.E. are NA.")
		} else {
			vect = e$vectors
			val[val<0] = 0
			vcov = vect%*%diag(val)%*%t(vect)
			warning("Variance was Eigenfixed.")
		}
	}

	sd2 = diag(vcov)
	sd2[sd2<0] = NA
	se = sqrt(sd2)

	# used to handle the case of bounded parameters
	params = names(x$coef)
	if(length(se) != length(params)){
		se = se[params]
	}
	names(se) = params

	# The coeftable is modified accordingly
	coeftable = x$coeftable

	# th z & p values
	zvalue <- x$coef/se
	pvalue <- 2*pnorm(-abs(zvalue))

	# update of se if bounded
	if(any(isBounded) & !keepBounded){
		se[!isBounded] = decimalFormat(se[!isBounded])
		boundText = attr(isBounded, "type")
		se[isBounded] = boundText
	}

	# modifs of the table
	coeftable[, 2] = se
	coeftable[, 3] = zvalue
	coeftable[, 4] = pvalue

	# coeftable[, 2:4] = cbind(se, x$coef/se, 2*pnorm(-abs(x$coef/se)))

	sd.dict = c("standard" = "Standard", "white"="White", "cluster"="Clustered", "twoway"="Two-way", "threeway"="Three-way", "fourway"="Four-way")
	attr(coeftable, "type") = attr(vcov, "type") = attr(se, "type") = as.vector(sd.dict[sd.val])
	x$cov.scaled = vcov
	x$coeftable = coeftable
	x$se = se

	return(x)
}

##

#' Facility to export the results of multiple \code{femlm} estimations in a Latex table.
#'
#' This function aggregates the results of multiple estimations and display them in the form of  one Latex table whose rownames are the variables and the columns contain the coefficients and standard-errors.
#'
#' @inheritParams summary.femlm
#'
#' @param ... Used to capture different \code{\link[FENmlm]{femlm}} objects. Note that any other type of element is discarded. Note that you can give a list of \code{\link[FENmlm]{femlm}} objects.
#' @param digits Integer. The number of digits to be displayed.
#' @param pseudo Logical, default is \code{TRUE}. Should the pseudo R2 be displayed?
#' @param title Character scalar. The title of the Latex table.
#' @param sdBelow Logical, default is \code{TRUE}. Should the standard-errors be displayed below the coefficients?
#' @param drop Character vector. This element is used if some variables are not to be displayed. This should be a regular expression (see \code{\link[base]{regexp}} help for more info). There can be more than one regular expression. Each variable satisfying the regular expression will be discarded.
#' @param order Character vector. This element is used if the user wants the variables to be ordered in a certain way. This should be a regular expression (see \code{\link[base]{regexp}} help for more info). There can be more than one regular expression. The variables satisfying the first regular expression will be placed first, then the order follows the sequence of regular expressions.
#' @param dict A named character vector. If provided, it changes the original variable names to the ones contained in the \code{dict}. Example: I want to change my variable named "a" to "$log(a)$" and "b3" to "$bonus^3$", then I used \code{dict=c(a="$log(a)$",b3="$bonus^3$")}.
#' @param file A character scalar. If provided, the Latex table will be saved in a file whose path is \code{file}.
#' @param append Logical, default is \code{TRUE}. Only used if option \code{file} is used. Should the Latex table be appended to the existing file?
#' @param convergence Logical, default is \code{TRUE}. Should the convergence state of the algorithm be displayed?
#' @param signifCode Named numeric vector, used to provide the significance codes with respect to the p-value of the coefficients. Default is \code{c("***"=0.01, "**"=0.05, "*"=0.10)}.
#' @param label Character scalar. The label of the Latex table.
#' @param aic Logical, default is \code{FALSE}. Should the AIC be displayed?
#' @param sqCor Logical, default is \code{FALSE}. Should the squared correlation be displayed?
#' @param subtitles Character vector of the same lenght as the number of models to be displayed. If provided, subtitles are added underneath the dependent variable name.
#' @param showClusterSize Logical, default is \code{FALSE}. If \code{TRUE} and clusters were used in the models, then the number "individuals" of per cluster is also displayed.
#' @param keepFactors Logical, default is \code{FALSE}. By default, when factor variables are contained in the estimation, they are printed as if they were a cluster variable. Put to \code{TRUE} to display all the coefficients of the factor variables.
#' @param bic Logical, default is \code{TRUE}.Should the BIC be reported?
#' @param loglik Logical, default is \code{TRUE}. Should the log-likelihood be reported?
#' @param yesNoCluster A character vector of lenght 2. Default is \code{c("Yes", "No")}. This is the message displayed when a given cluster is (or is not) included in a regression.
#'
#' @return
#' There is nothing returned, the result is only displayed on the console or saved in a file.
#'
#' @seealso
#' See also the main estimation function \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' n = 100
#' x = rnorm(n, 1, 5)**2
#' y = rnorm(n, -1, 5)**2
#' z = rpois(n, x*y) + rpois(n, 2)
#' base = data.frame(x, y, z)
#'
#' # Results of the Poisson..
#' est_poisson = femlm(z~log(x)+log(y), base, family="poisson")
#' # .. and of the Negative Binomial
#' est_negbin = femlm(z~log(x)+log(y), base, family="negbin")
#'
#' # We export the two results in one Latex table:
#' res2tex(est_poisson, est_negbin)
#'
#' # Changing the names & significance codes
#' res2tex(est_poisson, est_negbin, dict = c("log(x)" = "First variable (ln)"),
#'         signifCode = c("a" = 0.001, "$$" = 0.1))
#'
res2tex <- function(..., se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), cluster, digits=4, pseudo=TRUE, title, sdBelow=TRUE, drop, order, dict, file, append=TRUE, convergence=FALSE, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), label, aic=FALSE, sqCor = FALSE, subtitles, showClusterSize = FALSE, bic = TRUE, loglik = TRUE, yesNoCluster = c("Yes", "No"), keepFactors = FALSE){
	# drop: a vector of regular expressions
	# order: a vector of regular expressions
	# dict: a 'named' vector
	# file: a character string

	if(missing(se)){
		useSummary = FALSE
	} else {
		useSummary = TRUE
		sdType = match.arg(se)
	}

	info = results2formattedList(..., se=se, cluster=cluster, digits=digits, sdBelow=sdBelow, signifCode=signifCode, subtitles=subtitles, yesNoCluster=yesNoCluster, keepFactors=keepFactors, isTex = TRUE, useSummary=useSummary, sdType=sdType)

	# browser()

	n_models = length(info$depvar_list)
	# Getting the information
	se_type_list = info$se_type_list
	var_list = info$var_list
	coef_list = info$coef_list
	coef_below = info$coef_below
	sd_below = info$sd_below
	depvar_list = info$depvar_list
	obs_list = info$obs_list
	r2_list = info$r2_list
	aic_list = info$aic_list
	bic_list = info$bic_list
	loglik_list = info$loglik_list
	convergence_list = info$convergence_list
	sqCor_list = info$sqCor_list
	factorNames = info$factorNames
	isFactor = info$isFactor
	nbFactor = info$nbFactor

	if(!missing(subtitles)){
		isSubtitles = TRUE
	} else {
		isSubtitles = FALSE
	}

	#
	# prompting the infos gathered
	#

	# Starting the table
	myTitle = ifelse(!missing(title), title, "no title")
	if(!missing(label)) myTitle = paste0("\\label{", label, "} ", myTitle)
	start_table = paste0("\\begin{table}[htbp]\\centering\n\\caption{",  myTitle, "}\n")
	end_table = "\\end{table}"

	# intro and outro Latex tabular
	myAmpLine = paste0(paste0(rep(" ", length(depvar_list)+1), collapse="&"), "\\tabularnewline\n")
	intro_latex <- paste0("\\begin{tabular}{l", paste0(rep("c", n_models), collapse=""), "}\n",
								 myAmpLine,
								 "\\hline\n",
								 "\\hline\n")

	outro_latex <- "\\end{tabular}\n"

	# 1st lines => dep vars
	# first_line <- paste0("Variables&", paste0(depvar_list, collapse="&"), "\\\\\n\\hline\n\\hline\n")
	depvar_list = c(depvar_list, recursive = TRUE)
	if(!missing(dict)){
		if(!is.character(dict)|| is.null(names(dict))) stop("the arg. 'dict' must be a named character vector.")
		depvar_list = c(depvar_list, recursive = TRUE)
		qui = which(depvar_list%in%names(dict))
		who = depvar_list[qui]
		depvar_list[qui] = dict[who]
	}

	# We write the dependent variables properly, with multicolumn when necessary
	# to do that, we count the number of occurences of each variable (& we respect the order provided by the user)
	nb_multi = 1
	names_multi = depvar_list[1]

	if(n_models > 1){
		k = 1
		old_dep = depvar_list[1]
		for(i in 2:length(depvar_list)){
			if(depvar_list[i] == old_dep){
				nb_multi[k] = nb_multi[k] + 1
			} else {
				k = k + 1
				nb_multi[k] = 1
				names_multi[k] = old_dep = depvar_list[i]
			}
		}
	}

	# now the proper format
	first_line <- "Dependent Variables:"
	if(length(nb_multi) == 1) first_line = "Dependent Variable:"
	for(i in 1:length(nb_multi)){
		if(nb_multi[i] == 1){
			# no multi column
			first_line = paste0(first_line, "&", names_multi[i])
		} else {
			first_line = paste0(first_line, "&\\multicolumn{", nb_multi[i], "}{c}{", names_multi[i], "}")
		}
	}
	first_line = paste0(first_line, "\\\\\n")

	# Model line
	model_line = paste0("Model:&", paste0("(", 1:n_models, ")", collapse = "&"), "\\\\\n")

	# a simple line with only "variables" written in the first cell
	# variable_line = "Variables\\tabularnewline\n\\hline\n"
	variable_line = "\\hline\n\\emph{Variables}\\tabularnewline\n"


	# if(length(unique(depvar_list)) == 1){
	# 	first_line <- paste0("Dependent Variables:&\\multicolumn{",length(depvar_list), "}{c}{", depvar_list[1], "}\\\\\n", "Variables\\tabularnewline\n\\hline\n")
	# } else {
	# 	text_depvar = "Dependent Variables:&"
	# 	if(n_models == 1) text_depvar = "Dependent Variable:&"
	#
	# 	first_line <- paste0(text_depvar, paste0(depvar_list, collapse="&"), "\\\\\n",
	# 								"Variables\\tabularnewline\n\\hline\n")
	# }


	# Coefficients,  the tricky part
	coef_lines <- list()
	all_vars <- unique(c(var_list, recursive=TRUE))

	# dropping some coefs
	if(!missing(drop)){
		if(!is.character(drop)) stop("the arg. 'drop' must be a character vector of regular expression (see help regexp).")
		for(var2drop in drop) all_vars = all_vars[!grepl(var2drop, all_vars)]
	}

	# ordering the coefs
	if(!missing(order)){
		if(!is.character(order)) stop("the arg. 'order' must be a character vector of regular expression (see help regexp).")
		for(var2order in rev(order)){
			who = grepl(var2order, all_vars)
			all_vars = c(all_vars[who], all_vars[!who])
		}
	}

	# changing the names of the coefs
	aliasVars = all_vars
	names(aliasVars) = all_vars
	if(!missing(dict)){
		if(!is.character(dict)|| is.null(names(dict))) stop("the arg. 'dict' must be a named character vector.")
		# we add the overdispersion parameter
		dict[".theta"] = "Over-dispersion Parameter"
	} else {
		dict = c(".theta" = "Overdispersion Parameter")
	}

	# Changing the names
	qui = which(all_vars%in%names(dict))
	who = aliasVars[qui]
	aliasVars[qui] = dict[who]

	coef_mat <- all_vars
	for(m in 1:n_models) coef_mat <- cbind(coef_mat, coef_list[[m]][all_vars])
	coef_mat[is.na(coef_mat)] <- "  "
	if(sdBelow){
		coef_lines = c()
		for(v in all_vars){
			myCoef = mySd= myLine = c()
			for(m in 1:n_models){
				myCoef = c(myCoef, coef_below[[m]][v])
				mySd = c(mySd, sd_below[[m]][v])
			}
			myCoef[is.na(myCoef)] = "  "
			mySd[is.na(mySd)] = "  "
			myCoef = paste0(aliasVars[v], "&", paste0(myCoef, collapse="&"))
			mySd = paste0("  &", paste0(mySd, collapse="&"))
			myLines = paste0(myCoef, "\\\\\n", mySd, "\\\\\n")
			coef_lines = c(coef_lines, myLines)
		}
		coef_lines = paste0(coef_lines, collapse="")
	} else {
		coef_lines = paste0(paste0(apply(coef_mat, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")
	}

	# Factors (if needed)
	if(length(factorNames)>0){
		dumIntro = paste0("\\hline\n\\emph{Clusters}& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")

		for(m in 1:n_models) {
			quoi = isFactor[[m]][factorNames]
			quoi[is.na(quoi)] = yesNoCluster[2]
			isFactor[[m]] = quoi

			# We do the same for the number of items
			quoi = nbFactor[[m]][factorNames]
			quoi[is.na(quoi)] = "--"
			nbFactor[[m]] = quoi
		}

		allFactors = matrix(c(isFactor, recursive=TRUE), nrow = length(factorNames))
		# We change the names of the factors
		if(!missing(dict)){
			qui = which(factorNames %in% names(dict))
			factorNames[qui] = dict[factorNames[qui]]
		}

		allFactors = cbind(factorNames, allFactors)
		factor_lines <- paste0(paste0(apply(allFactors, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

		# For the number of items
		all_nb_Factors = matrix(c(nbFactor, recursive=TRUE), nrow = length(factorNames))
		factorNames_nbItems = paste0("# ", factorNames)
		all_nb_Factors = cbind(factorNames_nbItems, all_nb_Factors)
		nb_factor_lines <- paste0(paste0(apply(all_nb_Factors, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")


	} else {
		factor_lines = NULL
		dumIntro = NULL
	}

	# Subtitles
	if(isSubtitles){
		info_subtitles = paste0("  & ", paste(subtitles, collapse="&"), "\\\\\n")
	} else {
		info_subtitles = ""
	}

	# Fit statistics
	fit_part <- paste0("\\hline\n\\emph{Fit statistics}& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")
	# Misc
	# browser()
	# info_aic <- paste0("AIC & ", paste(addCommas(aic_list), collapse="&"), "\\\\\n")
	# info_loglik <- paste0("Log-Likelihood & ", paste(addCommas(loglik_list), collapse="&"), "\\\\\n")
	# info_bic <- paste0("BIC & ", paste(addCommas(bic_list), collapse="&"), "\\\\\n")
	info_aic <- paste0("AIC & ", paste(numberFormatLatex(aic_list), collapse="&"), "\\\\\n")
	info_loglik <- paste0("Log-Likelihood & ", paste(numberFormatLatex(loglik_list), collapse="&"), "\\\\\n")
	info_bic <- paste0("BIC & ", paste(numberFormatLatex(bic_list), collapse="&"), "\\\\\n")

	info_obs <- paste0("Observations& ", paste(addCommas(obs_list), collapse="&"), "\\\\\n")
	info_r2 <- paste0("Adj-pseudo $R^2$ &", paste(r2_list, collapse="&"), "\\\\\n")
	info_sqCor <- paste0("$R^2$ &", paste(sqCor_list, collapse="&"), "\\\\\n")
	info_convergence = paste0("Convergence &", paste(convergence_list, collapse="&"), "\\\\\n")

	# The standard errors
	isUniqueSD = length(unique(unlist(se_type_list))) == 1
	if(isUniqueSD){
		my_se = unique(unlist(se_type_list)) # it comes from summary
		# every model has the same type of SE
		nameSD = c("Standard"="Normal", "White"="White-corrected", "Clustered"="Clustered", "Two-way"="Two-way clustered")
		nb_col = length(obs_list)+1
		info_SD = paste0("\\hline\n\\hline\n\\multicolumn{", nb_col, "}{l}{\\emph{", nameSD[my_se], " standard-errors in parenthesis. Signif Codes: ", paste(names(signifCode), signifCode, sep=": ", collapse = ", "), "}}\\\\\n")
		info_muli_se = ""
	} else {
		info_muli_se = paste0("Standard-Error type& ", paste(se_type_list, collapse="&"), "\\\\\n")
		info_SD = "\\hline\n\\hline\n\\\\\n"
	}

	# Information on number of items

	supplemental_info = "\\global\\long\\def\\sym#1{\\ifmmode^{#1}\\else\\ensuremath{^{#1}}\\fi}\n"

	if(!pseudo) info_r2 <- ""
	if(!sqCor) info_sqCor <- ""
	if(!convergence) info_convergence = ""
	if(!aic) info_aic = ""
	if(!bic) info_bic = ""
	if(!loglik) info_loglik = ""
	if(!showClusterSize) nb_factor_lines = ""

	if(!missing(file)) sink(file = file, append = append)

	cat(paste0(supplemental_info,
				  start_table,
				  intro_latex,
				  first_line,
				  info_subtitles,
				  model_line,
				  variable_line,
				  coef_lines,
				  dumIntro,
				  factor_lines,
				  fit_part,
				  info_obs,
				  nb_factor_lines,
				  info_convergence,
				  info_muli_se,
				  info_r2,
				  info_sqCor,
				  info_aic,
				  info_loglik,
				  info_bic,
				  info_SD,
				  outro_latex,
				  end_table))

	if(!missing(file)) sink()

}

#' Facility to display the results of multiple \code{femlm} estimations.
#'
#' This function aggregates the results of multiple estimations and display them in the form of only one table whose rownames are the variables and the columns contain the coefficients and standard-errors.
#'
#' @inheritParams res2tex
#' @inheritParams summary.femlm
#'
#' @return
#' Returns a data.frame containing the formatted results.
#'
#' @seealso
#' See also the main estimation function \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#' n = 100
#' x = rnorm(n, 1, 5)**2
#' y = rnorm(n, -1, 5)**2
#' z = rpois(n, x*y) + rpois(n, 2)
#' base = data.frame(x, y, z)
#'
#' # Results of the Poisson..
#' est_poisson = femlm(z~log(x)+log(y), base, family="poisson")
#' # .. and of the Negative Binomial
#' est_negbin = femlm(z~log(x)+log(y), base, family="negbin")
#'
#' # We export the two results in one Latex table:
#' res2table(est_poisson, est_negbin)
#'
#' # Changing the names & significance codes
#' res2table(est_poisson, est_negbin, dict = c("log(x)" = "First variable (ln)"),
#'         signifCode = c("a" = 0.001, "$$" = 0.1))
#'
#'
res2table <- function(..., se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), cluster, digits=4, pseudo=TRUE, drop, order, convergence=FALSE, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), subtitles, keepFactors = FALSE){

	# Need to check for the presence of the se
	if(missing(se)){
		useSummary = FALSE
	} else {
		useSummary = TRUE
		sdType = match.arg(se)
	}

	info = results2formattedList(..., se=se, cluster=cluster, digits=digits, signifCode=signifCode, subtitles=subtitles, keepFactors=keepFactors, useSummary=useSummary, sdType=sdType)

	n_models = length(info$depvar_list)
	# Getting the information
	se_type_list = info$se_type_list
	var_list = info$var_list
	coef_list = info$coef_list
	coef_below = info$coef_below
	sd_below = info$sd_below
	depvar_list = info$depvar_list
	obs_list = info$obs_list
	r2_list = info$r2_list
	aic_list = info$aic_list
	bic_list = info$bic_list
	loglik_list = info$loglik_list
	convergence_list = info$convergence_list
	sqCor_list = info$sqCor_list
	factorNames = info$factorNames
	isFactor = info$isFactor
	nbFactor = info$nbFactor
	useSummary = info$useSummary

	if(!missing(subtitles)){
		isSubtitles = TRUE
	} else {
		isSubtitles = FALSE
	}

	# The coefficients

	all_vars <- unique(c(var_list, recursive=TRUE))

	# dropping some coefs
	if(!missing(drop)){
		if(!is.character(drop)) stop("the arg. 'drop' must be a character vector of regular expression (see help regexp).")
		for(var2drop in drop) all_vars = all_vars[!grepl(var2drop, all_vars)]
	}

	# ordering the coefs
	if(!missing(order)){
		if(!is.character(order)) stop("the arg. 'order' must be a character vector of regular expression (see help regexp).")
		for(var2order in rev(order)){
			who = grepl(var2order, all_vars)
			all_vars = c(all_vars[who], all_vars[!who])
		}
	}

	coef_mat <- all_vars
	for(m in 1:n_models) coef_mat <- cbind(coef_mat, coef_list[[m]][all_vars])
	coef_mat[is.na(coef_mat)] <- "  "
	res = coef_mat

	# Used to draw a line
	myLine = "______________________________________"
	longueur = apply(res, 2, function(x) max(nchar(as.character(x))))
	theLine = sapply(longueur, function(x) sprintf("%.*s", x, myLine))

	if(length(factorNames)>0){
		# TO REWRITE: from old code
		for(m in 1:n_models) {
			quoi = isFactor[[m]][factorNames]
			quoi[is.na(quoi)] = "No"
			isFactor[[m]] = quoi
		}
		allFactors = matrix(c(isFactor, recursive=TRUE), nrow = length(factorNames))
		allFactors = cbind(factorNames, allFactors)
		factor_lines <- paste0(paste0(apply(allFactors, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

		myLine = "-------------------------------"
		# res = rbind(res, theLine)
		res = rbind(res, c("Clusters:", sprintf("%.*s", longueur[-1], myLine)))
		factmat = matrix(c(strsplit(strsplit(factor_lines, "\n")[[1]], "&"), recursive = TRUE), ncol=n_models+1, byrow=TRUE)
		factmat[, ncol(factmat)]=gsub("\\", "", factmat[, ncol(factmat)], fixed = TRUE)
		res = rbind(res, factmat)
	}

	res <- rbind(res, theLine)
	res <- rbind(res, c("Observations", addCommas(obs_list)))
	if(!useSummary) res <- rbind(res, c("S.E. type", c(se_type_list, recursive = TRUE)))
	if(convergence) res <- rbind(res, c("Convergence", c(convergence_list, recursive = TRUE)))
	res <- rbind(res, c("Squared-Correlation", c(sqCor_list, recursive = TRUE)))
	res <- rbind(res, c("Adj-pseudo R2", c(r2_list, recursive = TRUE)))
	# res <- rbind(res, c("AIC", c(aic_list, recursive = TRUE)))
	res <- rbind(res, c("Log-Likelihood", numberFormatNormal(loglik_list)))
	res <- rbind(res, c("BIC", numberFormatNormal(bic_list)))

	# if subtitles
	if(isSubtitles){
		modelNames = subtitles
	} else {
		modelNames = paste0("model ", 1:n_models)
	}

	res <- as.data.frame(res)
	names(res) <- c("variables", modelNames)
	row.names(res) = res$variables
	res$variables = NULL

	# We rename theta when NB is used
	quiTheta = which(row.names(res) == ".theta")
	row.names(res)[quiTheta] = "Dispersion Parameter"

	return(res)
}

results2formattedList = function(..., se=c("standard", "white", "cluster", "twoway"), cluster, digits=4, pseudo=TRUE, sdBelow=TRUE, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), label, subtitles, yesNoCluster = c("Yes", "No"), keepFactors = FALSE, isTex = FALSE, useSummary, sdType){
	# This function is the core of the functions res2table and res2tex

	signifCode = sort(signifCode)
	if(any(signifCode<0) | any(signifCode>1)) stop("The argument 'signifCode' must lie between 0 and 1.")

	if(length(yesNoCluster) != 2) stop("The argument 'yesNoCluster' must be of length 2.")

	# To take care of old verions:
	allowed_types = c("femlm", "feNmlm")

	# We get all the models
	dots <- list(...)

	# for retro-compatibility:
	if("sd" %in% names(dots)){
		warning("The use of the argument 'sd' is deprecated, it is now replaced by the argument 'se'.")
		se = dots$sd
	}

	n = length(dots)
	all_models = list()
	k = 1
	for(i in 1:n){
		di = dots[[i]]

		if(any(allowed_types %in% class(di))){
			all_models[[k]] = di
			k = k+1
		} else if(length(class(di))==1 && class(di)=="list"){
			# we get into this list to get the FENmlm objects
			types = sapply(di, class)
			qui = which(types %in% allowed_types)
			for(m in qui){
				all_models[[k]] = di[[m]]
				k = k+1
			}
		}

	}

	if(length(all_models)==0) stop("Not any proper model (femlm) as argument!")

	n_models <- length(all_models)

	# we keep track of the SEs
	se_type_list = list()

	var_list <- coef_list <- coef_below <- sd_below <- list()
	depvar_list <- obs_list <- list()
	r2_list <- aic_list <- bic_list <- loglik_list <- convergence_list <- list()
	sqCor_list = list()

	# To take care of factors
	factorNames = c()
	isFactor = vector(mode = "list", n_models)
	nbFactor = vector(mode = "list", n_models) # the number of items per factor

	# if there are subtitles
	if(!missing(subtitles)){
		if(length(subtitles) != n_models){
			stop("If argument 'subtitles' is provided, it must be of the same length as the number of models.")
		} else {
			isSubtitles = TRUE
		}
	} else {
		isSubtitles = FALSE
	}

	for(m in 1:n_models){

		# If se is provided, we use summary
		if(useSummary){
			x = summary(all_models[[m]], se=sdType, cluster)
		} else {
			x = all_models[[m]]
		}
		se_type_list[[m]] = attr(x$se, "type")

		# variable dependante:
		depvar <- gsub(" ", "", as.character(x$fml)[[2]])

		a <- x$coeftable
		if(!is.data.frame(a)){
			class(a) <- NULL
			a = as.data.frame(a)
		}

		#
		# Formatting of the factors
		#

		# on enleve les facteurs des variables a garder
		if(!keepFactors){
			fact = rownames(a)
			qui_drop = grepl("factor(", fact, fixed = TRUE)
			a = a[!qui_drop, , FALSE]
			b = fact[qui_drop]
			c = sapply(b, function(x) strsplit(x, "factor(", fixed=TRUE)[[1]][2])
			d = sapply(c, function(x) strsplit(x, ")", fixed=TRUE)[[1]][1])
			factor_var = unique(d)

			# Now the number of items per factor
			if(length(factor_var) == 0){
				nbItems = character(0)
			} else {
				nbItems = addCommas(sapply(factor_var, function(x) 1+sum(grepl(x, b))))
			}
		} else {
			factor_var = c()
			nbItems = character(0)
		}

		# now the normal clusters
		if(!is.null(x$clusterNames)){
			factor_var = c(factor_var, x$clusterNames, recursive=TRUE)

			if(!is.null(x$clusterSize)){
				new_items = addCommas(as.vector(x$clusterSize))
				names(new_items) = names(x$clusterSize)
			} else {
				# for retro compatibility
				new_items = addCommas(sapply(x$id_dummies, function(y) length(unique(y))))
			}

			nbItems = c(nbItems, new_items)
		}

		nbFactor[[m]] = nbItems

		# Formatting

		lFactor = rep(yesNoCluster[1], length(factor_var))
		names(lFactor) = factor_var
		isFactor[[m]] = lFactor

		factorNames = unique(c(factorNames, factor_var, recursive=TRUE))

		#
		# END: cluster formatting
		#

		# on enleve les espaces dans les noms de variables
		var <- c(gsub(" ", "", row.names(a)))

		coef = as.character(round(a[, 1], digits))
		se = as.character(myRound(a[, 2], digits))
		if(isTex){
			pval = cut(a[, 4], breaks = c(-1, signifCode, 100), labels = c(paste0("\\sym{",names(signifCode),"}"), ""))
		} else {
			pval = cut(a[, 4], breaks = c(-1, signifCode, 100), labels = c(names(signifCode), ""))
		}

		# If the coefficient is bounded, we supress the 'stars'
		isBounded = grepl("bounded", se)
		if(any(isBounded)){
			pval[isBounded] = ""
		}

		structured_coef = c(paste0(coef, pval, " (", se, ")"))

		# saving the infos
		var_list[[m]] <- var
		names(structured_coef) <- var
		coef_list[[m]] <- structured_coef
		if(sdBelow){
			cb = c(paste0(coef, pval))
			sb = c(paste0("(", se, ")"))
			names(cb) = names(sb) = var
			coef_below[[m]] = cb
			sd_below[[m]] = sb
		}

		# La depvar
		depvar_list[[m]] <- depvar

		# statistics
		# Pseudo-R2 // AIC // BIC // N
		n <- x$n
		obs_list[[m]] <- n
		convergence_list[[m]] = x$convStatus

		K <- x$k #nb params
		ll <- x$loglik
		bic_list[[m]] <- round(-2*ll+K*log(n), 3)
		aic_list[[m]] <- round(-2*ll+2*K, 3)
		loglik_list[[m]] <- round(ll, 3)
		r2_list[[m]] <- round(x$pseudo_r2, 5)
		sqCor_list[[m]] <- round(x$sq.cor, 3)

	}

	res = list(se_type_list=se_type_list, var_list=var_list, coef_list=coef_list, coef_below=coef_below, sd_below=sd_below, depvar_list=depvar_list, obs_list=obs_list, r2_list=r2_list, aic_list=aic_list, bic_list=bic_list, loglik_list=loglik_list, convergence_list=convergence_list, sqCor_list=sqCor_list, factorNames=factorNames, isFactor=isFactor, nbFactor=nbFactor, useSummary=useSummary)

	return(res)
}

myPrintCoefTable = function(coeftable){
	# Simple function that does as the function coeftable but handles special cases
	# => to take care of the case when the coefficient is bounded

	if(!is.data.frame(coeftable)){
		class(coeftable) = NULL
		ct = as.data.frame(coeftable)
	} else {
		ct = coeftable
	}

	signifCode = c("***"=0.001, "** "=0.01, "*  "=0.05, ".  "=0.1)

	pvalues = ct[, 4]

	stars = cut(pvalues, breaks = c(-1, signifCode, 100), labels = c(names(signifCode), ""))
	stars[is.na(stars)] = ""

	whoIsLow = !is.na(pvalues) & pvalues < 2.2e-16

	for(i in 1:4){
		ct[, i] = decimalFormat(ct[, i])
	}

	# browser()

	ct[whoIsLow, 4] = "< 2.2e-16"
	ct[is.na(ct[, 4]), 4] = "NA"

	ct[, 5] = stars
	names(ct)[5] = ""

	print(ct)

	cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")

}


print.femlm.cluster = function(x, ...){
	# just to hide the information on the class and on the attributes
	print.default(x[1:length(x)], ...)
}


#' Extract the Fixed-Effects from a \code{femlm} estimation.
#'
#' This function retrives the fixed effects from a femlm estimation. It is useful only when there are more than one cluster.
#'
#' @param x A femlm object.
#'
#' @return
#' A list containig the vectors of the fixed effects.
#'
#' @seealso
#' See also the main estimation function \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Bilateral network
#' nb = 20
#' n = nb**2
#' k = nb
#' id1 = factor(rep(1:k, each=n/k))
#' id2 = factor(rep(1:(n/k), times=k))
#' d = rep(rnorm(k)**2, each=n/k)
#' x = rnorm(n, 1, 5)**2
#' y = rnorm(n, -1, 5)**2
#' z = rpois(n, x*y+rnorm(n, sd = 3)**2) + rpois(n, 2)
#' base = data.frame(x, y, z, id1, id2)
#'
#' # We want to use the ID's of each observation as a variable: we use the option cluster
#' est_poisson = femlm(z~log(x)+log(y), base, family="poisson", cluster=c("id1", "id2"))
#'
#' # To get the FE:
#' myFE = getFE(est_poisson)
#'
getFE = function(x){
	# x is a femlm object
	# This function retrieves the dummies

	# Preliminary stuff
	S = x$dummies
	if(is.null(S)) stop("There is no cluster to be retrieved.")

	family = x$family
	clustNames = x$clusterNames

	id_dummies = x$id_dummies

	Q = length(id_dummies)
	N = length(S)

	# either (we need to clean its attributes for unlist to be efficient)
	id_dummies_vect = list()
	for(i in 1:Q) id_dummies_vect[[i]] = as.vector(id_dummies[[i]])

	if(Q == 1){
		# This is the simplest case
		id = id_dummies_vect[[1]]

		myOrder = order(id)
		myDiff = c(1, diff(id[myOrder]))

		select = myOrder[myDiff == 1]

		cluster_values = list(S[select])

	} else {
		# We apply a Rcpp script to handle complicated cases (and we don't know beforehand if there are some)

		dumMat <- matrix(unlist(id_dummies_vect), N, Q) - 1
		orderCluster <- matrix(unlist(lapply(id_dummies_vect, order)), N, Q) - 1

		nbCluster = sapply(id_dummies, max)

		cluster_values <- RcppGetFE(Q, N, S, dumMat, nbCluster, orderCluster)
	}

	# now saving & adding information
	all_clust = list()
	for(i in 1:Q){
		cv = cluster_values[[i]]
		names(cv) = attr(id_dummies[[i]], "clust_names")
		class(cv) = "femlm.cluster"
		attr(cv, "name") = clustNames[i]
		attr(cv, "family") = family
		all_clust[[clustNames[i]]] = cv
	}

	class(all_clust) = "femlm.allClusters"

	# # Construction de S et calcul de la diff maximale
	# S_new = S*0
	# for(i in 1:Q) S_new = S_new + all_clust[[i]][id_dummies[[i]]]
	#
	# print(sprintf("Max diff(S, S_new) = %0.6f", sum(abs(S - S_new))))

	return(all_clust)
}



plot.femlm.cluster = function(x, n=5, ...){
	# It plots the n first and last most notable FEs

	# meta information
	clusterName = attr(x, "name")
	family = attr(x, "family")

	# we compare with the average of the coefficients
	x = sort(x) - mean(x)
	x_name = names(x)

	k = length(x)
	nb_show = min(k, 2*n+1)
	mid_value = floor(nb_show/2)
	xlim = c(1, nb_show)
	xlim = xlim + c(-1, 1) * diff(xlim)/30

	if(k <= 2*n+1){
		# no need of space to split the data
		x_values = 1:k
		y_values = x
		isSplit = FALSE
	} else {
		nb_show = nb_show - 1 # because we don't want to show the middle point
		x_values = c(1:n, (n+2):(2*n+1))
		y_values = c(head(x, n), tail(x, n))
		isSplit = TRUE
	}

	# very specific case where the axis confonds with the boxing
	ylim = range(y_values)
	if(pretty(range(y_values))[1] < min(y_values) || tail(pretty(range(y_values)), 1) > max(y_values)){
		ylim = range(y_values) + c(-1, 1) * diff(range(y_values))/30
	}

	plot(x_values, y_values, ann = FALSE, axes = FALSE, xlim=xlim, ylim = ylim, col = 0)

	# display
	box()
	y =axis(2)
	abline(h = y, lty = 4, col = "lightgrey")

	# name information & points
	points(x_values, y_values)
	text(head(x_values, mid_value), head(y_values, mid_value), head(x_name, mid_value), pos = 3)
	text(tail(x_values, nb_show-mid_value), tail(y_values, nb_show-mid_value), tail(x_name, nb_show-mid_value), pos = 1)

	axis(4, at=y, labels = signif(exp(y), 2))

	# browser()
	title(main = clusterName)
	title(xlab = "Centered Cluster Coefficients", line = 1)

	if(isSplit){
		abline(v = c(n+0.75, n+1.25), lty = 2)

		axis(1, at = c(n+0.75, n+1.25), col = "grey", labels = NA, lwd.ticks = 0)
		axis(3, at = c(n+0.75, n+1.25), col = "grey", labels = NA, lwd.ticks = 0)
	}

	coord = par("usr")
	axis(1, at = coord[1:2], labels = c("coef", "exp(coef)"))

}


plot.femlm.allClusters = function(x, ...){

	n = length(x)

	mfrow = as.character(c(11, 12, 22, 22, 32, 32, 33, 33))

	op = par(mfrow = as.numeric(strsplit(mfrow[n], "")[[1]]), mar = c(3, 3, 2.5, 3))

	for(i in 1:n){
		plot(x[[i]])
	}

	par(op)
}

#
# To compute clustered standard errors
#

vcovClust <- function (cluster, myBread, scores, dof_correction=FALSE, do.unclass=TRUE){
	# Internal function: no need for controls, they come beforehand
	#	- cluster: the vector of dummies
	#	- myBread: the naive variance covariance matrix
	# - scores
	#Note: if length(unique(cluster)) == n (i.e. White correction), then the dof are such that vcovClust is equivalent to vcovHC(res, type="HC1")

	n <- NROW(scores)
	k <- NCOL(scores)

	# Control for cluster type
	if(do.unclass){
		cluster <- unclass(as.factor(cluster))
	}

	Q <- max(cluster)
	#ind: the matrix used to create the sum of the scores per cluster (avoid loops)
	# ind = Matrix::Matrix(0, nid, n, sparse = TRUE)
	# myIndex = cbind(cluster, 1:n)
	# ind[myIndex] = 1
	#
	# # Compute the chocolate_bar and then the sandwich:
	# RightScores <- as.matrix(ind %*% scores)
	# # chocobar <- crossprod(RightScores)

	RightScores = cpp_tapply_sum(Q, scores, cluster)

	# Finite sample correction:
	if(dof_correction) dof  <- Q / (Q-1) * (n-1) / (n-k)
	else dof = 1

	return(crossprod(RightScores%*%myBread) * dof)

	# return(myBread %*% chocobar %*% myBread * dof)
}

addCommas_single = function(x){

	if (!is.finite(x)) return(as.character(x))

	s = sign(x)
	x = abs(x)
	decimal = x - floor(x)
	if (decimal > 0){
		dec_string = substr(decimal, 2, 4)
	} else {
		dec_string = ""
	}

	entier = sprintf("%.0f", floor(x))
	quoi = rev(strsplit(entier, "")[[1]])
	n = length(quoi)
	sol = c()
	for (i in 1:n) {
		sol = c(sol, quoi[i])
		if (i%%3 == 0 && i != n) sol = c(sol, ",")
	}
	res = paste0(ifelse(s == -1, "-", ""), paste0(rev(sol), collapse = ""),
					 dec_string)
	res
}

addCommas = function(x){
	sapply(x, addCommas_single)
}

myRound_single = function(x, digits=5){
	# There can be non numeric values...
	# we give away the non numeric ones and round the others

	if(is.na(x)){
		return(NA)
	}

	if(is.numeric(x)){
		res = round(x, digits)
	} else {

		if(!grepl("[[:digit:]]", x)){
			# means it is a character
			res = x
		} else {
			res = round(as.numeric(x), digits)
		}
	}

	res
}

myRound = function(x, digits=5){
	sapply(x, myRound_single, digits = digits)
}

decimalFormat_single = function(x){
	# for very small numbers: format 5.2e-08

	if(is.na(x) || !is.numeric(x)) return(x)

	xPower = log10(abs(x))

	if(xPower < -5){
		res = signif(x, 3)
	} else if(xPower < 0){
		res = round(x, 6)
	} else {
		res = round(x, max(1, 5 - ceiling(xPower)))
	}

	res
}

decimalFormat = function(x){
	sapply(x, decimalFormat_single)
}


numberFormat_single = function(x, type = "normal"){

	exponent = floor(log10(abs(x)))

	if(exponent<9) return(addCommas(x))

	# For numbers higher than 1e9 => we apply a specific formatting

	left_value = round(x*10**-exponent, 3)

	if(type == "latex"){
		res = paste0("$", left_value, "\\times 10^{", exponent, "}$")
	} else {
		res = paste0(left_value, "e+", exponent)
	}

	res
}

numberFormatLatex = function(x){
	sapply(x, numberFormat_single, type = "latex")
}

numberFormatNormal = function(x){
	sapply(x, numberFormat_single)
}


char2num = function(x){
	# we transform the data to numeric => faster analysis

	x_unik = unique(x)
	dict = 1:length(x_unik)
	names(dict) = x_unik
	x_num = dict[x]

	names(x_num) = NULL

	x_num
}

quickUnclassFactor = function(x){
	# does as unclass(as.factor(x))
	# but waaaaay quicker

	if(is.factor(x)){
		x = unclass(x)
	} else if(is.character(x)){
		res = char2num(x)
		return(res)
	}

	stopifnot(is.numeric(x))

	myOrder = order(x)
	x_sorted = x[myOrder]
	g = Rcpp_unclassFactor(x_sorted)
	res = g[order(myOrder)]

	return(res)
}


getItems = function(x){
	# to get the unique elements of x before quickunclassfactor
	# needs to be done because differs depending on the type of x

	if(is.character(x)){
		res = unique(x)
	} else if(is.factor(x)){
		res = levels(unique(x)[, drop=TRUE])
	} else {
		res = sort(unique(x))
	}

	return(res)
}

extractCluster = function(fml){
	# We extract the clusters (ie after the |, if there are any)

	x = as.character(fml)[3]

	x_split = strsplit(x, split = "|", fixed = TRUE)[[1]]

	# To return:
	# - cluster => vector of clusters
	# - fml_new => the fml clean of the clusters

	if(length(x_split) == 1){
		fml_new = fml
		cluster_vec = NULL
	} else if(length(x_split) == 2){
		# means it's OK (more there is a problem), less there is none
		# check
		cluster = x_split[2]
		if(grepl("\\(|\\)|\\*", cluster)){
			stop("If the formula shall contains the clusters, it should only consists of the variables names.\n(eg 'a~b|cluster_1+cluster_2' is OK, but 'a~b|5*cluster_1+factor(cluster_2)' reports an error.)")
		}

		fml_cluster = as.formula(paste0("~", cluster))
		cluster_vec = all.vars(fml_cluster)

		fml_new = as.formula(paste0(fml[2], "~", x_split[1]))
	} else {
		stop("The formula must not contain more than one pipe ('|').")
	}

	list(fml=fml_new, cluster=cluster_vec)
}


#' Trade data sample
#'
#' This data reports trade information between countries of the European Union (EU15).
#'
#' @usage
#' data(trade)
#'
#' @format
#' \code{trade} is a data frame with 38,325 observations and 6 variables named \code{Destination}, \code{Origin}, \code{Product}, \code{Year}, \code{dist_km} and \code{Euros}.
#'
#' \itemize{
#' \item{Origin}{2-digits codes of the countries of origin of the trade flow.}
#' \item{Destination}{2-digits codes of the countries of destination of the trade flow.}
#' \item{Products}{Number representing the product categories (from 1 to 20).}
#' \item{Year}{Years from 2007 to 2016}
#' \item{dist_km}{Geographic distance in km between the centers of the countries of origin and destination.}
#' \item{Euros}{The total amount in euros of the trade flow for the specific year/product category/origin-destination country pair.}
#'
#' }
#'
#' @source
#' This data has been extrated from Eurostat on October 2017.
#'
#'
#'
"trade"





































































