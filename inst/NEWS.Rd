\name{NEWS}
\alias{NEWS}

\title{News for \R Package \pkg{FENmlm}}

\section{Changes in version 2.2.0 (2018-08-03)}{
	\subsection{New features}{
		\itemize{
			\item [stats.methods] Several statistics methods have been implemented: nobs, AIB, BIC, logLik, coef, fitted, residuals, formula, model.matrix, predict, vcov, confint and update.
			\item [res2table] New arguments \sQuote{depvar} and \sQuote{family}: to control whether the dependent variables (resp. model families) should be displayed.
			\item [res2table] Now the \code{call} is used to produce more informative model names by default.
			\item [res2tex] New argument \sQuote{family} to control whether the families of the models should be displayed.
			\item [summary/vcov] Added the possibility to cluster the standard-errors with direct reference to a cluster variable used in the estimation. Example: \code{summary(est, se = "cluster", cluster = "var_2")}.
		}
	}
	\subsection{Improvements}{
		\itemize{
			\item [femlm] Formulas are more robustly handled with the package Formula.
			\item [femlm--negative binomial] Added a warning message when there is no overdispersion.
			\item [print.femlm] Information is now displayed more compactly
			\item Slightly faster internal functions.
		}
	}
	\subsection{Changing arguments}{
		\itemize{
			\item Some argument names in femlm have been modified for clarification:
			\itemize{
				\item [start => NL.start] Argument \code{start} now becomes \code{NL.start} to emphasize that it concerns the non-linear in parameters part only. This is the same for argument \code{start.init} which becomes \code{NL.start.init}.
				\item [debug => verbose] Argument \code{debug} of type logical now becomes \code{verbose} of type integer which can take values 0, 1 or 2 depending of the level of information to be reported.
			}
		}
	}
	\subsection{Documentation}{
		\itemize{
			\item Improved examples.
			\item A new reference is added which details the theory behind the method.
		}
	}
	\subsection{Bug}{
		\itemize{
			\item [femlm--data] Now the function properly handles matrices in argument \sQuote{data}.
		}
	}
	\subsection{Miscellaneous}{
		\itemize{
			\item OpenMP is handled differently (now used conditionnaly).
			\item CITATION is added.
			\item Some items in the femlm output have changed names to fit the \code{lm/glm} vocabuary more closely. This concerns: (old name => new name) coef => coefficients, expected.predictor => fitted.values, k => nparams, resids => residuals.
			\item A demo has been added.
			\item Warnings have been suppressed when the clustering of standard-errors was equivocal. They have been replaced by explicit mention of the cluster variables in the type of the standard-errors.
		}
	}
}

\section{Changes in version 2.1.2}{
	\subsection{Improvements}{
		\itemize{
			\item [femlm] Improved warning messages for overfitting problems.
			\item [femlm] Improved error-handling.
		}
	}
	\subsection{Documentation}{
		\itemize{
			\item Amelioration of the Vignette.
		}
	}
}

\section{Changes in version 2.1.1}{
	\subsection{New features}{
		\itemize{
			\item [femlm] Added argument \dQuote{precision.cluster} so that the user can set the level of precision manually (concerns the algorithm to obtain the cluster coefficients).
			\item [summary.femlm] For one-way clustering, the argument \dQuote{cluster} can now be a vector.
			\item A print method for the fixed-effects obtained with function \code{\link[FENmlm]{getFE}} has been implemented.
		}
	}
	\subsection{Documentation}{
		\itemize{
			\item A vignette has been created.
			\item Documentation has been added for fixed-effects plotting and for a new print method.
			\item Some examples have been updated.
		}
	}
	\subsection{Bug}{
		\itemize{
			\item [getFE] Corrected a small glitch that impeded the algorithm to find the right solution in presence of  complicated cases.
			\item [femlm] The function could not be evaluated at all when parallel::detectCores() gives NA. This has been corrected. Now when parallel::detectCores() gives NA, the number of cores cannot be set greater than 1.
		}
	}
	\subsection{Miscellaneous}{
		\itemize{
			\item [res2tex] In the exported Latex table the word \sQuote{Cluster} is changed to \sQuote{Fixed-Effects}.
			\item [summary.femlm] Using cluster=NULL does not cause problems any more.
			\item [print.femlm] Corrected a spelling mistake.
		}
	}

}

\section{Changes in version 2.1.0}{
	\subsection{New Features}{
		\itemize{
			\item New argument 'cores' => added multi cores with openMP. Particularly efficient for the Negative Binomial and the Logit models.
			\item Added a fixed-point acceleration algorithm for 3+ clusters and also to obtain the derivatives.
			\item Cluster coefficient precision is dynnamically adjusted for complex problems.
			\item  Faster matrix creation.
		}
	}
  \subsection{Other}{
		\itemize{
			\item [Design] Argument noWarning is dropped.
			\item [Design] clearer message when clusters are removed.
			\item [Bug] corrected undue information reporting.
			\item [Bug] corrected important bug that made some estimations fail.
			\item [Help] Added more accurate information on the returned objects.
		}
	}

}


\section{Changes in version 2.0.1}{
	\itemize{
	  \item [c++ code] Fixed overloading problem.
	  \item Improved help pages.
	  \item Code more robust to overfitting situations (in particular for the Logit case).
	  \item Added an acceleration algorithm to get the fixed-point for 2 clusters.
	  \item New collinearity warnings.
    }
  }

\section{Changes in version 2.0}{
  \subsection{New Features}{
    \itemize{
	  \item Algorithm is greatly improved.
      \item Two new likelihood models: Negative Binomial and Gaussian.
      \item Any number of clusters is supported.
      \item The functions res2table/res2tex to export the results are greatly enhanced.
      \item cluster-robust covariance: the method to compute it has been improved (more robust to rounding problems).
	  \item Computing standard errors: up to four-way clustering is supported.
    }
  }
  \subsection{New design}{
  	\itemize{
  		\item The function design has changed.
		\item There is only one function called femlm that performs \dQuote{normal} estimations and \dQuote{non.linear} in parameters estimations when adding the argument NL.fml.
		\item clusters can now directly be inserted in the formula without using the arument \dQuote{cluster}.
  }
  }
  \subsection{Argument name change}{
  	\itemize{
  		\item Name change: the argument \dQuote{dummy} becomes \dQuote{cluster}.
  }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item Bug related to the handling of the upper and lower bounds is now corrected.
    }
  }
}

