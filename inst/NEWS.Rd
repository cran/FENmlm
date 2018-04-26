\name{NEWS}
\alias{NEWS}

\title{News for \R Package \pkg{FENmlm}}

\section{Changes in version 2.1.1}{
	\subsection{New features}{
		\itemize{
			\item [femlm] Added argument 'precision.cluster' so that the user can set the level of precision manually (concerns the algorithm to obtain the cluster coefficients).
			\item [summary.femlm] For one-way clustering, the argument 'cluster' can now be a vector.
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
  		\item The function designing has changed.
		\item There is only one function called femlm that performs "normal" estimations and "non.linear" in parameters estimations when adding the argument NL.fml.
		\item clusters can now directly be inserted in the formula without using the arument "cluster"
  }
  }
  \subsection{Argument name change}{
  	\itemize{
  		\item Name change: the argument "dummy" becomes "cluster".
  }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item Bug related to the handling of the upper and lower bounds is now corrected.
    }
  }
}

