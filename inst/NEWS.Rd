\name{NEWS}
\alias{NEWS}

\title{News for \R Package \pkg{FENmlm}}

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

