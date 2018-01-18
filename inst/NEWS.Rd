\name{NEWS}
\alias{NEWS}

\title{News for \R Package \pkg{FENmlm}}

\section{Changes in version 1.1}{
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

