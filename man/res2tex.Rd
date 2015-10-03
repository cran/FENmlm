\name{res2tex}
\alias{res2tex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Facility to export the results of multiple \code{feNmlm} estimations in a Latex table.
}
\description{
This function aggregates the results of multiple estimations and display them in the form of  one Latex table whose rownames are the variables and the columns contain the coefficients and standard-errors.
}
\usage{
res2tex(..., sd = c("standard", "white", "cluster", "twoway"),
        cluster, digits = 4, pseudo = TRUE,
        title, sdBelow = TRUE, drop, order, dict, file, append = TRUE,
        convergence = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Used to capture different \code{feNmlm} objects. Note that any other type of element is discarded.
}
  \item{sd}{
A character scalar. The standard-error is computed similarly for each \code{feNmlm} object. Which kind of standard error should be prompted: \dQuote{standard} (default), \dQuote{White}, \dQuote{cluster} or \dQuote{twoway}?
}
  \item{cluster}{
A list of vectors. Used only if \code{sd = "cluster"} or \code{sd="twoway"}. The vectors should give the cluster of each observation. Note that if the estimation was run using \code{dummy}, the standard error is automatically clustered along the cluster given in \code{feNmlm}.
}
  \item{digits}{
Integer. The number of digits to be displayed.
}
  \item{pseudo}{
Logical. Should the pseudo R2 be displayed? (Default is \code{TRUE}.)
}
  \item{title}{
Character scalar. The title of the Latex table.
}
  \item{sdBelow}{
Logical. Should the standard-error be displayed below the coefficients? (Default is \code{TRUE}.)
}
  \item{drop}{
Character vector. This element is used if some variables are not to be displayed. This should be a regular expression (see regexp help for more info). There can be more than one regular expression. Each variable satisfying the regular expression will be discarded.
}
  \item{order}{
Character vector. This element is used if the user wants the variables to be ordered in a certain way. This should be a regular expression (see regexp help for more info). There can be more than one regular expression. The variables satisfying the first regular expression will be placed first, then the order follows the sequence of regular expressions.
}
  \item{dict}{
A named character vector. Default is missing. If provided, it changes the original variable names to the ones contained in the \code{dict}. Example: I want to change my variable named "a" to "$log(a)$" and "b3" to "$bonus^3$", then I used \code{dict=c(a="$log(a)$",b3="$bonus^3$")}.
}
  \item{file}{
A character scalar. Default is missing. If provided, the Latex table will be saved in a file whose path is \code{file}.
}
  \item{append}{
Logical. Only used if option \code{file} is used. Should the Latex table be appended to the existing file? (Default is \code{TRUE}.)
}
  \item{convergence}{
Logical. Should the convergence state of the algorithm be displayed? (Default is \code{TRUE}.)
}
}

\value{
There is nothing returned, the result is only displayed on the console or saved in a file.
}

\author{
Laurent Berge
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

n = 100
x = rnorm(n,1,5)**2
y = rnorm(n,-1,5)**2
z = rpois(n,x*y)
base = data.frame(x,y,z)

# Results of the Poisson..
est_poisson = femlm(z~log(x)+log(y),base,family="poisson")
# .. and of the Negative Binomial
est_negbin = femlm(z~log(x)+log(y),base,family="negbin")

# We export the two results in one Latex table:
res2tex(est_poisson,est_negbin)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
