\name{FENmlm-package}
\alias{FENmlm-package}
\alias{FENmlm}
\docType{package}
\title{
Fixed Effects Nonlinear Maximum Likelihood Models
}
\description{
Efficient estimation of multiple fixed-effects maximum likelihood models
    with, possibly, non-linear in parameters right hand sides. Standard-errors
    can easily be clustered. It also includes tools to seamlessly export (to Latex)
    the results of various estimations.
}
\details{
\tabular{ll}{
Package: \tab FENmlm\cr
Type: \tab Package\cr
Version: \tab 2.1.0\cr
Date: \tab 2018-02-12\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
This package intends to efficiently estimate fixed-effect maximum likelihood models. The function \code{\link[FENmlm]{femlm}} estimates fixed-effect unconditionnal maximum likelihood models with, possibly, non-linear in parameters right hand sides. The ML families currently available are: poisson, negative binomial, logit and Gaussian.

Several features are also included such as the possibility to easily compute different types of standard-errors (including multi-way clustering). It is possible to compare the results of severeal estimations by using the function \code{\link[FENmlm]{res2table}}, and to export them to Latex using \code{\link[FENmlm]{res2tex}}.
}
\author{
Laurent Berge

Maintainer: Laurent Berge <laurent.berge at uni.lu>
}

\keyword{ package }
