\name{FENmlm-package}
\alias{FENmlm-package}
\alias{FENmlm}
\docType{package}
\title{
Fixed Effects Nonlinear Maximum Likelihood Models
}
\description{
Efficient estimation of fixed-effect maximum likelihood models with, possibly, non-linear right hand sides.
}
\details{
\tabular{ll}{
Package: \tab FENmlm\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-10-02\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
This package intends to efficiently estimate fixed-effect maximum likelihood models. The function \var{feNmlm} perfoms estimates fixed-effect maximum likelihood models with non-linear right hand sides. The function \var{femlm} is similar but is restricted to linear right hand sides.
Several features are also included such as the possibility to easily compute different types of standard-errors (including one-way and two-way clustering). It is possible to compare the results of severeal estimations by using the function \var{res2table}, and to export them to Latex using \var{res2tex}.
}
\author{
Laurent Berge

Maintainer: Laurent Berge <laurent.berge at u-bordeaux.fr>
}

\keyword{ package }
