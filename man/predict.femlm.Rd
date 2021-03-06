% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MiscFuns.R
\name{predict.femlm}
\alias{predict.femlm}
\title{Predict method for femlm fits}
\usage{
\method{predict}{femlm}(object, newdata, type = c("response", "link"),
  ...)
}
\arguments{
\item{object}{An object of class \code{femlm}. Typically the result of a \code{\link[FENmlm]{femlm}} estimation.}

\item{newdata}{A data.frame containing the variables used to make the prediction. If not provided, the fitted expected (or linear if \code{type = "link"}) predictors are returned.}

\item{type}{Character either equal to \code{"response"} (default) or \code{"link"}. If \code{type="response"}, then the output is at the level of the response variable, i.e. it is the expected predictor \eqn{E(Y|X)}. If \code{"link"}, then the output is at the level of the explanatory variables, i.e. the linear predictor \eqn{X\cdot \beta}.}

\item{...}{Not currently used.}
}
\value{
It returns a numeric vector of length equal to the number of observations in argument \code{newdata}.
}
\description{
This function obtains prediction from a fitted model estimated with \code{\link[FENmlm]{femlm}}.
}
\examples{

# Estimation on iris data
res = femlm(Sepal.Length ~ Petal.Length | Species, iris)

# what would be the prediction if the data was all setosa?
newdata = data.frame(Petal.Length = iris$Petal.Length, Species = "setosa")
pred_setosa = predict(res, newdata = newdata)

# Let's look at it graphically
plot(c(1, 7), c(3, 11), type = "n", xlab = "Petal.Length",
     ylab = "Sepal.Length")

newdata = iris[order(iris$Petal.Length), ]
newdata$Species = "setosa"
lines(newdata$Petal.Length, predict(res, newdata))

# versicolor
newdata$Species = "versicolor"
lines(newdata$Petal.Length, predict(res, newdata), col=2)

# virginica
newdata$Species = "virginica"
lines(newdata$Petal.Length, predict(res, newdata), col=3)

# The original data
points(iris$Petal.Length, iris$Sepal.Length, col = iris$Species, pch = 18)
legend("topleft", lty = 1, col = 1:3, legend = levels(iris$Species))

}
\seealso{
\code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{update.femlm}}, \code{\link[FENmlm]{summary.femlm}}, \code{\link[FENmlm]{vcov.femlm}}, \code{\link[FENmlm]{getFE}}.
}
\author{
Laurent Berge
}
