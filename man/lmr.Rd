\name{lmr}
\alias{lmr}
\alias{coef.Lmr}
\alias{predict.Lmr}
\encoding{latin1}

\title{Linear regression models}

\description{
Linear regression models (uses function \code{\link{lm}}).
}

\usage{

lmr(X, Y, weights = NULL)

\method{coef}{Lmr}(object, ...) 

\method{predict}{Lmr}(object, X, ...)  

}

\arguments{

For \code{lmr}:

\item{X}{Training X-data (\eqn{n, p}).}

\item{Y}{Training Y-data (\eqn{n, q}).}

\item{weights}{Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

For \code{coef.Lmr} and \code{predict.Lmr}:

\item{object}{A fitted model, output of a call to the main functions.}

\item{X}{New X-data (\eqn{m, p}) to consider.}

\item{...}{Optional arguments. Not used.}

}

\value{

For \code{lmr}:

\item{coefficients}{coefficient matrix.}

\item{residuals}{residual matrix.}

\item{effects}{component relating to the linear fit, for use by extractor functions.}

\item{rank}{the numeric rank of the fitted linear model.}

\item{fitted.values}{the fitted mean values.}

\item{assign}{component relating to the linear fit, for use by extractor functions.}

\item{qr}{component relating to the linear fit, for use by extractor functions.}

\item{df.residual}{the residual degrees of freedom.}

\item{xlevels}{(only where relevant) a record of the levels of the factors used in fitting.}

\item{call}{the matched call.}

\item{terms}{the terms object used.}

\item{model}{the model frame used.}

For \code{coef.Lmr}: 

\item{int}{matrix (1,nlv) with the intercepts}

\item{B}{matrix (n,nlv) with the coefficients}

For \code{predict.Lmr}: 

\item{pred}{A list of matrices (\eqn{m, q}) with the Y predicted values for the new X-data}

}

\examples{

n <- 8 ; p <- 3
X <- matrix(rnorm(n * p, mean = 10), ncol = p, byrow = TRUE)
y <- rnorm(n)
Y <- cbind(y, rnorm(n))
Xtrain <- X[1:6, ] ; Ytrain <- Y[1:6, ]
Xtest <- X[7:8, ] ; Ytest <- Y[7:8, ]

fm <- lmr(Xtrain, Ytrain)
coef(fm)

predict(fm, Xtest)

pred <- predict(fm, Xtest)$pred
msep(pred, Ytest)

}

\keyword{datagen}