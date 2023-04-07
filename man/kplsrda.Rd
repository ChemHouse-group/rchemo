\name{kplsrda}
\alias{kplsrda}
\alias{predict.Kplsrda}
\encoding{latin1}

\title{KPLSR-DA models}

\description{

Discrimination (DA) based on kernel PLSR (KPLSR)

}

\usage{

kplsrda(X, y, weights = NULL, nlv, kern = "krbf", ...)

\method{predict}{Kplsrda}(object, X, ..., nlv = NULL) 

}

\arguments{

\item{X}{For main function: Training X-data (\eqn{n, p}). --- For auxiliary function: New X-data (\eqn{m, p}) to consider.}

\item{y}{Training class membership (\eqn{n}). \bold{Note:} If \code{y} is a factor, it is replaced by a character vector.}

\item{weights}{Weights (\eqn{n}) to apply to the training observations for the PLS2. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{nlv}{For main function: The number(s) of LVs to calculate. --- For auxiliary function: The number(s) of LVs to consider.}

\item{kern}{Name of the function defining the considered kernel for building the Gram matrix. See \code{\link{krbf}} for syntax, and other available kernel functions.}

\item{...}{Optional arguments to pass in the kernel function defined in \code{kern} (e.g. \code{gamma} for \code{\link{krbf}}).}

\item{object}{For auxiliary function: A fitted model, output of a call to the main functions.}

}

\details{

The training variable \eqn{y} (univariate class membership) is transformed to a dummy table containing \eqn{nclas} columns, where \eqn{nclas} is the number of classes present in \eqn{y}. Each column is a dummy variable (0/1). Then, a kernel PLSR (KPLSR) is run on the \eqn{X-}data and the dummy table, returning predictions of the dummy variables. For a given observation, the final prediction is the class corresponding to the dummy variable for which the prediction is the highest.

}

\value{

For \code{kplsrda}: 

\item{fm}{list with the kplsrda model: (\code{X}): the training X-data (\eqn{n, p}); (\code{Kt}): the Gram matrix; (\code{T}): X-scores matrix; (\code{C}): The Y-loading weights matrix; (\code{U}): intermediate output; (\code{R}): The PLS projection matrix (p,nlv); (\code{ymeans}): the centering vector of Y (q,1); (\code{weights}): vector of observation weights; (\code{kern}): kern function; (\code{dots}): Optional arguments.}

\item{lev}{y levels}

\item{ni}{number of observations by level of y}

For \code{predict.Kplsrda}: 

\item{pred}{predicted class for each observation}

\item{posterior}{calculated probability of belonging to a class for each observation}

}

\examples{

n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- sample(c(1, 4, 10), size = n, replace = TRUE)
m <- 5
Xtest <- Xtrain[1:m, ] ; ytest <- ytrain[1:m]

nlv <- 2
fm <- kplsrda(Xtrain, ytrain, nlv = nlv)
names(fm)
predict(fm, Xtest)

pred <- predict(fm, Xtest)$pred
err(pred, ytest)

predict(fm, Xtest, nlv = 0:nlv)$posterior
predict(fm, Xtest, nlv = 0)$posterior

predict(fm, Xtest, nlv = 0:nlv)$pred
predict(fm, Xtest, nlv = 0)$pred

}

\keyword{datagen}