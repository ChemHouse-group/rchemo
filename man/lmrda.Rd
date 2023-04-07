\name{lmrda}
\alias{lmrda}
\alias{predict.Lmrda}
\encoding{latin1}

\title{LMR-DA models}

\description{

Discrimination (DA) based on linear regression (LMR).

}

\usage{

lmrda(X, y, weights = NULL)

\method{predict}{Lmrda}(object, X, ...) 

}

\arguments{

\item{X}{For the main function: Training X-data (\eqn{n, p}). --- For the auxiliary function: New X-data (\eqn{m, p}) to consider.}

\item{y}{Training class membership (\eqn{n}). \bold{Note:} If \code{y} is a factor, it is replaced by a character vector.}

\item{weights}{Weights (\eqn{n}) to apply to the training observations for the PLS2. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{object}{For the auxiliary function: A fitted model, output of a call to the main functions.}

\item{...}{For the auxiliary function: Optional arguments. Not used.}

}

\details{

The training variable \eqn{y} (univariate class membership) is transformed to a dummy table containing \eqn{nclas} columns, where \eqn{nclas} is the number of classes present in \eqn{y}. Each column is a dummy variable (0/1). Then, a linear regression model (LMR) is run on the \eqn{X-}data and the dummy table, returning predictions of the dummy variables. For a given observation, the final prediction is the class corresponding to the dummy variable for which the prediction is the highest.

}

\value{

For \code{lrmda}:

\item{fm}{List with the outputs((\code{coefficients}): coefficient matrix; (\code{residuals}): residual matrix; (\code{fitted.values}): the fitted mean values; (\code{effects}): component relating to the linear fit, for use by extractor functions; (\code{weights}): Weights (\eqn{n}) applied to the training observations for the PLS2; (\code{rank}): the numeric rank of the fitted linear model; (\code{assign}): component relating to the linear fit, for use by extractor functions; (\code{qr}): component relating to the linear fit, for use by extractor functions; (\code{df.residual}): the residual degrees of freedom;  (\code{xlevels}): (only where relevant) a record of the levels of the factors used in fitting; (\code{call}): the matched call; (\code{terms}): the terms object used; (\code{model}): the model frame used).}

\item{lev}{y levels.}

\item{ni}{number of observations by level of y.}

For \code{predict.Lrmda}:

\item{pred}{predicted classes of observations.}

\item{posterior}{posterior probability of belonging to a class for each observation.}

}

\examples{

n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- sample(c(1, 4, 10), size = n, replace = TRUE)
m <- 5
Xtest <- Xtrain[1:m, ] ; ytest <- ytrain[1:m]

fm <- lmrda(Xtrain, ytrain)
names(fm)
predict(fm, Xtest)

coef(fm$fm)

pred <- predict(fm, Xtest)$pred
err(pred, ytest)

}

\keyword{datagen}