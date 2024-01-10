\name{mbplsrda}
\alias{mbplsrda}
\alias{mbplslda}
\alias{mbplsqda}
\alias{predict.Mbplsrda}
\alias{predict.Mbplsprobda}
\encoding{latin1}

\title{multi-block PLSDA models}

\description{

Multi-block discrimination (DA) based on PLS.

The training variable \eqn{y} (univariate class membership) is firstly transformed to a dummy table containing \eqn{nclas} columns, where \eqn{nclas} is the number of classes present in \eqn{y}. Each column is a dummy variable (0/1). Then, a PLS2 is implemented on the \eqn{X-}data and the dummy table, returning latent variables (LVs) that are used as dependent variables in a DA model.

- \code{mbplsrda}: Usual "PLSDA". A linear regression model predicts the Y-dummy table from the PLS2 LVs. This corresponds to the PLSR2 of the X-data and the Y-dummy table. For a given observation, the final prediction is the class corresponding to the dummy variable for which the prediction is the highest.

- \code{mbplslda} and \code{mbplsqda}: Probabilistic LDA and QDA are run over the PLS2 LVs, respectively.

}

\usage{

mbplsrda(Xlist, y, scaling = c("centered", "pareto", "ctreduced")[1], blockscaling = TRUE, weights = NULL, nlv)

mbplslda(Xlist, y, scaling = c("centered", "pareto", "ctreduced")[1], blockscaling = TRUE, weights = NULL, nlv, prior = c("unif", "prop"))

mbplsqda(Xlist, y, scaling = c("centered", "pareto", "ctreduced")[1], blockscaling = TRUE, weights = NULL, nlv, prior = c("unif", "prop"))

\method{predict}{Mbplsrda}(object, Xlist, ..., nlv = NULL) 

\method{predict}{Mbplsprobda}(object, Xlist, ..., nlv = NULL) 

}

\arguments{

\item{Xlist}{For the main functions: list of training X-data (\eqn{n}rows). --- For the auxiliary functions: list of new X-data (\eqn{n} rows), with the same variables than the training X-data.}

\item{y}{Training class membership (\eqn{n}). \bold{Note:} If \code{y} is a factor, it is replaced by a character vector.}

\item{scaling}{vector (of length Xlist) of variable scaling for each datablock, among "centered", "pareto", "ctreduced" (for centered and scaled to unit variance). If "pareto" or "ctreduced", uncorrected standard deviation is used. Note: After binary transformation, Y-variables are always mean-centered.}

\item{blockscaling}{logical. If TRUE, the scaling factor (computed on the training) is the "norm" of the block, i.e. the square root of the sum of the variances of each column of the block.}

\item{weights}{Weights (\eqn{n}) to apply to the training observations for the PLS2. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{nlv}{The number(s) of LVs to calculate.}

\item{prior}{The prior probabilities of the classes. Possible values are "unif" (default; probabilities are set equal for all the classes) or "prop" (probabilities are set equal to the observed proportions of the classes in \code{y}).}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main functions.}

\item{...}{For the auxiliary functions: Optional arguments. Not used.}

}

\value{

For \code{mbplsrda}, \code{mbplslda}, \code{mbplsqda}:

\item{fm}{list with the model: (\code{T}): X-scores matrix; (\code{P}): X-loading matrix;(\code{R}): The PLS projection matrix (p,nlv); (\code{W}): X-loading weights matrix ;(\code{C}): The Y-loading weights matrix; (\code{TT}): the X-score normalization factor; (\code{xmeans}): the centering vector of X (p,1);  (\code{ymeans}): the centering vector of Y (q,1); (\code{weights}): vector of observation weights; (\code){scaling}: scaling; (\code){blockscaling}: block scaling; (\code){Xnorms}: "norm" of each block; (\code{U}): intermediate output.}

\item{lev}{classes}

\item{ni}{number of observations in each class}

For \code{predict.Mbplsrda}, \code{predict.Mbplsprobda}: 

\item{pred}{predicted class for each observation}

\item{posterior}{calculated probability of belonging to a class for each observation}

}

\note{

The first example concerns MB-PLSDA, and the second one concerns MB-PLS LDA.
\code{fm} are PLS1 models, and \code{zfm} are PLS2 models.

}

\examples{

## EXAMPLE OF MB-PLSDA

n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
Xtrainlist <- list(Xtrain[,1:3], Xtrain[,4:8])

ytrain <- sample(c(1, 4, 10), size = n, replace = TRUE)

Xtest <- Xtrain[1:5, ] ; ytest <- ytrain[1:5]
Xtestlist <- list(Xtest[,1:3], Xtest[,4:8])

nlv <- 5
fm <- mbplsrda(Xtrainlist, ytrain, scaling = "ctreduced", nlv = nlv)
names(fm)

predict(fm, Xtestlist)
predict(fm, Xtestlist, nlv = 0:2)$pred

pred <- predict(fm, Xtestlist)$pred
err(pred, ytest)

zfm <- fm$fm
transform(zfm, Xtestlist)
transform(zfm, Xtestlist, nlv = 1)
summary(zfm, Xtrainlist)
coef(zfm)
coef(zfm, nlv = 0)
coef(zfm, nlv = 2)

## EXAMPLE OF MB-PLS LDA

n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
Xtrainlist <- list(Xtrain[,1:3], Xtrain[,4:8])

ytrain <- sample(c(1, 4, 10), size = n, replace = TRUE)

Xtest <- Xtrain[1:5, ] ; ytest <- ytrain[1:5]
Xtestlist <- list(Xtest[,1:3], Xtest[,4:8])

nlv <- 5
fm <- mbplslda(Xtrainlist, ytrain, scaling = "centered", nlv = nlv)
predict(fm, Xtestlist)
predict(fm, Xtestlist, nlv = 1:2)$pred

zfm <- fm[[1]][[1]]
class(zfm)
names(zfm)
summary(zfm, Xtrainlist)
transform(zfm, Xtestlist)
coef(zfm)

## EXAMPLE OF MB-PLS QDA

n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
Xtrainlist <- list(Xtrain[,1:3], Xtrain[,4:8])

ytrain <- sample(c(1, 4, 10), size = n, replace = TRUE)

Xtest <- Xtrain[1:5, ] ; ytest <- ytrain[1:5]
Xtestlist <- list(Xtest[,1:3], Xtest[,4:8])

nlv <- 5
fm <- mbplsqda(Xtrainlist, ytrain, scaling = "centered", nlv = nlv)
predict(fm, Xtestlist)
predict(fm, Xtestlist, nlv = 1:2)$pred

zfm <- fm[[1]][[1]]
class(zfm)
names(zfm)
summary(zfm, Xtrainlist)
transform(zfm, Xtestlist)
coef(zfm)

}

\keyword{datagen}