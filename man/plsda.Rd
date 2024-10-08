\name{plsrda}
\alias{plsrda}
\alias{plslda}
\alias{plsqda}
\alias{predict.Plsrda}
\alias{predict.Plsprobda}
\encoding{latin1}

\title{PLSDA models}

\description{

Discrimination (DA) based on PLS.

The training variable \eqn{y} (univariate class membership) is firstly transformed to a dummy table containing \eqn{nclas} columns, where \eqn{nclas} is the number of classes present in \eqn{y}. Each column is a dummy variable (0/1). Then, a PLS2 is implemented on the \eqn{X-}data and the dummy table, returning latent variables (LVs) that are used as dependent variables in a DA model.

- \code{plsrda}: Usual "PLSDA". A linear regression model predicts the Y-dummy table from the PLS2 LVs. This corresponds to the PLSR2 of the X-data and the Y-dummy table. For a given observation, the final prediction is the class corresponding to the dummy variable for which the prediction is the highest.

- \code{plslda} and \code{plsqda}: Probabilistic LDA and QDA are run over the PLS2 LVs, respectively.

}

\usage{

plsrda(X, y, weights = NULL, nlv, 
Xscaling = c("none","pareto","sd")[1], Yscaling = c("none","pareto","sd")[1])

plslda(X, y, weights = NULL, nlv, prior = c("unif", "prop"), 
Xscaling = c("none","pareto","sd")[1], Yscaling = c("none","pareto","sd")[1])

plsqda(X, y, weights = NULL, nlv, prior = c("unif", "prop"), 
Xscaling = c("none","pareto","sd")[1], Yscaling = c("none","pareto","sd")[1])

\method{predict}{Plsrda}(object, X, ..., nlv = NULL) 

\method{predict}{Plsprobda}(object, X, ..., nlv = NULL) 

}

\arguments{

\item{X}{For the main functions: Training X-data (\eqn{n, p}). --- For the auxiliary functions: New X-data (\eqn{m, p}) to consider.}

\item{y}{Training class membership (\eqn{n}). \bold{Note:} If \code{y} is a factor, it is replaced by a character vector.}

\item{weights}{Weights (\eqn{n}) to apply to the training observations for the PLS2. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{nlv}{The number(s) of LVs to calculate.}

\item{prior}{The prior probabilities of the classes. Possible values are "unif" (default; probabilities are set equal for all the classes) or "prop" (probabilities are set equal to the observed proportions of the classes in \code{y}).}

\item{Xscaling}{X variable scaling among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{Yscaling}{Y variable scaling, once converted to binary variables, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main functions.}

\item{...}{For the auxiliary functions: Optional arguments. Not used.}

}

\value{

For \code{plsrda}, \code{plslda}, \code{plsqda}:

\item{fm}{list with the model: (\code{T}): X-scores matrix; (\code{P}): X-loading matrix;(\code{R}): The PLS projection matrix (p,nlv); (\code{W}): X-loading weights matrix ;(\code{C}): The Y-loading weights matrix; (\code{TT}): the X-score normalization factor; (\code{xmeans}): the centering vector of X (p,1);  (\code{ymeans}): the centering vector of Y (q,1); (\code{xscales}): the scaling vector of X (p,1);  (\code{yscales}): the scaling vector of Y (q,1); (\code{weights}): vector of observation weights; (\code{U}): intermediate output.}

\item{lev}{classes}

\item{ni}{number of observations in each class}

For \code{predict.Plsrda}, \code{predict.Plsprobda}: 

\item{pred}{predicted class for each observation}

\item{posterior}{calculated probability of belonging to a class for each observation}

}

\seealso{

\code{\link{plsr_plsda_allsteps}} function to help determine the optimal number of latent variables, perform a permutation test, calculate model parameters and predict new observations.

}

\note{

The first example concerns PLSDA, and the second one concerns PLS LDA.
\code{fm} are PLS1 models, and \code{zfm} are PLS2 models to predict the disjunctive matrix.

}

\examples{

## EXAMPLE OF PLSDA

n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- sample(c(1, 4, 10), size = n, replace = TRUE)

Xtest <- Xtrain[1:5, ] ; ytest <- ytrain[1:5]

nlv <- 5
fm <- plsrda(Xtrain, ytrain, Xscaling = "sd", nlv = nlv)
names(fm)

predict(fm, Xtest)
predict(fm, Xtest, nlv = 0:2)$pred

pred <- predict(fm, Xtest)$pred
err(pred, ytest)

zfm <- fm$fm
transform(zfm, Xtest)
transform(zfm, Xtest, nlv = 1)
summary(zfm, Xtrain)
coef(zfm)
coef(zfm, nlv = 0)
coef(zfm, nlv = 2)

## EXAMPLE OF PLS LDA

n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- sample(c(1, 4, 10), size = n, replace = TRUE)
Xtest <- Xtrain[1:5, ] ; ytest <- ytrain[1:5]

nlv <- 5
fm <- plslda(Xtrain, ytrain, Xscaling = "sd", nlv = nlv)
predict(fm, Xtest)
predict(fm, Xtest, nlv = 1:2)$pred

zfm <- fm$fm[[1]]
class(zfm)
names(zfm)
summary(zfm, Xtrain)
transform(zfm, Xtest[1:2, ])
coef(zfm)

}

\keyword{datagen}