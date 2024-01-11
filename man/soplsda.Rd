\name{soplsrda}
\alias{soplsrda}
\alias{soplslda}
\alias{soplsqda}
\alias{transform.Soplsrda} 
\alias{transform.Soplsprobda}
\alias{predict.Soplsrda}
\alias{predict}.Soplsprobda} 

\encoding{UTF-8}

\title{Block dimension reduction by SO-PLS-DA}

\description{

Function \code{soplsrda} implements dimension reductions of pre-selected blocks of variables (= set of columns) of a reference (= training) matrix, by sequential orthogonalization-PLS (said "SO-PLS") in a context of discrimination. 

The block reduction consists in calculating latent variables (= scores) for each block, each block being sequentially orthogonalized to the information computed from the previous blocks.

The function allows giving a priori weights to the rows of the reference matrix in the calculations.

In\code{soplslda} and \code{soplsqda}, probabilistic LDA and QDA are run over the PLS2 LVs, respectively.

}

\usage{

soplsrda(Xlist, y, scaling = c("none", "pareto", "sd")[1], weights = NULL, nlv, ...)

soplslda(Xlist, y, scaling = c("none", "pareto", "sd")[1], weights = NULL, nlv, ...)

soplsqda(Xlist, y, scaling = c("none", "pareto", "sd")[1], weights = NULL, nlv, ...)


\method{transform}{Soplsrda}(object, Xlist) 

\method{transform}{Soplsprobda}(object, Xlist) 

\method{predict}{Soplsrda}(object, Xlist) 

\method{predict}{Soplsprobda}(object, Xlist) 

}

\arguments{

\item{X_list}{A list of matrices or data frames of reference (= training) observations.}

\item{y}{Training class membership (\eqn{n}). \bold{Note:} If \code{y} is a factor, it is replaced by a character vector.}

\item{scaling}{vector (of length Xlist) of variable scaling for each datablock, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "ctreduced", uncorrected standard deviation is used.}

\item{weights}{a priori weights to the rows of the reference matrix in the calculations.}

\item{nlv}{A vector of same length as the number of blocks defining the number of scores to calculate for each block, or a single number. In this last case, the same number of scores is used for all the blocks.}

\item{prior}{The prior probabilities of the classes. Possible values are "unif" (default; probabilities are set equal for all the classes) or "prop" (probabilities are set equal to the observed proportions of the classes in \code{y}).}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main functions.}

\item{...}{Other arguments to pass in functions \code{\link{pls}}.}

}

\value{

For \code{soplsrda}, \code{soplslda}, \code{soplsqda}:

\item{fm}{list with the PLS models: (\code{T}): X-scores matrix; (\code{P}): X-loading matrix;(\code{R}): The PLS projection matrix (p,nlv); (\code{W}): X-loading weights matrix ;(\code{C}): The Y-loading weights matrix; (\code{TT}): the X-score normalization factor; (\code{xmeans}): the centering vector of X (p,1);  (\code{ymeans}): the centering vector of Y (q,1); (\code{weights}): vector of observation weights; (\code){scaling}: scaling; (\code){Xnorms}: "norm" of each block; (\code{U}): intermediate output.}

\item{lev}{classes}

\item{ni}{number of observations in each class}

For \code{transform.Soplsrda}, \code{transform.Soplsprobda}:  the LVs Calculated for the new matrices list \eqn{Xlist} from the model.

For \code{predict.Soplsrda}, \code{predict.Soplsprobda}: 

\item{pred}{predicted class for each observation}

\item{posterior}{calculated probability of belonging to a class for each observation}

}

\references{

- Biancolillo et al. , 2015. Combining SO-PLS and linear discriminant analysis for
multi-block classification. Chemometrics and Intelligent Laboratory Systems, 141, 58-67.

- Biancolillo, A. 2016. Method development in the area of multi-block analysis focused on food analysis. PhD. University of copenhagen.

- Menichelli et al., 2014. SO-PLS as an exploratory tool for path modelling. Food Quality and Preference, 36, 122-134.

- Tenenhaus, M., 1998. La régression PLS: théorie et pratique. Editions Technip, Paris, France.

}


\examples{

N <- 10 ; p <- 12
set.seed(1)
X <- matrix(rnorm(N * p, mean = 10), ncol = p, byrow = TRUE)
y <- matrix(sample(c("1", "4", "10"), size = N, replace = TRUE), ncol=1)
colnames(X) <- paste("x", 1:ncol(X), sep = "")
set.seed(NULL)

n <- nrow(X)

X_list <- list(X[,1:4], X[,5:7], X[,9:ncol(X)])
X_list_2 <- list(X[1:2,1:4], X[1:2,5:7], X[1:2,9:ncol(X)])

# EXEMPLE WITH SO-PLS-RDA
ncomp <- 2
fm <- soplsrda(X_list, y, nlv = ncomp)
predict(fm,X_list_2)
transform(fm,X_list_2)

ncomp <- c(2, 0, 3)
fm <- soplsrda(X_list, y, nlv = ncomp)
predict(fm,X_list_2)
transform(fm,X_list_2)

ncomp <- 0
fm <- soplsrda(X_list, y, nlv = ncomp)
predict(fm,X_list_2)
transform(fm,X_list_2)

# EXEMPLE WITH SO-PLS-LDA
ncomp <- 2
weights <- rep(1 / n, n)
#w <- 1:n
soplslda(X_list, y, scaling = "none", nlv = ncomp, weights = weights)
soplslda(X_list, y, scaling = "pareto", nlv = ncomp, weights = weights)
soplslda(X_list, y, scaling = "sd", nlv = ncomp, weights = weights)

fm <- soplslda(X_list, y, scaling = c("none","pareto","sd"), nlv = ncomp, weights = weights)
predict(fm,X_list_2)
transform(fm,X_list_2)

}

\keyword{datagen}