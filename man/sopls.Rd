\name{sopls}
\alias{sopls}

\encoding{latin1}

\title{Block dimension reduction by SO-PLS}

\description{

Function \code{blocksopls} implements dimension reductions of pre-selected blocks of variables (= set of columns) of a reference (= training) matrix, by sequential orthogonalization-PLS (said "SO-PLS"). 

SO-PLS is described for instance in Menichelli et al. (2014), Biancolillo et al. (2015) and Biancolillo (2016). 

The block reduction consists in calculating latent variables (= scores) for each block, each block being sequentially orthogonalized to the information computed from the previous blocks.

The function allows giving a priori weights to the rows of the reference matrix in the calculations.

}

\usage{

sopls(X_list, Y, nlv_vect, blocknames = NULL, weights = NULL, ...)

}

\arguments{

\item{X_list}{A list of matrices or data frames of reference (= training) observations.}

\item{Y}{A \eqn{n x q} matrix or data frame, or a vector of length \eqn{n}, of reference (= training) responses.}

\item{nlv_vect}{A vector of same length as the number of blocks defining the number of scores to calculate for each block, or a single number. In this last case, the same number of scores is used for all the blocks.}

\item{blocknames}{A list of same length as the number of blocks, with block names.}

\item{weights}{a priori weights to the rows of the reference matrix in the calculations.}

\item{...}{Other arguments to pass in functions \code{\link{pls}}.}

}

\value{

\item{W}{A list of matrices with the loadings weights calculated from the X-blocks.}

\item{T}{A list of matrices with the scores calculated from the X-blocks.}

\item{C}{A list of matrices with the Y-loading weights after each addition of X-block}

\item{BCoef}{A list of matrices with the coefficients calculated for the X-blocks.}

\item{VIP}{A list of matrices with the VIP values calculated for the variables of the X-blocks.}

\item{Ymeans}{A list of matrices with the mean Y values calculated before each addition of X-block.}

\item{Fit}{A matrice \eqn{n x q} with the calculated fitted values.}

}

\references{

- Biancolillo et al. , 2015. Combining SO-PLS and linear discriminant analysis for
multi-block classification. Chemometrics and Intelligent Laboratory Systems, 141, 58-67.

- Biancolillo, A. 2016. Method development in the area of multi-block analysis focused on food analysis. PhD. University of copenhagen.

- Menichelli et al., 2014. SO-PLS as an exploratory tool for path modelling. Food Quality and Preference, 36, 122-134.

}


\examples{

N <- 10 ; p <- 12
set.seed(1)
X <- matrix(rnorm(N * p, mean = 10), ncol = p, byrow = TRUE)
Y <- matrix(rnorm(N * 2, mean = 10), ncol = 2, byrow = TRUE)
colnames(X) <- paste("x", 1:ncol(X), sep = "")
colnames(Y) <- paste("y", 1:ncol(Y), sep = "")
set.seed(NULL)
X
Y

n <- nrow(X)

X_list <- list(X[,1:4], X[,5:7], X[,9:ncol(X)])

ncomp <- 2
sopls(X_list, Y, nlv_vect = ncomp)

ncomp <- c(2, 0, 3)
sopls(X_list, Y, nlv_vect = ncomp)

ncomp <- 0
sopls(X_list, Y, nlv_vect = ncomp)

ncomp <- 2
weights <- rep(1 / n, n)
#w <- 1:n
sopls(X_list, Y, nlv_vect = ncomp, weights = weights)

}

\keyword{datagen}