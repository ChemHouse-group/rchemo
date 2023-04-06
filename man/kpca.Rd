\name{kpca}
\alias{kpca}
\alias{transform.Kpca}
\alias{summary.Kpca}
\encoding{UTF-8}

\title{KPCA}

\description{
Kernel PCA (Scholkopf et al. 1997, Scholkopf & Smola 2002, Tipping 2001) by SVD factorization of the weighted Gram matrix \eqn{D^(1/2) * Phi(X) * Phi(X)' * D^(1/2)}. \eqn{D} is a (\eqn{n, n}) diagonal matrix of weights for the observations (rows of \eqn{X}).
}

\usage{

kpca(X, weights = NULL, nlv, kern = "krbf", ...)

\method{transform}{Kpca}(object, X, ..., nlv = NULL)  

\method{summary}{Kpca}(object, ...)  

}

\arguments{

\item{X}{For the main function: Training X-data (\eqn{n, p}). --- For the auxiliary functions: New X-data (\eqn{m, p}) to consider.}

\item{weights}{Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{nlv}{The number of PCs to calculate.}

\item{kern}{Name of the function defining the considered kernel for building the Gram matrix. See \code{\link{krbf}} for syntax, and other available kernel functions.}

\item{...}{Optional arguments to pass in the kernel function defined in \code{kern} (e.g. \code{gamma} for \code{\link{krbf}}).}

\item{object}{--- For the auxiliary functions: A fitted model, output of a call to the main functions.}

}

\value{

For \code{kpca}:

\item{X}{Training X-data (\eqn{n, p}).}

\item{Kt}{Gram matrix}

\item{T}{X-scores matrix.}

\item{P}{X-loadings matrix.}

\item{sv}{vector of singular values}

\item{eig}{vector of eigenvalues.}

\item{weights}{vector of observation weights.}

\item{kern}{kern function.}

\item{dots}{Optional arguments.}

For \code{transform.Kpca}: X-scores matrix for new X-data.

For \code{summary.Kpca}: 

\item{explvar}{explained variance matrix.}

}

\references{

Scholkopf, B., Smola, A., Müller, K.-R., 1997. Kernel principal component analysis, in: Gerstner, W., Germond, A., Hasler, M., Nicoud, J.-D. (Eds.), Artificial Neural Networks — ICANN 97, Lecture Notes in Computer Science. Springer, Berlin, Heidelberg, pp. 583-588. https://doi.org/10.1007/BFb0020217

Scholkopf, B., Smola, A.J., 2002. Learning with kernels: support vector machines, regularization, optimization, and beyond, Adaptive computation and machine learning. MIT Press, Cambridge, Mass.

Tipping, M.E., 2001. Sparse kernel principal component analysis. Advances in neural information processing systems, MIT Press. http://papers.nips.cc/paper/1791-sparse-kernel-principal-component-analysis.pdf

}

\examples{

## EXAMPLE 1

n <- 5 ; p <- 4
X <- matrix(rnorm(n * p), ncol = p)

nlv <- 3
kpca(X, nlv = nlv, kern = "krbf")

fm <- kpca(X, nlv = nlv, kern = "krbf", gamma = .6)
fm$T
transform(fm, X[1:2, ])
transform(fm, X[1:2, ], nlv = 1)
summary(fm)

## EXAMPLE 2

n <- 5 ; p <- 4
X <- matrix(rnorm(n * p), ncol = p)
nlv <- 3
pcasvd(X, nlv = nlv)$T
kpca(X, nlv = nlv, kern = "kpol")$T

}

\keyword{datagen}