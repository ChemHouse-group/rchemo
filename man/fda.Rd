\name{fda}
\alias{fda}
\alias{fdasvd}
\alias{transform.Fda}
\alias{summary.Fda}

\title{Factorial discriminant analysis}

\description{
Factorial discriminant analysis (FDA). The functions maximize the compromise \eqn{p'Bp / p'Wp}, i.e. \eqn{max p'Bp} with constraint \eqn{p'Wp = 1}. Vectors \eqn{p} are the linear discrimant coefficients "LD".

- \code{fda}: Eigen factorization of \eqn{W^(-1)B}

- \code{fdasvd}: Weighted SVD factorization of the matrix of the class centers.

If \eqn{W} is singular, W^(-1) is replaced by a MP pseudo-inverse.
}

\usage{

fda(X, y, nlv = NULL)

fdasvd(X, y, nlv = NULL)

\method{transform}{Fda}(object, X, ..., nlv = NULL) 

\method{summary}{Fda}(object, ...) 

}

\arguments{

For \code{fda} and \code{fdasvd}:

\item{X}{For the main functions: Training X-data (\eqn{n, p}).}

\item{y}{Training class membership (\eqn{n}). \bold{Note:} If \code{y} is a factor, it is replaced by a character vector.}

\item{nlv}{The number(s) of LVs to calculate.}

For \code{transform.Fda}

\item{object}{A fitted model, output of a call to the main function.}

\item{X}{New X-data (\eqn{m, p}) to consider.}

\item{...}{Optional arguments. Not used.}

\item{nlv}{The number(s) of LVs to consider.}

For \code{summary.Fda}

\item{object}{A fitted model, output of a call to the main function.}

\item{...}{Optional arguments. Not used.}

}

\value{

For \code{fda} and \code{fdasvd}:

\item{T}{X-scores matrix (n,nlv).}

\item{P}{X-loadings matrix (p,nlv) = coefficients of the linear discriminant function = "LD" of function lda of package MASS.}

\item{Tcenters}{projection of the class centers in the score space.}

\item{eig}{vector of eigen values}

\item{sstot}{total variance}

\item{W}{unbiased within covariance matrix}

\item{xmeans}{means of the X variables}

\item{lev}{y levels}

\item{ni}{number of observations per level of the y variable}

For \code{transform.Fda}: scores of the new X-data in the model.

For \code{summary.Fda}:

\item{explvar}{Explained variance by PCA of the class centers in transformed scale.}

}

\references{

Saporta G., 2011. Probabilités analyse des données et statistique. Editions Technip, Paris, France.

}

\examples{

data(iris)

X <- iris[, 1:4]
y <- iris[, 5]
table(y)

fm <- fda(X, y)
headm(fm$T)

transform(fm, X[1:3, ])

summary(fm)
plotxy(fm$T, group = y, ellipse = TRUE, 
    zeroes = TRUE, pch = 16, cex = 1.5, ncol = 2)
points(fm$Tcenters, pch = 8, col = "blue", cex = 1.5)

}

\keyword{datagen}