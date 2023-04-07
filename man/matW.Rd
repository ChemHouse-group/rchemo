\name{matW}
\alias{matW}
\alias{matB}
\encoding{latin1}

\title{Between and within covariance matrices}

\description{
Calculation of within (\code{matW}) and between (\code{matB}) covariance matrices for classes of observations.
}

\usage{

matW(X, y)

matB(X, y)

}

\arguments{

\item{X}{Data (\eqn{n, p}) on whch are calculated the covariances.}

\item{y}{Class membership (\eqn{n, 1}).}

}

\details{The denominator in the variance calculations is \eqn{n}.}

\value{

For (\code{matW}):

\item{W}{within covariance matrix.}

\item{Wi}{list of covariance matrices for each class.}

\item{lev}{classes}

\item{ni}{number of observations in each per class}

For (\code{matB}):

\item{B}{between covariance matrix.}

\item{ct}{matrix of class centers.}

\item{lev}{classes}

\item{ni}{number of observations in each per class}

}


\examples{

n <- 8 ; p <- 3
X <- matrix(rnorm(n * p), ncol = p)
y <- sample(1:2, size = n, replace = TRUE)
X
y

matW(X, y)

matB(X, y)

matW(X, y)$W + matB(X, y)$B
(n - 1) / n * cov(X)

}

\keyword{datagen}