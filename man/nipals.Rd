\name{nipals}
\alias{nipals}
\encoding{latin1}

\title{nipals}

\description{
Nipals to compute the first score and loading vectors of a matrix.
}

\usage{nipals(X, tol = 1e-6, maxit = 100)}

\arguments{

\item{X}{Datamatrix (\eqn{n, p}).}

\item{tol}{Tolerance for testing convergence of the NIPALS iterations for each PC.}

\item{maxit}{Maximum number of NIPALS iterations for each PC.}

}

\value{

\item{u}{left singular vector (u * sv = scores).}

\item{v}{right singular vector (loadings).}

\item{sv}{singular value.}

}

\references{

K.R. Gabriel, S. Zamir, Lower rank approximation of matrices by least squares with any choice of weights, Technometrics 21 (1979) 489–498.

}

\examples{

n <- 8 ; p <- 6
X <- matrix(rnorm(n * p, mean = 10), ncol = p, byrow = TRUE)
nipals(X)

}

\keyword{datagen}