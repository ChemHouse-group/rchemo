\name{orthog}
\alias{orthog}
\encoding{latin1}

\title{Orthogonalization of a matrix to another matrix}

\description{

Function \code{orthog} orthogonalizes a matrix \eqn{Y} to a matrix \eqn{X}. The row observations can be weighted. 

The function uses function \code{\link{lm}}.

}

\usage{

orthog(X, Y, weights = NULL)

}

\arguments{

\item{X}{A \eqn{n x p} matrix or data frame.}

\item{Y}{A \eqn{n x q} matrix or data frame to orthogonalize to \eqn{X}.}

\item{weights}{A vector of length \eqn{n} defining a priori weights to apply to the observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

}

\value{

\item{Y}{The \eqn{Y} matrix orthogonalized to \eqn{X}.}

\item{b}{The regression coefficients used for orthogonalization.}

}

\examples{

n <- 8 ; p <- 3
set.seed(1)
X <- matrix(rnorm(n * p, mean = 10), ncol = p, byrow = TRUE)
Y <- matrix(rnorm(n * 2, mean = 10), ncol = 2, byrow = TRUE)
colnames(Y) <- c("y1", "y2")
set.seed(NULL)
X
Y

res <- orthog(X, Y)
res$Y
crossprod(res$Y, X)
res$b

# Same as:
fm <- lm(Y ~ X)
Y - fm$fitted.values
fm$coef

#### WITH WEIGHTS

w <- 1:n
fm <- lm(Y ~ X, weights = w)
Y - fm$fitted.values
fm$coef

res <- orthog(X, Y, weights = w)
res$Y
t(res$Y) %*% diag(w) %*% X
res$b

}

\keyword{datagen}