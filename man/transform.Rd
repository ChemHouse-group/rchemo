\name{transform}
\alias{transform}
\encoding{latin1}

\title{Generic transform function}

\description{Transformation of the X-data by a fitted model.}

\usage{

transform(object, X, ...)

}

\arguments{

\item{object}{A fitted model, output of a call to a fitting function.}

\item{X}{New X-data to consider.}

\item{...}{Optional arguments.}

}

\value{

the transformed X-data

}

\examples{

n <- 6 ; p <- 4
X <- matrix(rnorm(n * p), ncol = p)
y <- rnorm(n)

fm <- pcaeigen(X, nlv = 3)

fm$T
transform(fm, X[1:2, ], nlv = 2)

}

\examples{

n <- 6 ; p <- 4
X <- matrix(rnorm(n * p), ncol = p)
y <- rnorm(n)

fm <- plskern(X, y, nlv = 3)
fm$T
transform(fm, X[1:2, ], nlv = 2)

}

\keyword{datagen}