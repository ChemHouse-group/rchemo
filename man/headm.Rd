\name{headm}
\alias{headm}
\encoding{latin1}

\title{Display of the first part of a data set}

\description{

Function \code{headm} displays the first part and the dimension of a data set.

}

\usage{

headm(X)

}

\arguments{

\item{X}{A matrix or dataframe.}

}

\value{

first 6 rows and columns of a dataset, number of rows, number of columns, dataset class.

}



\examples{

n <- 1000
p <- 200
X <- matrix(rnorm(n * p), nrow = n)

headm(X)

}

\keyword{datagen}