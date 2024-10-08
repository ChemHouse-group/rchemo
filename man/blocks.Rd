\name{blockscal}
\alias{blockscal}
\alias{mblocks}
\alias{hconcat}
\encoding{latin1}

\title{Block autoscaling}

\description{

Functions managing blocks of data.

- \code{blockscal}: Autoscales a list of blocks (i.e. sets of columns) of a training X-data, and eventually the blocks of new X-data. The scaling factor (computed on the training) is the "norm" of the block, i.e. the square root of the sum of the variances of each column of the block. 

- \code{mblocks}: Makes a list of blocks from X-data.

- \code{hconcat}: Concatenates horizontally the blocks of a list. 

}

\usage{

blockscal(Xtrain, X = NULL, weights = NULL)

mblocks(X, blocks)

hconcat(X)

}

\arguments{

\item{Xtrain}{A list of blocks of training X-data}

\item{X}{For \code{blockscal}: A list of blocks of new X-data. For \code{mblocks}: X-data. For \code{hconcat}: a list of blocks of X-data.}

\item{blocks}{A list (of same length as the number of blocks) giving the column numbers in \code{X}.}

\item{weights}{Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

}

\value{

For \code{mblocks}: a list of blocks of X-data.

For \code{hconcat}: a matrix concatenating a list of data blocks.

For \code{blockscal}: 

\item{Xtrain}{A list of blocks of training X-data, after block autoscaling.}

\item{X}{A list of blocks of new X-data, after block autoscaling.}  

\item{disp}{The scaling factor (computed on the training).}

}

\note{
The second example is equivalent to MB-PLSR
}

\examples{

n <- 10 ; p <- 10
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- rnorm(n)
m <- 2

Xtest <- matrix(rnorm(m * p), ncol = p)
colnames(Xtest) <- paste("v", 1:p, sep = "")
Xtrain
Xtest

blocks <- list(1:2, 4, 6:8)
zXtrain <- mblocks(Xtrain, blocks = blocks)
zXtest <- mblocks(Xtest, blocks = blocks)

zXtrain
blockscal(zXtrain, zXtest)

res <- blockscal(zXtrain, zXtest)
hconcat(res$Xtrain)
hconcat(res$X)

## example of equivalence with MB-PLSR

n <- 10 ; p <- 10
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- rnorm(n)
m <- 2

Xtest <- matrix(rnorm(m * p), ncol = p)
colnames(Xtest) <- paste("v", 1:p, sep = "")
Xtrain
Xtest

blocks <- list(1:2, 4, 6:8)
X1 <- mblocks(Xtrain, blocks = blocks)
X1 <- lapply(1:length(X1), function(x) scale(X1[[x]]))
res <- blockscal(X1)
zXtrain <- hconcat(res$Xtrain)

nlv <- 3
fm <- plskern(zXtrain, ytrain, nlv = nlv)

}

\keyword{datagen}