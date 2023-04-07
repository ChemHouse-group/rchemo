\name{segmkf}
\alias{segmkf}
\alias{segmts}
\encoding{latin1}

\title{Segments for cross-validation}

\description{

Build segments of observations for K-Fold or "test-set" cross-validation (CV). 

The CV can eventually be randomly repeated. For each repetition:

- \bold{K-fold CV} - Function \code{segmkf} returns the \eqn{K} segments.

- \bold{Test-set CV} - Function \code{segmts} returns a segment (of a given length) randomly sampled in the dataset.

\bold{CV of blocks} 

Argument \code{y} allows sampling \bold{blocks of observations} instead of observations. This can be required when there are repetitions in the data. In such a situation, CV should account for the repetition level (if not, the error rates are in general highly underestimated). For implementing such a CV, object \code{y} must be a a vector (\eqn{n}) defining the blocks, in the same order as in the data.

In any cases (\code{y = NULL} or not), the functions return a list of vector(s). Each vector contains the indexes of the observations defining the segment.

}

\usage{

segmkf(n, y = NULL, K = 5, 
    type = c("random", "consecutive", "interleaved"), nrep = 1) 

segmts(n, y = NULL, m, nrep) 

}

\arguments{

\item{n}{The total number of row observations in the dataset. If \code{y = NULL}, the CV is implemented on \code{1:n}. If \code{y != NULL}, blocks of observations (defined in \code{y}) are sampled instead of observations (but indexes of observations are returned).}

\item{y}{A vector (\eqn{n}) defining the blocks. Default to \code{NULL}.}

\item{K}{For \code{segmkf}.The number of folds (i.e. segments) in the K-fold CV.}

\item{type}{For \code{segmkf}.The type K-fold CV. Possible values are "random" (default), "consecutive" and "interleaved".}

\item{m}{For \code{segmts}. If \code{y = NULL}, the number of observations in the segment. If not, the number of blocks in the segment.}

\item{nrep}{The number of replications of the repeated CV. Default to \code{nrep = 1}.}

}

\value{

The segments (lists of indexes).

}


\examples{

Kfold <- segmkf(n = 10, K = 3)

interleavedKfold <- segmkf(n = 10, K = 3, type = "interleaved")

LeaveOneOut <- segmkf(n = 10, K = 10)

RepeatedKfold <- segmkf(n = 10, K = 3, nrep = 2)

repeatedTestSet <- segmts(n = 10, m = 3, nrep = 5)

n <- 10
y <- rep(LETTERS[1:5], 2)
y

Kfold_withBlocks <- segmkf(n = n, y = y, K = 3, nrep = 1)
z <- Kfold_withBlocks 
z
y[z$rep1$segm1]
y[z$rep1$segm2]
y[z$rep1$segm3]

TestSet_withBlocks <- segmts(n = n, y = y, m = 3, nrep = 1)
z <- TestSet_withBlocks
z
y[z$rep1$segm1]

}

\keyword{datagen}