\name{knnda}
\alias{knnda}
\alias{predict.Knnda}

\encoding{latin1}

\title{KNN-DA}

\description{

KNN weighted discrimination. For each new observation to predict, a number of \eqn{k} nearest neighbors is selected and the prediction is calculated by the most frequent class in \eqn{y} in this neighborhood.

}

\usage{

knnda(X, y,
    nlvdis, diss = c("eucl", "mahal"),
    h, k)

\method{predict}{Knnda}(object, X, ...)  

}

\arguments{

\item{X}{For the main function: Training X-data (\eqn{n, p}). --- For the auxiliary functions: New X-data (\eqn{m, p}) to consider.}

\item{y}{Training class membership (\eqn{n}). \bold{Note:} If \code{y} is a factor, it is replaced by a character vector.}

\item{nlvdis}{The number of LVs to consider in the global PLS used for the dimension reduction before calculating the dissimilarities. If \code{nlvdis = 0}, there is no dimension reduction. (see details)}

\item{diss}{The type of dissimilarity used for defining the neighbors. Possible values are "eucl" (default; Euclidean distance), "mahal" (Mahalanobis distance), or "correlation". Correlation dissimilarities are calculated by sqrt(.5 * (1 - rho)).}

\item{h}{A scale scalar defining the shape of the weight function. Lower is \eqn{h}, sharper is the function. See \code{\link{wdist}}.}

\item{k}{The number of nearest neighbors to select for each observation to predict.}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main function.}

\item{...}{For the auxiliary functions: Optional arguments. Not used.}

}

\details{

In function \code{knnda}, the dissimilarities used for computing the neighborhood and the weights can be calculated from the original X-data or after a dimension reduction (argument \code{nlvdis}). In the last case, global PLS scores are computed from \eqn{(X, Y)} and the dissimilarities are calculated on these scores. For high dimension X-data, the dimension reduction is in general required for using the Mahalanobis distance.   

}

\value{

For \code{knnda}:list with input arguments.

For \code{predict.Knnda}:

\item{pred}{prediction calculated for each observation by the most frequent class in \eqn{y} in its neighborhood.}

\item{listnn}{list with the neighbors used for each observation to be predicted}

\item{listd}{list with the distances to the neighbors used for each observation to be predicted}

\item{listw}{list with the weights attributed to the neighbors used for each observation to be predicted}

}

\references{
Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer.
}

\examples{

n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- sample(c(1, 4, 10), size = n, replace = TRUE)

m <- 5
Xtest <- Xtrain[1:m, ] ; ytest <- ytrain[1:m]

nlvdis <- 5 ; diss <- "mahal"
h <- 2 ; k <- 10
fm <- knnda(
    Xtrain, ytrain, 
    nlvdis = nlvdis, diss = diss,
    h = h, k = k
    )
res <- predict(fm, Xtest)
names(res)
res$pred
err(res$pred, ytest)

}

