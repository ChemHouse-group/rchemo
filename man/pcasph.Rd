\name{pcasph}
\alias{pcasph}
\encoding{latin1}

\title{Spherical PCA}

\description{



}

\usage{

pcasph(X, weights = NULL, nlv)

}

\arguments{

\item{X}{For the main functions and auxiliary function \code{summary}: Training X-data (\eqn{n, p}). --- For the other auxiliary functions: New X-data (\eqn{m, p}) to consider.}

\item{weights}{Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{nlv}{The number of PCs to calculate.}

\item{object}{A fitted model, output of a call to the main functions.}

\item{...}{Optional arguments.}

\bold{Specific for the NIPALS algorithm}

\item{gs}{Logical indicating if a Gram-Schmidt orthogonalization is implemented or not (default to \code{TRUE}).}

\item{tol}{Tolerance for testing convergence of the NIPALS iterations for each PC.}

\item{maxit}{Maximum number of NIPALS iterations for each PC.}


}

\value{

A list of outputs, such as:

\item{T}{The score matrix (\eqn{n, nlv}).}

\item{P}{The loadings matrix (\eqn{p, nlv}).}

\item{R}{The projection matrix (= \eqn{P}) (\eqn{p, nlv}).}

\item{sv}{The singular values (\eqn{min(n, p), 1}) except for NIPALS = (\eqn{nlv, 1}).}

\item{eig}{The eigenvalues (\code{= sv^2}) (\eqn{min(n, p), 1}) except for NIPALS = (\eqn{nlv, 1}).}

\item{xmeans}{The centering vector of \eqn{X} (\eqn{p, 1}).}

\item{niter}{Numbers of iterations of the NIPALS.}

\item{conv}{Logical indicating if the NIPALS converged before reaching the maximal number of iterations.}

}

\references{

Daszykowski, M., Kaczmarek, K., Vander Heyden, Y., Walczak, B., 2007. 
Robust statistics in data analysis - A review. Chemometrics and Intelligent 
Laboratory Systems 85, 203-219. https://doi.org/10.1016/j.chemolab.2006.06.016

Locantore N., Marron J.S., Simpson D.G., Tripoli N., Zhang J.T., Cohen K.L.
Robust principal component analysis for functional data, Test 8 (1999) 1–7

Maronna, R., 2005. Principal components and orthogonal regression based on 
robust scales, Technometrics, 47:3, 264-273, DOI: 10.1198/004017005000000166

}

\examples{

n <- 6 ; p <- 4
Xtrain <- matrix(rnorm(n * p), nrow = n)
s <- c(3, 4, 7, 10, 11, 15, 21:24)   
zX <- replace(Xtrain, s, NA)
Xtrain
zX
m <- 2
Xtest <- matrix(rnorm(m * p), nrow = m)

pcasvd(Xtrain, nlv = 3)
pcaeigen(Xtrain, nlv = 3)
pcaeigenk(Xtrain, nlv = 3)
pcanipals(Xtrain, nlv = 3)
pcanipalsna(Xtrain, nlv = 3)
pcanipalsna(zX, nlv = 3)

fm <- pcaeigen(Xtrain, nlv = 3)
fm$T
transform(fm, Xtest)
transform(fm, Xtest, nlv = 2)

pcaeigen(Xtrain, nlv = 3)$T
pcaeigen(Xtrain, nlv = 3, weights = 1:n)$T

## Scores

Ttrain <- fm$T
Ttest <- transform(fm, Xtest)
T <- rbind(Ttrain, Ttest)
group <- c(rep("Training", nrow(Ttrain)), rep("Test", nrow(Ttest)))
i <- 1
plotxy(T[, i:(i+1)], group = group, pch = 16, zeroes = TRUE, cex = 1.3)

## Loadings

plotxy(fm$P, zeroes = TRUE, label = TRUE, cex = 2, col = "red3")

## Summary

summary(fm, Xtrain)

res <- summary(fm, Xtrain)
plotxy(res$cor.circle, zeroes = TRUE, label = TRUE, cex = 2, col = "red3",
    circle = TRUE, ylim = c(-1, 1))

}

\keyword{datagen}