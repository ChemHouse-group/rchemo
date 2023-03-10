\name{dkplsr}
\alias{dkplsr}
\alias{transform.Dkpls}
\alias{coef.Dkpls}
\alias{predict.Dkplsr}
\encoding{latin1}

\title{Direct KPLSR Models}

\description{
Direct kernel PLSR (DKPLSR) (Bennett & Embrechts 2003). The method builds kernel Gram matrices and then runs a usual PLSR algorithm on them. This is faster (but not equivalent) to the "true" NIPALS KPLSR algorithm such as described in Rosipal & Trejo (2001). 
}

\usage{

dkplsr(X, Y, weights = NULL, nlv, kern = "krbf", ...)

\method{transform}{Dkpls}(object, X, ..., nlv = NULL)  

\method{coef}{Dkpls}(object, ..., nlv = NULL)  

\method{predict}{Dkplsr}(object, X, ..., nlv = NULL)  

}

\arguments{

For \code{dkplsr}: 

\item{X}{Matrix with the training X-data (\eqn{n, p}).}

\item{Y}{Matrix with the training Y-data (\eqn{n, q}).}

\item{weights}{vector of weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{nlv}{The number(s) of LVs to calculate.}

\item{kern}{Name of the function defining the considered kernel for building the Gram matrix. See \code{\link{krbf}} for syntax, and other available kernel functions (\code{krbf}, \code{kpol}, \code{ktanh}).}

\item{...}{Optional arguments to pass in the kernel function defined in \code{kern} (e.g. \code{gamma} for \code{\link{krbf}}, \code{gamma} and \code{coef0} for \code{ktanh}, \code{gamma} and \code{coef0} and \code{degree} for \code{kpol}).}



For \code{transform.dkplsr} and \code{predict.dkplsr}: 

\item{object}{A fitted model, output of a call to the main function.}

\item{X}{A matrix with new X-data (\eqn{m, p}) to consider.}

\item{...}{Optional arguments to pass to methods.}

\item{nlv}{The number(s) of LVs to consider.}


For \code{coef.dkplsr}: 

\item{object}{A fitted model, output of a call to the main function.}

\item{...}{Optional arguments to pass to methods.}

\item{nlv}{The number(s) of LVs to consider.}

}

\value{

For \code{dkplsr}: 

\item{X}{Matrix with the training X-data (\eqn{n, p}).}

\item{fm}{List with the outputs of the PLSR ((\code{T}): the X-score matrix (n,nlv); (\code{P}): the X-loadings matrix (p,nlv); (\code{R}): The PLS projection matrix (p,nlv); (\code{W}): The X-loading weights matrix (p,nlv); (\code{C}): The Y-loading weights matrix; (\code{TT}): the X-score normalization factor, (\code{xmeans}): the centering vector of X (p,1); (\code{ymeans}): the centering vector of Y (q,1); (\code{weights}): the weights vector of X-variables (p,1)}

\item{K}{kernel Gram matrix} 

\item{kern}{kernel function} 

\item{dots}{Optional arguments passed in the kernel function}

For \code{transform.dkplsr} : A matrix (\eqn{m, nlv}) with the projection of the new X-data on the X-scores

For \code{predict.dkplsr}: 

\item{pred}{A list of matrices (\eqn{m, q}) with the Y predicted values for the new X-data}

\item{K}{kernel Gram matrix (\eqn{m, nlv}), with values for the new X-data}

For \code{coef.dkplsr}: 

\item{int}{matrix (1,nlv) with the intercepts}

\item{B}{matrix (n,nlv) with the coefficients}

}


\references{

Bennett, K.P., Embrechts, M.J., 2003. An optimization perspective on kernel partial least squares regression, in: Advances in Learning Theory: Methods, Models and Applications, NATO Science Series III: Computer & Systems Sciences. IOS Press Amsterdam, pp. 227-250.

Rosipal, R., Trejo, L.J., 2001. Kernel Partial Least Squares Regression in Reproducing Kernel Hilbert Space. Journal of Machine Learning Research 2, 97-123.

}

\note{
The second example concerns the fitting of the function sinc(x) described in Rosipal & Trejo 2001 p. 105-106
}

\examples{

n <- 6 ; p <- 4
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- rnorm(n)
Ytrain <- cbind(y1 = ytrain, y2 = 100 * ytrain)
m <- 3
Xtest <- Xtrain[1:m, , drop = FALSE] 
Ytest <- Ytrain[1:m, , drop = FALSE] ; ytest <- Ytest[1:m, 1]

nlv <- 2
fm <- dkplsr(Xtrain, Ytrain, nlv = nlv, kern = "krbf", gamma = .8)
transform(fm, Xtest)
transform(fm, Xtest, nlv = 1)
coef(fm)
coef(fm, nlv = 1)

predict(fm, Xtest)
predict(fm, Xtest, nlv = 0:nlv)$pred

pred <- predict(fm, Xtest)$pred
msep(pred, Ytest)

nlv <- 2
fm <- dkplsr(Xtrain, Ytrain, nlv = nlv, kern = "kpol", degree = 2, coef0 = 10)
predict(fm, Xtest, nlv = nlv)

}

\examples{

x <- seq(-10, 10, by = .2)
x[x == 0] <- 1e-5
n <- length(x)
zy <- sin(abs(x)) / abs(x)
y <- zy + rnorm(n, 0, .2)
plot(x, y, type = "p")
lines(x, zy, lty = 2)
X <- matrix(x, ncol = 1)

nlv <- 3
fm <- dkplsr(X, y, nlv = nlv)
pred <- predict(fm, X)$pred
plot(X, y, type = "p")
lines(X, zy, lty = 2)
lines(X, pred, col = "red")

}

\keyword{datagen}