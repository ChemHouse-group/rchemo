\name{dkrr}
\alias{dkrr}
\alias{coef.Dkrr}
\alias{predict.Dkrr}
\encoding{latin1}

\title{Direct KRR Models}

\description{

Direct kernel ridge regression (DKRR), following the same approcah as for DKPLSR (Bennett & Embrechts 2003). The method builds kernel Gram matrices and then runs a RR algorithm on them. This is not equivalent to the "true" KRR (= LS-SVM) algorithm. 

}

\usage{

dkrr(X, Y, weights = NULL, lb = 1e-2, kern = "krbf", ...)

\method{coef}{Dkrr}(object, ..., lb = NULL)  

\method{predict}{Dkrr}(object, X, ..., lb = NULL)  

}

\arguments{

For \code{dkrr}:

\item{X}{For the main function: Training X-data (\eqn{n, p}). --- For the auxiliary functions: New X-data (\eqn{m, p}) to consider.}

\item{Y}{Training Y-data (\eqn{n, q}).}

\item{weights}{Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{lb}{A value of regularization parameter \eqn{lambda}.}

\item{kern}{Name of the function defining the considered kernel for building the Gram matrix. See \code{\link{krbf}} for syntax, and other available kernel functions.}

\item{object}{A fitted model, output of a call to the main function.}

\item{...}{Optional arguments to pass in the kernel function defined in \code{kern} (e.g. \code{gamma} for \code{\link{krbf}}).}

For \code{coeff.Dkrr}:

\item{object}{A fitted model, output of a call to the main function.}

\item{lb}{A value of regularization parameter \eqn{lambda}.}

For \code{predict.Dkrr}:

\item{object}{A fitted model, output of a call to the main function.}

\item{X}{A matrix with new X-data (\eqn{m, p}) to consider.}

\item{lb}{A value of regularization parameter \eqn{lambda}.}

}

\value{

For \code{dkrr}: 

\item{X}{Matrix with the training X-data (\eqn{n, p}).}

\item{fm}{List with the outputs of the RR ((\code{V}): eigenvector matrix of the correlation matrix (n,n); (\code{TtDY}): the Y-scores matrix (n,q); (\code{sv}): eigen values of the matrix (1,n); (\code{lb}): value of regularization parameter \eqn{lambda}; (\code{xmeans}): the centering vector of X (p,1); (\code{ymeans}): the centering vector of Y (q,1); (\code{weights}): the weights vector of X-variables (p,1)}

\item{K}{kernel Gram matrix} 

\item{kern}{kernel function} 

\item{dots}{Optional arguments passed in the kernel function}

For \code{predict.Dkrr}: 

\item{pred}{A list of matrices (\eqn{m, q}) with the Y predicted values for the new X-data}

\item{K}{kernel Gram matrix (\eqn{m, nlv}), with values for the new X-data}

For \code{coef.Dkrr}: 

\item{int}{matrix (1,nlv) with the intercepts}

\item{B}{matrix (n,nlv) with the coefficients}

\item{df}{}
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

lb <- 2
fm <- dkrr(Xtrain, Ytrain, lb = lb, kern = "krbf", gamma = .8)
coef(fm)
coef(fm, lb = .6)
predict(fm, Xtest)
predict(fm, Xtest, lb = c(0.1, .8))

pred <- predict(fm, Xtest)$pred
msep(pred, Ytest)

lb <- 2
fm <- dkrr(Xtrain, Ytrain, lb = lb, kern = "kpol", degree = 2, coef0 = 10)
predict(fm, Xtest)

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

fm <- dkrr(X, y, lb = .01, gamma = .5)
pred <- predict(fm, X)$pred
plot(X, y, type = "p")
lines(X, zy, lty = 2)
lines(X, pred, col = "red")

}

\keyword{datagen}