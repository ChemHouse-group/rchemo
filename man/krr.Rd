\name{krr}
\alias{krr}
\alias{coef.Krr}
\alias{predict.Krr}
\encoding{latin1}

\title{KRR (LS-SVMR)}

\description{

Kernel ridge regression models (KRR = LS-SVMR) (Suykens et al. 2000, Bennett & Embrechts 2003, Krell 2018). 

}

\usage{

krr(X, Y, weights = NULL, lb = 1e-2, kern = "krbf", ...)

\method{coef}{Krr}(object, ..., lb = NULL)  

\method{predict}{Krr}(object, X, ..., lb = NULL)  

}

\arguments{

\item{X}{For main function: Training X-data (\eqn{n, p}). --- For auxiliary function: New X-data (\eqn{m, p}) to consider.}

\item{Y}{Training Y-data (\eqn{n, q}).}

\item{weights}{Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{lb}{A value of regularization parameter \eqn{lambda}. If \code{lb = 0}, a pseudo-inverse is used in the RR.} 

\item{kern}{Name of the function defining the considered kernel for building the Gram matrix. See \code{\link{krbf}} for syntax, and other available kernel functions.}

\item{...}{Optional arguments to pass in the kernel function defined in \code{kern} (e.g. \code{gamma} for \code{\link{krbf}}).}

\item{object}{--- For auxiliary function: A fitted model, output of a call to the main function.}

}

\value{

For \code{krr}: 

\item{X}{Training X-data (\eqn{n, p}).}

\item{K}{Gram matrix}

\item{Kt}{Gram matrix}

\item{U}{intermediate output.}

\item{UtDY}{intermediate output.}

\item{sv}{singular values of the matrix (1,n)}

\item{lb}{value of regularization parameter \eqn{lambda}}

\item{ymeans}{the centering vector of Y (q,1)}

\item{weights}{the weights vector of X-variables (p,1)}

\item{kern}{kern function.}

\item{dots}{Optional arguments.}

For \code{coef.Krr}:

\item{int}{matrix (1,nlv) with the intercepts}

\item{alpha}{matrix (n,nlv) with the coefficients}

\item{df}{model complexity (number of degrees of freedom)}

For \code{predict.Krr}:

\item{pred}{A list of matrices (\eqn{m, q}) with the Y predicted values for the new X-data}

}

\references{

Bennett, K.P., Embrechts, M.J., 2003. An optimization perspective on kernel partial least squares regression, in: Advances in Learning Theory: Methods, Models and Applications, NATO Science Series III: Computer & Systems Sciences. IOS Press Amsterdam, pp. 227-250.

Cawley, G.C., Talbot, N.L.C., 2002. Reduced Rank Kernel Ridge Regression. Neural Processing Letters 16, 293-302. https://doi.org/10.1023/A:1021798002258

Krell, M.M., 2018. Generalizing, Decoding, and Optimizing Support Vector Machine Classification. arXiv:1801.04929.

Saunders, C., Gammerman, A., Vovk, V., 1998. Ridge Regression Learning Algorithm in Dual Variables, in: In Proceedings of the 15th International Conference on Machine Learning. Morgan Kaufmann, pp. 515-521.

Suykens, J.A.K., Lukas, L., Vandewalle, J., 2000. Sparse approximation using least squares support vector machines. 2000 IEEE International Symposium on Circuits and Systems. Emerging Technologies for the 21st Century. Proceedings (IEEE Cat No.00CH36353). https://doi.org/10.1109/ISCAS.2000.856439

Welling, M., n.d. Kernel ridge regression. Department of Computer Science, University of Toronto, Toronto, Canada. https://www.ics.uci.edu/~welling/classnotes/papers_class/Kernel-Ridge.pdf

}

\note{

KRR is close to the particular SVMR setting the \eqn{epsilon} coefficient to zero (no marges excluding observations). The difference is that a L2-norm optimization is done, instead L1 in SVM.   

The second example concerns the fitting of the function sinc(x) described in Rosipal & Trejo 2001 p. 105-106

}

\examples{

## EXAMPLE 1

n <- 6 ; p <- 4
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- rnorm(n)
Ytrain <- cbind(y1 = ytrain, y2 = 100 * ytrain)
m <- 3
Xtest <- Xtrain[1:m, , drop = FALSE] 
Ytest <- Ytrain[1:m, , drop = FALSE] ; ytest <- Ytest[1:m, 1]

lb <- 2
fm <- krr(Xtrain, Ytrain, lb = lb, kern = "krbf", gamma = .8)
coef(fm)
coef(fm, lb = .6)
predict(fm, Xtest)
predict(fm, Xtest, lb = c(0.1, .6))

pred <- predict(fm, Xtest)$pred
msep(pred, Ytest)

lb <- 2
fm <- krr(Xtrain, Ytrain, lb = lb, kern = "kpol", degree = 2, coef0 = 10)
predict(fm, Xtest)

## EXAMPLE 2

x <- seq(-10, 10, by = .2)
x[x == 0] <- 1e-5
n <- length(x)
zy <- sin(abs(x)) / abs(x)
y <- zy + rnorm(n, 0, .2)
plot(x, y, type = "p")
lines(x, zy, lty = 2)
X <- matrix(x, ncol = 1)

fm <- krr(X, y, lb = .1, gamma = .5)
pred <- predict(fm, X)$pred
plot(X, y, type = "p")
lines(X, zy, lty = 2)
lines(X, pred, col = "red")

}

\keyword{datagen}