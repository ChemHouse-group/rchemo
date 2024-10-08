\name{svmr}
\alias{svmr}
\alias{svmda}
\alias{predict.Svm}
\alias{summary.Svm}
\encoding{latin1}

\title{SVM Regression and Discrimination}

\description{

SVM models with Gaussian (RBF) kernel.

\code{svmr}: SVM regression (SVMR).

\code{svmda}: SVM discrimination (SVMC).

The SVM models are fitted with parameterization \eqn{'C'}, not the \eqn{'nu'} parameterization. 

The RBF kernel is defined by: exp(-gamma * |x - y|^2).

For tuning the model, usual preliminary ranges are for instance:

- cost = 10^(-5:15)

- epsilon = seq(.1, .3, by = .1)

- gamma = 10^(-6:3)

The functions uses function \code{\link[e1071]{svm}} of package \code{e1071} (Meyer et al. 2021) available on CRAN (e1071 uses the tool box LIVSIM; Chang & Lin, http://www.csie.ntu.edu.tw/~cjlin/libsvm). 

}

\usage{

svmr(X, y, cost = 1, epsilon = .1, gamma = 1, scale = FALSE)

svmda(X, y, cost = 1, epsilon = .1, gamma = 1, scale = FALSE)

\method{predict}{Svm}(object, X, ...)  

\method{summary}{Svm}(object, ...)  

}

\arguments{

\item{X}{For the main functions: Training X-data (\eqn{n, p}). --- For the auxiliary functions: New X-data (\eqn{m, p}) to consider.}

\item{y}{Training Y-data (\eqn{n}).}

\item{cost}{The cost of constraints violation \eqn{cost} parameter. See \code{\link[e1071]{svm}}.}

\item{epsilon}{The \eqn{epsilon} parameter in the insensitive-loss function. See \code{\link[e1071]{svm}}.}

\item{gamma}{The \eqn{gamma} parameter in the RBF kernel.}

\item{scale}{Logical. If \code{TRUE}, \code{X} and \code{Y} are scaled internally.}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main function.}

\item{...}{For the auxiliary functions: Optional arguments.}

}

\value{

For \code{svmr} and \code{svmda}:

\item{fm}{list of outputs such as:
\code{call}; \code{type}; \code{kernel}; \code{cost}; \code{degree}; \code{gamma}; \code{coef0}; \code{nu}; \code{epsilon}; \code{sparse}; \code{scaled}; \code{x.scale}; \code{y.scale}; \code{nclasses}; \code{levels}; \code{tot.nSV}; \code{nSV}; \code{labels}; \code{SV}: The resulting support vectors (possibly scaled); \code{index}: The index of the resulting support vectors in the data matrix. Note that this index refers to the preprocessed data (after the possible effect of na.omit and subset); \code{rho}: The negative intercept; \code{compprob}; \code{probA, probB}: numeric vectors of length k(k-1)/2, k number of classes, containing the parameters of the logistic distributions fitted to the decision values of the binary classifiers (1 / (1 + exp(a x + b))); \code{sigma}: In case of a probabilistic regression model, the scale parameter of the hypothesized (zero-mean) laplace distribution estimated by maximum likelihood; \code{coefs}: The corresponding coefficients times the training labels; \code{na.action}; \code{fitted}; \code{decision.values}; \code{residuals}; \code{isnum}.}

For \code{predict.Svm}:

\item{pred}{predictions for each observation.}

For \code{summary.Svm}:display of call, parameters, and number of support vectors.

}

\references{

Meyer, M. 2021 Support Vector Machines - The Interface to libsvm in package e1071. FH Technikum Wien, Austria, David.Meyer@R-Project.org. https://cran.r-project.org/web/packages/e1071/vignettes/svmdoc.pdf

Chang, cost.-cost. & Lin, cost.-J. (2001). LIBSVM: a library for support vector machines. Software available at http://www.csie.ntu.edu.tw/~cjlin/libsvm. Detailed documentation (algorithms, formulae, . . . ) can be found in http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.ps.gz

}

\note{

The first example illustrates SVMR.
The second one is the example of fitting the function sinc(x) described in Rosipal & Trejo 2001 p. 105-106.
The third one illustrates SVMC.

}

\examples{

## EXAMPLE 1 (SVMR)

n <- 50 ; p <- 4
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- rnorm(n)
m <- 3
Xtest <- Xtrain[1:m, , drop = FALSE] 
ytest <- ytrain[1:m]

fm <- svmr(Xtrain, ytrain)
predict(fm, Xtest)

pred <- predict(fm, Xtest)$pred
msep(pred, ytest)

summary(fm)

## EXAMPLE 2 

x <- seq(-10, 10, by = .2)
x[x == 0] <- 1e-5
n <- length(x)
zy <- sin(abs(x)) / abs(x)
y <- zy + rnorm(n, 0, .2)
plot(x, y, type = "p")
lines(x, zy, lty = 2)
X <- matrix(x, ncol = 1)

fm <- svmr(X, y, gamma = .5)
pred <- predict(fm, X)$pred
plot(X, y, type = "p")
lines(X, zy, lty = 2)
lines(X, pred, col = "red")

## EXAMPLE 3 (SVMC)

n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- sample(c("a", "10", "d"), size = n, replace = TRUE)
m <- 5
Xtest <- Xtrain[1:m, ] ; ytest <- ytrain[1:m]

cost <- 100 ; epsilon <- .1 ; gamma <- 1 
fm <- svmda(Xtrain, ytrain,
    cost = cost, epsilon = epsilon, gamma = gamma)
predict(fm, Xtest)

pred <- predict(fm, Xtest)$pred
err(pred, ytest)

summary(fm)

}

\keyword{datagen}