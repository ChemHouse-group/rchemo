\name{vip}
\alias{vip}
\encoding{UTF-8}

\title{Variable Importance in Projection (VIP)}

\description{

\code{vip} calculates the Variable Importance in Projection (VIP) for a PLS model. 

}

\usage{

vip(object, X, Y = NULL, nlv = NULL)

}

\arguments{

\item{object}{A fitted model, output of a call to a fitting function among \code{plskern}, \code{plsnipals}, \code{plsrannar}, \code{plsrda}, \code{plslda}), \code{plsqda}).}

\item{X}{X-data involved in the fitted model}

\item{Y}{Y-data involved in the fitted model.
If \code{Y} is NULL (default value), the VIP calculation is based on the proportion of Y-variance explained by the components, as proposed by Mehmood et al (2012, 2020).
If \code{Y} is not NULL, the VIP calculation is based on the redundancy, as proposed by Tenenhaus (1998).}

\item{nlv}{Number of components (LVs) to consider.}

}

\value{

matrix (\eqn{(q,nlv)}) with VIP values, for models with 1 to nlv latent variables.

}

\references{

Mehmood, T.,Liland, K.H.,Snipen, L.,Sæbø, S., 2012. A review of variable selection methods in Partial Least Squares Regression. Chemometrics and Intelligent Laboratory Systems, 118, 62-69. 

Mehmood, T., Sæbø, S.,Liland, K.H., 2020. Comparison of variable selection methods in partial least squares regression. Journal of Chemometrics, 34, e3226.

Tenenhaus, M., 1998. La régression PLS: théorie et pratique. Editions Technip, Paris, France.

}


\examples{

## EXAMPLE OF PLS

n <- 50 ; p <- 4
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- rnorm(n)
Ytrain <- cbind(y1 = ytrain, y2 = 100 * ytrain)
m <- 3
Xtest <- Xtrain[1:m, , drop = FALSE] 
Ytest <- Ytrain[1:m, , drop = FALSE] ; ytest <- Ytest[1:m, 1]

nlv <- 3
fm <- plskern(Xtrain, Ytrain, nlv = nlv)
vip(fm, Xtrain, Ytrain, nlv = nlv)
vip(fm, Xtrain, nlv = nlv)

fm <- plskern(Xtrain, ytrain, nlv = nlv)
vip(fm, Xtrain, ytrain, nlv = nlv)
vip(fm, Xtrain, nlv = nlv)

## EXAMPLE OF PLSDA

n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- sample(c("1", "4", "10"), size = n, replace = TRUE)

Xtest <- Xtrain[1:5, ] ; ytest <- ytrain[1:5]

nlv <- 5
fm <- plsrda(Xtrain, ytrain, nlv = nlv)
vip(fm, Xtrain, ytrain, nlv = nlv)

}

\keyword{datagen}