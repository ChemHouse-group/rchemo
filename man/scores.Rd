\name{mse}
\alias{mse}
\alias{msep}
\alias{rmsep}
\alias{sep}
\alias{bias}
\alias{cor2}
\alias{r2}
\alias{rpd}
\alias{rpdr}
\alias{err}
\alias{residreg}
\alias{residcla}
\encoding{UTF-8}

\title{Residuals and prediction error rates}

\description{

Residuals and prediction error rates (MSEP, SEP, etc. or classification error rate)  for models with quantitative or qualitative responses. 

}

\usage{

residreg(pred, Y)
residcla(pred, y)

msep(pred, Y)
rmsep(pred, Y)
sep(pred, Y)
bias(pred, Y)
cor2(pred, Y)
r2(pred, Y)
rpd(pred, Y)
rpdr(pred, Y)
mse(pred, Y, digits = 3)

err(pred, y)

}

\arguments{

\item{pred}{Prediction (\eqn{m, q}); output of a function \code{predict}.}

\item{Y}{Observed response (\eqn{m, q}).}

\item{y}{Observed response (\eqn{m, 1}).}

\item{digits}{Number of digits for the numerical outputs.}

}

\details{

The rate \eqn{R2} is calculated by \eqn{R2 = 1 - MSEP(current model) / MSEP(null model)}, where  \eqn{MSEP = Sum((y_i - pred_i)^2)/n} and "null model" is the overall mean of \eqn{y}. For predictions over CV or Test sets, and/or for non linear models, it can be different from the square of the correlation coefficient (\eqn{cor2}) between the observed values and the predictions. 

Function \code{sep} computes the SEP, referred to as "corrected SEP" (SEP_c) in Bellon et al. 2010. SEP is the standard deviation of the residuals. There is the relation: \eqn{MSEP = BIAS^2 + SEP^2}.

Function \code{rpd} computes the ratio of the "deviation" (sqrt of the mean of the squared residuals for the null model when it is defined by the simple average) to the "performance" (sqrt of the mean of the squared residuals for the current model, i.e. RMSEP), i.e. \eqn{RPD = SD / RMSEP = RMSEP(null model) / RMSEP} (see eg. Bellon et al. 2010).

Function \code{rpdr} computes a robust RPD. 

}


\value{

Residuals or prediction error rates.

}

\references{

Bellon-Maurel, V., Fernandez-Ahumada, E., Palagos, B., Roger, J.-M., McBratney, A., 2010. Critical review of chemometric indicators commonly used for assessing the quality of the prediction of soil attributes by NIR spectroscopy. TrAC Trends in Analytical Chemistry 29, 1073–1081. https://doi.org/10.1016/j.trac.2010.05.006

}

\examples{
n <- 6 ; p <- 4
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- rnorm(n)
Ytrain <- cbind(y1 = ytrain, y2 = 100 * ytrain)
m <- 3
Xtest <- Xtrain[1:m, , drop = FALSE] 
Ytest <- Ytrain[1:m, , drop = FALSE] 
ytest <- Ytest[1:m, 1]
nlv <- 3
fm <- plskern(Xtrain, Ytrain, nlv = nlv)
pred <- predict(fm, Xtest)$pred

residreg(pred, Ytest)
msep(pred, Ytest)
rmsep(pred, Ytest)
sep(pred, Ytest)
bias(pred, Ytest)
cor2(pred, Ytest)
r2(pred, Ytest)
rpd(pred, Ytest)
rpdr(pred, Ytest)
mse(pred, Ytest, digits = 3)
}

\examples{
n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- sample(c(1, 4, 10), size = n, replace = TRUE)
Xtest <- Xtrain[1:5, ]
ytest <- ytrain[1:5]
nlv <- 5
fm <- plsrda(Xtrain, ytrain, nlv = nlv)
pred <- predict(fm, Xtest)$pred

residcla(pred, ytest)
err(pred, ytest)

}

\keyword{datagen}