\name{cglsr}
\alias{cglsr}
\alias{coef.Cglsr}
\alias{predict.Cglsr}
\encoding{latin1}

\title{CG Least Squares Models}

\description{

Conjugate gradient algorithm (CG) for the normal equations (CGLS algorithm 7.4.1, Bjorck 1996, p.289)

}

\usage{

cglsr(X, y, nlv, reorth = TRUE, filt = FALSE)

\method{coef}{Cglsr}(object, ..., nlv = NULL)  

\method{predict}{Cglsr}(object, X, ..., nlv = NULL)  

}

\arguments{

\item{X}{For the main function: Training X-data (\eqn{n, p}). --- For auxiliary functions: New X-data (\eqn{m, p}) to consider. }

\item{y}{Univariate training Y-data (\eqn{n, 1}).}

\item{nlv}{The number(s) of CG iterations.}

\item{reorth}{Logical. If \code{TRUE}, a Gram-Schmidt reorthogonalization of the normal equation residual vectors is done.}

\item{filt}{Logical. If \code{TRUE}, the filter factors are computed (output \code{F}).}

\item{object}{For auxiliary functions: A fitted model, output of a call to the main functions.}

\item{...}{For auxiliary functions: Optional arguments. Not used.}

}

\details{

The code for re-orthogonalization (Hansen 1998) and filter factors (Vogel 1987, Hansen 1998) computations is a transcription (with few adaptations) of the matlab function `cgls` (Saunders et al. https://web.stanford.edu/group/SOL/software/cgls/; Hansen 2008).

The filter factors can be used to compute the model complexity of CGLSR and PLSR models (see \code{\link{dfplsr_cg}}).

Data \eqn{X} and \eqn{y} are internally centered. 

Missing values are not allowed.

}

\value{	

For \code{cglsr}: 

  \item{B}{matrix with the model coefficients for the fix nlv.}
  
	\item{gnew}{squared norm of the s vector}
	
	\item{xmeans}{variable means for the training X-data}
	
	\item{ymeans}{variable means for the training Y-data}
	
	\item{F}{If \code{filt = TRUE}, the filter factors}
	
	
For \code{coef.Cglsr} : 

\item{int}{intercept value.}

\item{B}{matrix with the model coefficients.}

	
For \code{predict.Cglsr} : 

  \item{pred}{list of matrices, with the predicted values for each number \code{nlv} of CG iterations}
	
}


\references{

Bjorck, A., 1996. Numerical Methods for Least Squares Problems, Other Titles in Applied Mathematics. 
Society for Industrial and Applied Mathematics. https://doi.org/10.1137/1.9781611971484

Hansen, P.C., 1998. Rank-Deficient and Discrete Ill-Posed Problems, Mathematical Modeling and Computation. 
Society for Industrial and Applied Mathematics. https://doi.org/10.1137/1.9780898719697

Hansen, P.C., 2008. Regularization Tools version 4.0 for Matlab 7.3. 
Numer Algor 46, 189-194. https://doi.org/10.1007/s11075-007-9136-9

Manne R. Analysis of two partial-least-squares algorithms for multivariate calibration. Chemometrics Intell.
Lab. Syst. 1987; 2: 187-197.

Phatak A, De Hoog F. Exploiting the connection between
PLS, Lanczos methods and conjugate gradients: alternative proofs of some properties of PLS. J. Chemometrics
2002; 16: 361-367.

Vogel, C. R.,  "Solving ill-conditioned linear systems using the conjugate gradient method", 
Report, Dept. of Mathematical Sciences, Montana State University, 1987.

}

\examples{

z <- ozone$X
u <- which(!is.na(rowSums(z)))
X <- z[u, -4]
y <- z[u, 4]
dim(X)
headm(X)
Xtest <- X[1:2, ]
ytest <- y[1:2]

nlv <- 10
fm <- cglsr(X, y, nlv = nlv)

coef(fm)
coef(fm, nlv = 1)

predict(fm, Xtest)
predict(fm, Xtest, nlv = 1:3)

pred <- predict(fm, Xtest)$pred
msep(pred, ytest)

cglsr(X, y, nlv = 5, filt = TRUE)$F


}

\keyword{datagen}