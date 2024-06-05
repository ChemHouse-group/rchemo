\name{mbplsr}
\alias{mbplsr}
\alias{transform.Mbplsr}
\alias{summary.Mbplsr}
\alias{coef.Mbplsr}
\alias{predict.Mbplsr}
\encoding{UTF-8}

\title{multi-block PLSR algorithms}

\description{

Algorithm fitting a multi-block PLS1 or PLS2 model between dependent variables \eqn{Xlist} and responses \eqn{Y}, based on the "Improved kernel algorithm #1" proposed by Dayal and MacGregor (1997). 

For weighted versions, see for instance Schaal et al. 2002, Siccard & Sabatier 2006, Kim et al. 2011 and Lesnoff et al. 2020.

\bold{Auxiliary functions}

\code{transform} Calculates the LVs for any new matrix \eqn{X} from the model.

\code{summary} returns summary information for the model.

\code{coef} Calculates b-coefficients from the model, adjuted for raw data.

\code{predict} Calculates the predictions for any new matrix \eqn{X} from the model.

}

\usage{

mbplsr(Xlist, Y, blockscaling = TRUE, weights = NULL, nlv, 
Xscaling = c("none", "pareto", "sd")[1], Yscaling = c("none", "pareto", "sd")[1])

\method{transform}{Mbplsr}(object, Xlist, ..., nlv = NULL)  

\method{summary}{Mbplsr}(object, Xlist, ...)  

\method{coef}{Mbplsr}(object, ..., nlv = NULL) 

\method{predict}{Mbplsr}(object, Xlist, ..., nlv = NULL)  

}

\arguments{

\item{Xlist}{For the main functions: list of training X-data (\eqn{n}rows). --- For the auxiliary functions: list of new X-data, with the same variables than the training X-data.}

\item{Y}{Training Y-data (\eqn{n, q}).}

\item{blockscaling}{logical. If TRUE, the scaling factor (computed on the training) is the "norm" of the block, i.e. the square root of the sum of the variances of each column of the block.}

\item{weights}{Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{nlv}{For the main functions: The number(s) of LVs to calculate. --- For the auxiliary functions: The number(s) of LVs to consider.}

\item{Xscaling}{vector (of length Xlist) of variable scaling for each datablock, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{Yscaling}{character. variable scaling for the Y-block, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main functions.}

\item{...}{For the auxiliary functions: Optional arguments. Not used.}

}

\value{

A list of outputs, such as

\item{T}{The X-score matrix (\eqn{n, nlv}).}

\item{P}{The X-loadings matrix (\eqn{p, nlv}).}

\item{W}{The X-loading weights matrix (\eqn{p, nlv}).}

\item{C}{The Y-loading weights matrix (C = t(Beta), where Beta is the scores regression coefficients matrix).}

\item{R}{The PLS projection matrix (\eqn{p, nlv}).}

\item{xmeans}{The list of centering vectors of \eqn{Xlist}.}

\item{ymeans}{The centering vector of \eqn{Y} (\eqn{q, 1}).}

\item{xscales}{The list of \eqn{Xlist} variable standard deviations.}

\item{yscales}{The vector of \eqn{Y} variable standard deviations (\eqn{q, 1}).}

\item{weights}{Weights applied to the training observations.}

\item{TT}{the X-score normalization factor.}

\item{blockscaling}{block scaling.}

\item{Xnorms}{"norm" of each block, i.e. the square root of the sum of the variances of each column of each block, computed on the training, and used as scaling factor}.

\item{U}{intermediate output.}

For \code{transform.Mbplsr}: X-scores matrix for new Xlist-data.

For \code{summary.Mbplsr}:

\item{explvarx}{matrix of explained variances.}

For \code{coef.Mbplsr}: 

\item{int}{matrix (1,nlv) with the intercepts}

\item{B}{matrix (n,nlv) with the coefficients}

For \code{predict.Mbplsr}: 

\item{pred}{A list of matrices (\eqn{m, q}) with the Y predicted values for the new Xlist-data}

}

\references{

Andersson, M., 2009. A comparison of nine PLS1 algorithms. Journal of Chemometrics 23, 518-529.

Dayal, B.S., MacGregor, J.F., 1997. Improved PLS algorithms. Journal of Chemometrics 11, 73-85.

Hoskuldsson, A., 1988. PLS regression methods. Journal of Chemometrics 2, 211-228. https://doi.org/10.1002/cem.1180020306

Kim, S., Kano, M., Nakagawa, H., Hasebe, S., 2011. Estimation of active pharmaceutical ingredients content using locally weighted partial least squares and statistical wavelength selection. Int. J. Pharm., 421, 269-274.

Lesnoff, M., Metz, M., Roger, J.M., 2020. Comparison of locally weighted PLS strategies for regression and discrimination on agronomic NIR Data. Journal of Chemometrics. e3209. https://onlinelibrary.wiley.com/doi/abs/10.1002/cem.3209

Rannar, S., Lindgren, F., Geladi, P., Wold, S., 1994. A PLS kernel algorithm for data sets with many variables and fewer objects. Part 1: Theory and algorithm. Journal of Chemometrics 8, 111-125. https://doi.org/10.1002/cem.1180080204

Schaal, S., Atkeson, C., Vijayamakumar, S. 2002. Scalable techniques from nonparametric statistics for the real time robot learning. Applied Intell., 17, 49-60.

Sicard, E. Sabatier, R., 2006. Theoretical framework for local PLS1 regression and application to a rainfall data set. Comput. Stat. Data Anal., 51, 1393-1410.

Tenenhaus, M., 1998. La régression PLS: théorie et pratique. Editions Technip, Paris, France.

Wold, S., Sjostrom, M., Eriksson, l., 2001. PLS-regression: a basic tool for chemometrics. Chem. Int. Lab. Syst., 58, 109-130.

}

\examples{

n <- 10 ; p <- 10
Xtrain <- matrix(rnorm(n * p), ncol = p)
ytrain <- rnorm(n)

m <- 2
Xtest <- matrix(rnorm(m * p), ncol = p)

colnames(Xtrain) <- colnames(Xtest) <- paste("v", 1:p, sep = "")

Xtrain
Xtest

blocks <- list(1:2, 4, 6:8)
X1 <- mblocks(Xtrain, blocks = blocks)
X2 <- mblocks(Xtest, blocks = blocks)

nlv <- 3
fm <- mbplsr(Xlist = X1, Y = ytrain, Xscaling = c("sd","none","none"), 
blockscaling = TRUE, weights = NULL, nlv = nlv)

summary(fm, X1)
coef(fm)
transform(fm, X2)
predict(fm, X2)

}

\keyword{datagen}