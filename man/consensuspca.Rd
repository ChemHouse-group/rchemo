\name{consensuspca}
\alias{consensuspca}
\alias{transform.Consensuspca}
\alias{summary.Consensuspca}

\encoding{UTF-8}

\title{consensus PCA}

\description{

Algorithms fitting a consensus PCA of a list of matrices \eqn{Xlist}. A chosen PCA algorithm is applied on the X-matrix obtained after variable scaling, blockscaling, block concatenation.

\bold{Auxiliary functions}

\code{transform} Calculates the principal components for any new matrix \eqn{X} from the model.

\code{summary} returns summary information for the model.

}

\usage{

consensuspca(Xlist, blockscaling = TRUE, weights = NULL, nlv, 
      Xscaling = c("none", "pareto", "sd")[1], 
      algo = c("svd","eigen","eigenk","nipals","nipalsna","sph")[1], 
      gs = TRUE, tol = .Machine$double.eps^0.5, maxit = 200)

\method{transform}{Consensuspca}(object, X, ..., nlv = NULL) 

\method{summary}{Consensuspca}(object, X, ...) 

}

\arguments{

\item{Xlist}{For the main function: list of training X-data (\eqn{n}rows).}

\item{X}{For the auxiliary functions: list of new X-data, with the same variables than the training X-data.}

\item{blockscaling}{logical. If TRUE, the scaling factor (computed on the training) is the "norm" of the block, i.e. the square root of the sum of the variances of each column of the block.}

\item{weights}{excepted for "nipalsna". Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{nlv}{For the main functions: The number(s) of PCs to calculate. --- For the auxiliary functions: The number(s) of PCs to consider.}

\item{Xscaling}{vector (of length Xlist) of variable scaling for each datablock, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{algo}{For \code{consensuspca}, the algorithm for the PCA among : "svd" (SVD factorization of \eqn{D^(1/2) * X}, using function \code{\link{svd}}), "eigen" (Eigen factorization of \eqn{X' * D * X}, using function \code{\link{eigen}}), "eigenk" (Eigen factorization of \eqn{D^(1/2) * X * X' D^(1/2)}, using function \code{\link{eigen}}), "nipals" (Eigen factorization of \eqn{X' * D * X} using NIPALS), "nipalsna" (Eigen factorization of \eqn{X' * D * X} using NIPALS allowing missing data in \eqn{X}), "sph" (Robust spherical PCA)}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main functions.}

\item{...}{For the auxiliary functions: Optional arguments. Not used.}

\bold{Specific for the NIPALS algorithm}

\item{gs}{Logical indicating if a Gram-Schmidt orthogonalization is implemented or not (default to \code{TRUE}).}

\item{tol}{Tolerance for testing convergence of the NIPALS iterations for each PC.}

\item{maxit}{Maximum number of NIPALS iterations for each PC.}

}

\value{

For \code{consensuspca}: 

\item{T}{The X-score matrix (\eqn{n, nlv}).}

\item{P}{The X-loadings matrix (\eqn{p, nlv}).}

\item{sv}{The singular values (\eqn{min(n, p), 1}) except for NIPALS = (\eqn{nlv, 1}).}

\item{eig}{The eigenvalues (\code{= sv^2}) (\eqn{min(n, p), 1}) except for NIPALS = (\eqn{nlv, 1}).}

\item{xmeans}{The list of centering vectors of \eqn{Xlist}.}

\item{xscales}{The list of \eqn{Xlist} variable standard deviations.}

\item{weights}{Weights applied to the training observations.}

\item{blockscaling}{block scaling.}

\item{Xnorms}{"norm" of each block, i.e. the square root of the sum of the variances of each column of each block, computed on the training, and used as scaling factor}.

\item{niter}{Numbers of iterations of the NIPALS.}

\item{conv}{Logical indicating if the NIPALS converged before reaching the maximal number of iterations.}



For \code{transform.Consensuspca}: X-scores matrix for new Xlist-data.


For \code{summary.Consensuspca}:

\item{explvarx}{matrix of explained variances.}

\item{contr_ind}{observation contributions.}

\item{contr_var}{variable contributions.}

\item{coord_var}{variable coordinates.}

\item{cor_circle}{variable coordinates on the correlation circle.}

}

\references{

Mangamana, E.T., Cariou, V., Vigneau, E., Glèlè Kakaï, R.L., Qannari, E.M., 2019. Unsupervised multiblock data analysis: A unified approach and extensions. Chemometrics and 
Intelligent Laboratory Systems 194, 103856. https://doi.org/10.1016/j.chemolab.2019.103856

Westerhuis, J.A., Kourti, T., MacGregor, J.F., 1998. Analysis of multiblock and hierarchical PCA and PLS models. Journal of Chemometrics 12, 301–321. 
https://doi.org/10.1002/(SICI)1099-128X(199809/10)12:5<301::AID-CEM515>3.0.CO;2-S

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
fm <- consensuspca(Xlist = X1, Xscaling = c("sd","none","none"), 
blockscaling = TRUE, weights = NULL, nlv = nlv)

summary(fm, X1)
transform(fm, X2)

}

\keyword{datagen}