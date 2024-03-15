\name{sopls}
\alias{soplsr}

\encoding{UTF-8}

\title{Block dimension reduction by SO-PLS}

\description{

Function \code{soplsr} implements dimension reductions of pre-selected blocks of variables (= set of columns) of a reference (= training) matrix, by sequential orthogonalization-PLS (said "SO-PLS"). 

SO-PLS is described for instance in Menichelli et al. (2014), Biancolillo et al. (2015) and Biancolillo (2016). 

The block reduction consists in calculating latent variables (= scores) for each block, each block being sequentially orthogonalized to the information computed from the previous blocks.

The function allows giving a priori weights to the rows of the reference matrix in the calculations.

\bold{Auxiliary functions}

\code{transform} Calculates the LVs for any new matrices list \eqn{Xlist} from the model.

\code{predict} Calculates the predictions for any new matrices list \eqn{Xlist} from the model.

}

\usage{

soplsr(Xlist, Y, scaling = c("none", "pareto", "sd")[1], weights = NULL, nlv, ...)

\method{transform}{Soplsr}(object, Xlist, ...)  

\method{predict}{Soplsr}(object, Xlist, ...)  

}

\arguments{

\item{X_list}{A list of matrices or data frames of reference (= training) observations.--- For the auxiliary functions: list of new X-data, with the same variables than the training X-data.}

\item{Y}{A \eqn{n x q} matrix or data frame, or a vector of length \eqn{n}, of reference (= training) responses.}

\item{scaling}{vector (of length Xlist) of variable scaling for each datablock, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "ctreduced", uncorrected standard deviation is used.}

\item{weights}{a priori weights to the rows of the reference matrix in the calculations.}

\item{nlv}{A vector of same length as the number of blocks defining the number of scores to calculate for each block, or a single number. In this last case, the same number of scores is used for all the blocks.}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main functions.}

\item{...}{For the auxiliary functions: Optional arguments. Not used.}

}

\value{

\item{fm}{A list of the plsr models.}

\item{T}{A matrix with the concatenated scores calculated from the X-blocks.}

\item{pred}{A matrice \eqn{n x q} with the calculated fitted values.}

\item{xmeans}{list of vectors of X-mean values.}

\item{ymeans}{vector of Y-mean values.}

\item{xsds}{list of vectors of X-mean standard deviations.}

\item{ysds}{vector of Y-mean standard deviations.}

\item{b}{A list of X-loading weights, used in the orthogonalization step.}

\item{weights}{Weights applied to the training observations.}

\item{Xscaling}{vector of variable scaling for each X-block.}

\item{Yscaling}{scaling for the Y-block.}

\item{nlv}{vector of numbers of latent variables from each X-block.}

}

\references{

- Biancolillo et al. , 2015. Combining SO-PLS and linear discriminant analysis for
multi-block classification. Chemometrics and Intelligent Laboratory Systems, 141, 58-67.

- Biancolillo, A. 2016. Method development in the area of multi-block analysis focused on food analysis. PhD. University of copenhagen.

- Menichelli et al., 2014. SO-PLS as an exploratory tool for path modelling. Food Quality and Preference, 36, 122-134.

- Tenenhaus, M., 1998. La régression PLS: théorie et pratique. Editions Technip, Paris, France.

}


\examples{

N <- 10 ; p <- 12
set.seed(1)
X <- matrix(rnorm(N * p, mean = 10), ncol = p, byrow = TRUE)
Y <- matrix(rnorm(N * 2, mean = 10), ncol = 2, byrow = TRUE)
colnames(X) <- paste("varx", 1:ncol(X), sep = "")
colnames(Y) <- paste("vary", 1:ncol(Y), sep = "")
rownames(X) <- rownames(Y) <- paste("obs", 1:nrow(X), sep = "")
set.seed(NULL)
X
Y

n <- nrow(X)

X_list <- list(X[,1:4], X[,5:7], X[,9:ncol(X)])
X_list_2 <- list(X[1:2,1:4], X[1:2,5:7], X[1:2,9:ncol(X)])

ncomp <- 2
fm <- soplsr(X_list, Y, nlv = ncomp)
transform(fm, X_list_2)
predict(fm, X_list_2)

mse(predict(fm, X_list), Y)

# VIP calculation based on the proportion of Y-variance explained by the components
vip(fm$fm[[1]], X_list[[1]], Y = NULL, nlv = ncomp)
vip(fm$fm[[2]], X_list[[2]], Y = NULL, nlv = ncomp)
vip(fm$fm[[3]], X_list[[3]], Y = NULL, nlv = ncomp)

ncomp <- c(2, 0, 3)
fm <- soplsr(X_list, Y, nlv = ncomp)
transform(fm, X_list_2)
predict(fm, X_list_2)
mse(predict(fm, X_list), Y)

ncomp <- 0
fm <- soplsr(X_list, Y, nlv = ncomp)
transform(fm, X_list_2)
predict(fm, X_list_2)

ncomp <- 2
weights <- rep(1 / n, n)
#w <- 1:n
fm <- soplsr(X_list, Y, scaling = c("sd","pareto","none"), nlv = ncomp, weights = weights)
transform(fm, X_list_2)
predict(fm, X_list_2)

}

\keyword{datagen}