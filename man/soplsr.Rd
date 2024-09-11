\name{sopls}
\alias{soplsr}
\alias{transform.Soplsr} 
\alias{predict.Soplsr} 
\alias{soplsrcv}

\encoding{UTF-8}

\title{Block dimension reduction by SO-PLS}

\description{

Function \code{soplsr} implements dimension reductions of pre-selected blocks of variables (= set of columns) of a reference (= training) matrix, by sequential orthogonalization-PLS (said "SO-PLS"). 

Function \code{soplsrcv} perfoms repeteated cross-validation of an SO-PLS model in order to choose the optimal lv combination from the different blocks.

SO-PLS is described for instance in Menichelli et al. (2014), Biancolillo et al. (2015) and Biancolillo (2016). 

The block reduction consists in calculating latent variables (= scores) for each block, each block being sequentially orthogonalized to the information computed from the previous blocks.

The function allows giving a priori weights to the rows of the reference matrix in the calculations.

\bold{Auxiliary functions}

\code{transform} Calculates the LVs for any new matrices list \eqn{Xlist} from the model.

\code{predict} Calculates the predictions for any new matrices list \eqn{Xlist} from the model.

}

\usage{

soplsr(Xlist, Y, Xscaling = c("none", "pareto", "sd")[1], 
Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlv)

soplsrcv(Xlist, Y, Xscaling = c("none", "pareto", "sd")[1], 
Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlvlist = list(), 
nbrep = 30, cvmethod = "kfolds", seed = 123, samplingk = NULL, nfolds = 7, 
optimisation = c("global","sequential")[1], 
selection = c("localmin","globalmin","1std")[1], majorityvote = FALSE)


\method{transform}{Soplsr}(object, X, ...)  

\method{predict}{Soplsr}(object, X, ...)  

}

\arguments{

\item{Xlist}{A list of matrices or data frames of reference (= training) observations.}

\item{X}{For the auxiliary functions: list of new X-data, with the same variables than the training X-data.}

\item{Y}{A \eqn{n x q} matrix or data frame, or a vector of length \eqn{n}, of reference (= training) responses.}

\item{Xscaling}{vector (of length Xlist) of variable scaling for each datablock, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{Yscaling}{variable scaling for the Y-block, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{weights}{a priori weights to the rows of the reference matrix in the calculations.}

\item{nlv}{A vector of same length as the number of blocks defining the number of scores to calculate for each block, or a single number. In this last case, the same number of scores is used for all the blocks.}

\item{nlvlist}{A list of same length as the number of X-blocks. Each component of the list gives the number of PLS components of the corresponding X-block to test.}

\item{nbrep}{An integer, setting the number of CV repetitions. Default value is 30.}

\item{cvmethod}{"kfolds" for k-folds cross-validation, or "loo" for leave-one-out.}

\item{seed}{a numeric. Seed used for the repeated resampling, and if cvmethod is "kfolds" and samplingk is not NULL.}

\item{samplingk}{A vector of length n. The elements are the values of a qualitative variable used for stratified partition creation. If NULL, the first observation is set in the first fold, the second observation in the second fold, etc...}

\item{nfolds}{An integer, setting the number of partitions to create. Default value is 7.}

\item{optimisation}{"global" or "sequential" optimisation of the number of components. If "sequential", the optimal lv number is found for the first X-block, then for the 2nd one, etc...}

\item{selection}{a character indicating the selection method to use to choose the optimal combination of components, among "localmin","globalmin","1std". If "localmin": the optimal combination corresponds to the first local minimum of the mean CV rmse. If "globalmin" : the optimal combination corresponds to the minimum mean CV rmse. If "1std" (one standard errror rule): it corresponds to the first combination after which the mean cross-validated rmse does not decrease significantly.}

\item{majorityvote}{only if optimisation is "global" or one X-block. If majorityvote is TRUE, the optimal combination is chosen for each Y variable, with the chosen selection, before a majority vote. If majorityvote is "FALSE, the optimal combination is simply chosen with the chosen selection.}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main functions.}

\item{...}{For the auxiliary functions: Optional arguments. Not used.}

}

\value{

For \code{soplsr}:

\item{fm}{A list of the plsr models.}

\item{T}{A matrix with the concatenated scores calculated from the X-blocks.}

\item{pred}{A matrice \eqn{n x q} with the calculated fitted values.}

\item{xmeans}{list of vectors of X-mean values.}

\item{ymeans}{vector of Y-mean values.}

\item{xscales}{list of vectors of X-scaling values.}

\item{yscales}{vector of Y-scaling values.}

\item{b}{A list of X-loading weights, used in the orthogonalization step.}

\item{weights}{Weights applied to the training observations.}

\item{nlv}{vector of numbers of latent variables from each X-block.}

For \code{transform.Soplsr}:  the LVs calculated for the new matrices list \eqn{Xlist} from the model.

For \code{predict.Soplsr}: predicted values for each observation

For \code{soplsrcv}:

\item{lvcombi}{matrix or list of matrices, of tested component combinations.}

\item{optimcombi}{the number of PLS components of each X-block allowing the optimisation of the mean rmseCV.}

\item{rmseCV_byY}{matrix or list of matrices of mean and sd of cross-validated RMSE in the model for each combination and each response variable.}

\item{ExplVarCV_byY}{matrix or list of matrices of mean and sd of cross-validated explained variances in the model for each combination and each response variable.}

\item{rmseCV}{matrix or list of matrices of mean and sd of cross-validated RMSE in the model for each combination and response variables.}

\item{ExplVarCV}{matrix or list of matrices of mean and sd of cross-validated explained variances in the model for each combination and response variables.}

}

\references{

- Biancolillo et al. , 2015. Combining SO-PLS and linear discriminant analysis for
multi-block classification. Chemometrics and Intelligent Laboratory Systems, 141, 58-67.

- Biancolillo, A. 2016. Method development in the area of multi-block analysis focused on food analysis. PhD. University of copenhagen.

- Menichelli et al., 2014. SO-PLS as an exploratory tool for path modelling. Food Quality and Preference, 36, 122-134.

- Tenenhaus, M., 1998. La régression PLS: théorie et pratique. Editions Technip, Paris, France.

}

\seealso{

\code{\link{soplsr_soplsda_allsteps}} function to help determine the optimal number of latent variables, perform a permutation test, calculate model parameters and predict new observations.

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

soplsrcv(X_list, Y, Xscaling = c("none", "pareto", "sd")[1], 
Yscaling = c("none", "pareto", "sd")[1], weights = NULL, 
nlvlist=list(0:1, 1:2, 0:1), nbrep=1, cvmethod="loo", seed = 123, samplingk=NULL,
optimisation="global", selection="localmin", majorityvote=FALSE)


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
fm <- soplsr(X_list, Y, Xscaling = c("sd","pareto","none"), nlv = ncomp, weights = weights)
transform(fm, X_list_2)
predict(fm, X_list_2)

}

\keyword{datagen}