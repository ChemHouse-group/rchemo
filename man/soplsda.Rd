\name{soplsrda}
\alias{soplsrda}
\alias{soplslda}
\alias{soplsqda}
\alias{transform.Soplsrda} 
\alias{transform.Soplsprobda}
\alias{predict.Soplsrda}
\alias{predict.Soplsprobda} 
\alias{soplsrdacv}
\alias{soplsldacv}
\alias{soplsqdacv}

\encoding{UTF-8}

\title{Block dimension reduction by SO-PLS-DA}

\description{

Function \code{soplsrda} implements dimension reductions of pre-selected blocks of variables (= set of columns) of a reference (= training) matrix, by sequential orthogonalization-PLS (said "SO-PLS") in a context of discrimination. 

Function \code{soplsrdacv} perfoms repeteated cross-validation of an SO-PLS-RDA model in order to choose the optimal lv combination from the different blocks.

The block reduction consists in calculating latent variables (= scores) for each block, each block being sequentially orthogonalized to the information computed from the previous blocks.

The function allows giving a priori weights to the rows of the reference matrix in the calculations.

In\code{soplslda} and \code{soplsqda}, probabilistic LDA and QDA are run over the PLS2 LVs, respectively.

}

\usage{

soplsrda(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], 
Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlv)

soplslda(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], 
Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlv, 
prior = c("unif", "prop"))

soplsqda(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], 
Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlv, 
prior = c("unif", "prop"))

soplsrdacv(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], 
Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlvlist=list(), 
nbrep=30, cvmethod="kfolds", seed = 123, samplingk = NULL, nfolds = 7, 
optimisation = c("global","sequential")[1], 
criterion = c("err","rmse")[1], selection = c("localmin","globalmin","1std")[1])

soplsldacv(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], 
Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlvlist=list(), 
prior = c("unif", "prop"), nbrep = 30, cvmethod = "kfolds", seed = 123, samplingk = NULL, 
nfolds = 7, optimisation = c("global","sequential")[1], 
criterion = c("err","rmse")[1], selection = c("localmin","globalmin","1std")[1])

soplsqdacv(Xlist, y, Xscaling = c("none", "pareto", "sd")[1], 
Yscaling = c("none", "pareto", "sd")[1], weights = NULL, nlvlist = list(), 
prior = c("unif", "prop"), nbrep = 30, cvmethod = "kfolds", seed = 123, samplingk = NULL, 
nfolds = 7, optimisation = c("global","sequential")[1], 
criterion = c("err","rmse")[1], selection = c("localmin","globalmin","1std")[1])

\method{transform}{Soplsrda}(object, X, ...) 

\method{transform}{Soplsprobda}(object, X, ...) 

\method{predict}{Soplsrda}(object, X, ...) 

\method{predict}{Soplsprobda}(object, X, ...) 

}

\arguments{

\item{Xlist}{For the main functions: A list of matrices or data frames of reference (= training) observations.}

\item{X}{For the auxiliary functions: list of new X-data, with the same variables than the training X-data.}

\item{y}{Training class membership (\eqn{n}). \bold{Note:} If \code{y} is a factor, it is replaced by a character vector.}

\item{Xscaling}{vector (of length Xlist) of variable scaling for each datablock, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{Yscaling}{variable scaling for the Y-block, among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{weights}{a priori weights to the rows of the reference matrix in the calculations.}

\item{nlv}{A vector of same length as the number of blocks defining the number of scores to calculate for each block, or a single number. In this last case, the same number of scores is used for all the blocks.}

\item{nlvlist}{A list of same length as the number of X-blocks. Each component of the list gives the number of PLS components of the corresponding X-block to test.}

\item{nbrep}{An integer, setting the number of CV repetitions. Default value is 30.}

\item{cvmethod}{"kfolds" for k-folds cross-validation, or "loo" for leave-one-out.}

\item{seed}{a numeric. Seed used for the repeated resampling, and if cvmethod is "kfolds" and samplingk is not NULL.}

\item{samplingk}{Optional. A vector of length n. The elements are the values of a qualitative variable used for stratified partition creation.}

\item{nfolds}{An integer, setting the number of partitions to create. Default value is 7.}

\item{optimisation}{"global" or "sequential" optimisation of the number of components. If "sequential", the optimal lv number is found for the first X-block, then for the 2nd one, etc...}

\item{criterion}{optimisation criterion among "rmse" and "err" (for classification error rate)}

\item{selection}{a character indicating the selection method to use to choose the optimal combination of components, among "localmin","globalmin","1std". If "localmin": the optimal combination corresponds to the first local minimum of the mean CV rmse or error rate. If "globalmin" : the optimal combination corresponds to the minimum mean CV rmse or error rate. If "1std" (one standard error rule) : it corresponds to the first combination after which the mean cross-validated rmse or error rate does not decrease significantly.}


\item{prior}{The prior probabilities of the classes. Possible values are "unif" (default; probabilities are set equal for all the classes) or "prop" (probabilities are set equal to the observed proportions of the classes in \code{y}).}

\item{object}{For the auxiliary functions: A fitted model, output of a call to the main functions.}

\item{...}{For the auxiliary functions: Optional arguments. Not used.}

}

\value{

For \code{soplsrda}, \code{soplslda}, \code{soplsqda}:

\item{fm}{list with the PLS models: (\code{T}): X-scores matrix; (\code{P}): X-loading matrix;(\code{R}): The PLS projection matrix (p,nlv); (\code{W}): X-loading weights matrix ;(\code{C}): The Y-loading weights matrix; (\code{TT}): the X-score normalization factor; (\code{xmeans}): the centering vector of X (p,1);  (\code{ymeans}): the centering vector of Y (q,1); (\code{weights}): vector of observation weights; (\code{Xscales}): X scaling values; (\code{Yscales}): Y scaling values; (\code{U}): intermediate output.}

\item{lev}{classes}

\item{ni}{number of observations in each class}

For \code{transform.Soplsrda}, \code{transform.Soplsprobda}:  the LVs Calculated for the new matrices list \eqn{Xlist} from the model.

For \code{predict.Soplsrda}, \code{predict.Soplsprobda}: 

\item{pred}{predicted class for each observation}

\item{posterior}{calculated probability of belonging to a class for each observation}

For \code{soplsrdacv}, \code{soplsldacv}, \code{soplsqdacv}:

\item{lvcombi}{matrix or list of matrices, of tested component combinations.}

\item{optimCombiLine}{number of the combination line corresponding to the optimal one. In the case of a sequential optimisation, it is the number of the combination line in the model with all the X-blocks.}

\item{optimcombi}{the number of PLS components of each X-block allowing the optimisation of the mean rmseCV.}

\item{optimExplVarCV}{cross-validated explained variance for the optimal soplsda model.}

\item{rmseCV}{matrix or list of matrices of mean and sd of cross-validated rmse in the model for each combination and response variables.}

\item{ExplVarCV}{matrix or list of matrices of mean and sd of cross-validated explained variances in the model for each combination and response variables.}

\item{errCV}{matrix or list of matrices of mean and sd of cross-validated classification error rates in the model for each combination and response variables.}

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
y <- matrix(sample(c("1", "4", "10"), size = N, replace = TRUE), ncol=1)
colnames(X) <- paste("x", 1:ncol(X), sep = "")
set.seed(NULL)

n <- nrow(X)

X_list <- list(X[,1:4], X[,5:7], X[,9:ncol(X)])
X_list_2 <- list(X[1:2,1:4], X[1:2,5:7], X[1:2,9:ncol(X)])

# EXEMPLE WITH SO-PLS-RDA
soplsrdacv(X_list, y, Xscaling = c("none", "pareto", "sd")[1], 
Yscaling = c("none", "pareto", "sd")[1], weights = NULL,
nlvlist=list(0:1, 1:2, 0:1), nbrep=1, cvmethod="loo", seed = 123, 
samplingk = NULL, nfolds = 3, optimisation = "global", 
criterion = c("err","rmse")[1], selection = "localmin")

ncomp <- 2
fm <- soplsrda(X_list, y, nlv = ncomp)
predict(fm,X_list_2)
transform(fm,X_list_2)

ncomp <- c(2, 0, 3)
fm <- soplsrda(X_list, y, nlv = ncomp)
predict(fm,X_list_2)
transform(fm,X_list_2)

ncomp <- 0
fm <- soplsrda(X_list, y, nlv = ncomp)
predict(fm,X_list_2)
transform(fm,X_list_2)

# EXEMPLE WITH SO-PLS-LDA
ncomp <- 2
weights <- rep(1 / n, n)
#w <- 1:n
soplslda(X_list, y, Xscaling = "none", nlv = ncomp, weights = weights)
soplslda(X_list, y, Xscaling = "pareto", nlv = ncomp, weights = weights)
soplslda(X_list, y, Xscaling = "sd", nlv = ncomp, weights = weights)

fm <- soplslda(X_list, y, Xscaling = c("none","pareto","sd"), nlv = ncomp, weights = weights)
predict(fm,X_list_2)
transform(fm,X_list_2)

}

\keyword{datagen}