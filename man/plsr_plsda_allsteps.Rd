\name{plsr_plsda_allsteps}
\encoding{latin1}

\title{PLSR or PLSDA analysis steps}

\description{

Help determine the optimal number of latent variables by cross-validation, perform a permutation test, calculate model parameters and predict new observations, for plsr (\code{\link{plskern}}), plsrda (\code{\link{plsrda}}), plslda (\code{\link{plslda}}) or plsqda (\code{\link{plsqda}}) models.

}

\usage{

plsr_plsda_allsteps(import = c("R","ChemFlow","W4M")[1],
                   X, Xname = NULL, Xscaling = c("none","pareto","sd")[1], 
                   Y, Yscaling = "none", weights = NULL,
                   newX = NULL, newXname = NULL,
                   
                   method = c("plsr", "plsrda","plslda","plsqda")[1],
                   prior = c("unif", "prop")[1],
                   
                   step = c("nlvtest","permutation","model","prediction")[1],
                   nlv, 
                   modeloutput = c("scores","loadings","coef","vip"), 
                   
                   cvmethod = c("kfolds","loo")[1], 
                   nbrep = 30, 
                   seed = 123, 
                   samplingk = NULL, 
                   nfolds = 10, 
                   npermut = 30, 
                   
                   criterion = c("err","rmse")[1], 
                   selection = c("localmin","globalmin","1std")[1],
                   
                   outputfilename = NULL)
}

\arguments{

\item{import}{If "R", X and Y are in the global environment, and the observation names are in rownames. If "ChemFlow", X and Y are tabulated tables (.txt), and the observation names are in the first column. If "W4M", X and Y are tabulated tables (.txt), and the observation names are in the headers of X, and in the first column of Y.}

\item{X}{Training X-data (\eqn{n, p}). 

\item{Xname}{name of the X-matrix}

\item{Xscaling}{X variable scaling among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{Y}{Training Y-data (\eqn{n, q}) for plsr models, and (\eqn{n, 1}) for plsrda, plslda or plsqda models.}

\item{Yscaling}{Y variable scaling among "none" (mean-centering only), "pareto" (mean-centering and pareto scaling), "sd" (mean-centering and unit variance scaling). If "pareto" or "sd", uncorrected standard deviation is used.}

\item{weights}{Weights (\eqn{n, 1}) to apply to the training observations. Internally, weights are "normalized" to sum to 1. Default to \code{NULL} (weights are set to \eqn{1 / n}).}

\item{newX}{New X-data (\eqn{m, p}) to consider.}

\item{newXname}{name of the newX-matrix}

\item{method}{method to apply among "plsr", "plsrda","plslda","plsqda"}

\item{prior}{for plslda or plsqda models : The prior probabilities of the classes. Possible values are "unif" (default; probabilities are set equal for all the classes) or "prop" (probabilities are set equal to the observed proportions of the classes in \code{y}).}

\item{step}{step of the analysis among "nlvtest" (cross-validation to help determine the optimal number of latent variables), "permutation" (permutation test),"model" (model calculation),"prediction" (prediction of newX-data or X-data if any))}

\item{nlv}{number of latent variables to test if step is "nlvtest"; number of latent variables of the model if step is not "nlvtest".}

\item{modeloutput}{if step is "model": outputs among "scores", "loadings", "coef" (regression coefficients), "vip" (Variable Importance in Projection; the VIP calculation being based on the proportion of Y-variance explained by the components, as proposed by Mehmood et al (2012, 2020).)}

\item{cvmethod}{if step is "nlvtest" or "permutation": "kfolds" for k-folds cross-validation, or "loo" for leave-one-out.}
                   
\item{nbrep}{if step is "nlvtest" and cvmethod is "kfolds": An integer, setting the number of CV repetitions. Default value is 30. Must me set to 1 if cvmethod is "loo"}
                   
\item{seed}{if step is "nlvtest" and cvmethod is "kfolds", or if step is "permutation: a numeric. Seed used for the repeated resampling} 
                   
\item{samplingk}{A vector of length n. The elements are the values of a qualitative variable used for stratified partition creation. If NULL, the first observation is set in the first fold, the second observation in the second fold, etc...}
                   
\item{nfolds}{if cvmethod is "kfolds". An integer, setting the number of partitions to create. Default value is 10.} 
                   
\item{npermut}{if step is "permutation": An integer, setting the number of Y-Block with permutated responses to create. Default value is 30.} 
                   
\item{criterion}{if step is "nlvtest" or "permutation" and method is "plsrda", "plslda" or "plsqda": optimisation criterion among "rmse" and "err" (for classification error rate)))}
                   
\item{selection}{if step is "nlvtest": a character indicating the selection method to use to choose the optimal combination of components, among "localmin","globalmin","1std". If "localmin": the optimal combination corresponds to the first local minimum of the mean CV rmse or error rate. If "globalmin" : the optimal combination corresponds to the minimum mean CV rmse or error rate. If "1std" (one standard error rule) : it corresponds to the first combination after which the mean cross-validated rmse or error rate does not decrease significantly.}
                   
\item{outputfilename}{character: If not NULL, name of the tabular file, in which the function outputs have to be written.)

}

\value{

If step is "nlvtest": table with rmsecv or cross-validated classification error rates. The suggested optimal number of latent variables is indicated by the binary "optimum" variable.

If step is "permutation": table with the dissimilarity between the original and the permutated Y-block, and the rmsecv or cross-validated classification error rates obtained with the permutated Y-block by the model and the given number of latent variables.

If step is "model": tables of scores, loadings, regression coefficients, and vip values, depending of the "modeloutput" parameter.

If step is "prediction": table of predicted scores and predicted classes or values.

}

\examples{

n <- 50 ; p <- 8
Xtrain <- matrix(rnorm(n * p), ncol = p)
colnames(Xtrain) <- paste0("V",1:p)
ytrain <- sample(c(1, 4, 10), size = n, replace = TRUE)

Xtest <- Xtrain[1:5, ] ; ytest <- ytrain[1:5]

resnlvtestplsrda <- plsr_plsda_allsteps(X = Xtrain, Xname = NULL, Xscaling = c("none","pareto","sd")[1], 
                   Y = ytrain, Yscaling = "none", weights = NULL,
                   newX = Xtest, newXname = NULL,
                   
                   method = c("plsr", "plsrda","plslda","plsqda")[2],
                   prior = c("unif", "prop")[1],
                   
                   step = c("nlvtest","permutation","model","prediction")[1],
                   nlv = 5, 
                   modeloutput = c("scores","loadings","coef","vip"), 
                   
                   cvmethod = c("kfolds","loo")[2], 
                   nbrep = 1, 
                   seed = 123, 
                   samplingk = NULL, 
                   nfolds = 10, 
                   npermut = 5, 
                   
                   criterion = c("err","rmse")[1], 
                   selection = c("localmin","globalmin","1std")[1],
                   
                   outputfilename = NULL)
                   
respermutationplsrda <- plsr_plsda_allsteps(X = Xtrain, Xname = NULL, Xscaling = c("none","pareto","sd")[1], 
                   Y = ytrain, Yscaling = "none", weights = NULL,
                   newX = Xtest, newXname = NULL,
                   
                   method = c("plsr", "plsrda","plslda","plsqda")[2],
                   prior = c("unif", "prop")[1],
                   
                   step = c("nlvtest","permutation","model","prediction")[2],
                   nlv = 2, 
                   modeloutput = c("scores","loadings","coef","vip"), 
                   
                   cvmethod = c("kfolds","loo")[2], 
                   nbrep = 1, 
                   seed = 123, 
                   samplingk = NULL, 
                   nfolds = 10, 
                   npermut = 5, 
                   
                   criterion = c("err","rmse")[1], 
                   selection = c("localmin","globalmin","1std")[1],
                   
                   outputfilename = NULL)
                   
plotxy(respermutationplsrda, pch=16)
abline (h = respermutationplsrda[respermutationplsrda[,"permut_dyssimilarity"]==0,"res_permut"])

resmodelplsrda <- plsr_plsda_allsteps(X = Xtrain, Xname = NULL, Xscaling = c("none","pareto","sd")[1], 
                   Y = ytrain, Yscaling = "none", weights = NULL,
                   newX = Xtest, newXname = NULL,
                   
                   method = c("plsr", "plsrda","plslda","plsqda")[2],
                   prior = c("unif", "prop")[1],
                   
                   step = c("nlvtest","permutation","model","prediction")[3],
                   nlv = 2, 
                   modeloutput = c("scores","loadings","coef","vip"), 
                   
                   cvmethod = c("kfolds","loo")[2], 
                   nbrep = 1, 
                   seed = 123, 
                   samplingk = NULL, 
                   nfolds = 10, 
                   npermut = 5, 
                   
                   criterion = c("err","rmse")[1], 
                   selection = c("localmin","globalmin","1std")[1],
                   
                   outputfilename = NULL)
                   

resmodelplsrda$scores
resmodelplsrda$loadings
resmodelplsrda$coef
resmodelplsrda$vip

respredictionplsrda <- plsr_plsda_allsteps(X = Xtrain, Xname = NULL, Xscaling = c("none","pareto","sd")[1], 
                   Y = ytrain, Yscaling = "none", weights = NULL,
                   newX = Xtest, newXname = NULL,
                   
                   method = c("plsr", "plsrda","plslda","plsqda")[2],
                   prior = c("unif", "prop")[1],
                   
                   step = c("nlvtest","permutation","model","prediction")[4],
                   nlv = 2, 
                   modeloutput = c("scores","loadings","coef","vip"), 
                   
                   cvmethod = c("kfolds","loo")[2], 
                   nbrep = 1, 
                   seed = 123, 
                   samplingk = NULL, 
                   nfolds = 10, 
                   npermut = 5, 
                   
                   criterion = c("err","rmse")[1], 
                   selection = c("localmin","globalmin","1std")[1],
                   
                   outputfilename = NULL)
                   
}

\keyword{datagen}