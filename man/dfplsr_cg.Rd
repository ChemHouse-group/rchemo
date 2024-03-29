\name{dfplsr_cg}
\alias{dfplsr_cg}
\alias{dfplsr_cov}
\alias{dfplsr_div}
\encoding{latin1}

\title{Degrees of freedom of Univariate PLSR Models}

\description{

Computation of the model complexity \eqn{df} (number of degrees of freedom) of univariate PLSR models (with intercept). See Lesnoff et al. 2021 for an illustration.

(1) Estimation from the CGLSR algorithm (Hansen, 1998).

- \code{dfplsr_cov}

(2) Monte Carlo estimation (Ye, 1998 and Efron, 2004). Details in relation with the functions are given in Lesnoff et al. 2021.

- \code{dfplsr_cov}: The covariances are computed  by parametric bootstrap (Efron, 2004, Eq. 2.16). The residual variance \eqn{sigma^2} is estimated from a low-biased model.

- \code{dfplsr_div}: The divergencies \eqn{dy_fit/dy} are computed by perturbation analysis(Ye, 1998 and Efron, 2004). This is a Stein unbiased risk estimation (SURE) of \eqn{df}.

}

\usage{

dfplsr_cg(X, y, nlv, reorth = TRUE)

dfplsr_cov(
    X, y, nlv, algo = NULL,
    maxlv = 50, B = 30, print = FALSE, ...)

dfplsr_div(
    X, y, nlv, algo = NULL,
    eps = 1e-2, B = 30, print = FALSE, ...) 

}

\arguments{

\item{X}{A \eqn{n x p} matrix or data frame of training observations.}
  
\item{y}{A vector of length \eqn{n} of training responses.}
  
\item{nlv}{The maximal number of latent variables (LVs) to consider in the model.}
  
\item{reorth}{For \code{dfplsr_cg}: Logical. If \code{TRUE}, a Gram-Schmidt reorthogonalization of the normal equation residual vectors is done.}
  
\item{algo}{a PLS algorithm. Default to  \code{NULL} (\code{\link{plskern}} is used).}

\item{maxlv}{For \code{dfplsr_cov}: dDmension of the PLSR model (nb. LVs) used for parametric bootstrap.}

\item{eps}{For \code{dfplsr_div}: The \eqn{epsilon} quantity used for scaling the perturbation analysis.}

\item{B}{For \code{dfplsr_cov}: Number of bootstrap replications. For \code{dfplsr_div}: number of observations in the data receiving perturbation (the maximum is \eqn{n}).}

\item{print}{Logical. If \code{TRUE}, fitting information are printed.}

\item{...}{Optionnal arguments to pass in the function defined in \code{algo}.}

}

\details{

Missing values are not allowed.

The example below reproduces the numerical illustration given by Kramer & Sugiyama 2011 on the Ozone data (Fig. 1, center).
The \eqn{pls.model} function from the R package "plsdof" v0.2-9 (Kramer & Braun 2019) is used for \eqn{df} calculations (\eqn{df.kramer}), and automatically scales the X matrix before PLS. The example scales also X for consistency when using the other functions.

For the Monte Carlo estimations, \code{B} Should be increased for more stability


}

\value{

A list of outputs :

\item{df}{vector with the model complexity for the models with \eqn{a = 0, 1, ..., nlv} components.}

\item{cov}{For \code{dfplsr_cov}: vector with covariances, computed  by parametric bootstrap.}

}

\references{

Efron, B., 2004. The Estimation of Prediction Error. Journal of the American Statistical Association 99,
619-632. https://doi.org/10.1198/016214504000000692

Hastie, T., Tibshirani, R.J., 1990. Generalized Additive Models, Monographs on statistics and applied
probablity. Chapman and Hall/CRC, New York, USA.

Hastie, T., Tibshirani, R., Friedman, J., 2009. The elements of statistical learning: data mining,
inference, and prediction, 2nd ed. Springer, NewYork.

Hastie, T., Tibshirani, R., Wainwright, M., 2015. Statistical Learning with Sparsity: The Lasso and
Generalizations. CRC Press

Kramer, N., Braun, M.L., 2007. Kernelizing PLS, degrees of freedom, and efficient model selection, in: Proceedings of the 24th International Conference on Machine Learning, ICML 07. Association for Computing Machinery, New York, NY, USA, pp. 441-448. https://doi.org/10.1145/1273496.1273552

Kramer, N., Sugiyama, M., 2011. The Degrees of Freedom of Partial Least Squares Regression. Journal of the American Statistical Association 106, 697-705. https://doi.org/10.1198/jasa.2011.tm10107

Kramer, N., Braun, M. L. 2019. plsdof: Degrees of Freedom and Statistical Inference for Partial Least Squares Regression. R package version 0.2-9. https://cran.r-project.org

Lesnoff, M., Roger, J.M., Rutledge, D.N., 2021. Monte Carlo methods for estimating Mallow's Cp and AIC criteria for PLSR models. Illustration on agronomic spectroscopic NIR data. Journal of Chemometrics, 35(10), e3369. https://doi.org/10.1002/cem.3369

Stein, C.M., 1981. Estimation of the Mean of a Multivariate Normal
Distribution. The Annals of Statistics 9, 1135-1151.

Ye, J., 1998. On Measuring and Correcting the Effects of Data Mining and Model Selection. Journal of
the American Statistical Association 93, 120-131. https://doi.org/10.1080/01621459.1998.10474094

Zou, H., Hastie, T., Tibshirani, R., 2007. On the degrees of freedom of the lasso. The Annals of
Statistics 35, 2173-2192. https://doi.org/10.1214/009053607000000127

}

\examples{

## EXAMPLE 1

data(ozone)

z <- ozone$X
u <- which(!is.na(rowSums(z)))
X <- z[u, -4]
y <- z[u, 4]
dim(X)

Xs <- scale(X)

nlv <- 12
res <- dfplsr_cg(Xs, y, nlv = nlv)

df.kramer <- c(1.000000, 3.712373, 6.456417, 11.633565, 12.156760, 11.715101, 12.349716,
  12.192682, 13.000000, 13.000000, 13.000000, 13.000000, 13.000000)

znlv <- 0:nlv
plot(znlv, res$df, type = "l", col = "red",
     ylim = c(0, 15),
     xlab = "Nb components", ylab = "df")
lines(znlv, znlv + 1, col = "grey40")
points(znlv, df.kramer, pch = 16)
abline(h = 1, lty = 2, col = "grey")
legend("bottomright", legend=c("dfplsr_cg","Naive df","df.kramer"), col=c("red","grey40","black"),
lty=c(1,1,0), pch=c(NA,NA,16), bty="n")

## EXAMPLE 2

data(ozone)

z <- ozone$X
u <- which(!is.na(rowSums(z)))
X <- z[u, -4]
y <- z[u, 4]
dim(X)

Xs <- scale(X)

nlv <- 12
B <- 50 
u <- dfplsr_cov(Xs, y, nlv = nlv, B = B)
v <- dfplsr_div(Xs, y, nlv = nlv, B = B)

df.kramer <- c(1.000000, 3.712373, 6.456417, 11.633565, 12.156760, 11.715101, 12.349716,
  12.192682, 13.000000, 13.000000, 13.000000, 13.000000, 13.000000)

znlv <- 0:nlv
plot(znlv, u$df, type = "l", col = "red",
     ylim = c(0, 15),
     xlab = "Nb components", ylab = "df")
lines(znlv, v$df, col = "blue")                 
lines(znlv, znlv + 1, col = "grey40")
points(znlv, df.kramer, pch = 16)
abline(h = 1, lty = 2, col = "grey")
legend("bottomright", legend=c("dfplsr_cov","dfplsr_div","Naive df","df.kramer"), 
col=c("blue","red","grey40","black"),
lty=c(1,1,1,0), pch=c(NA,NA,NA,16), bty="n")

}

\keyword{datagen}