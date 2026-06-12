\name{Zhang2023}
\alias{Zhang2023}
\docType{data}
\encoding{latin1}

\title{Zhang2023}

\description{

A selection of reduced non-targeted metabolomics datasets from the article of Zhang et al.(2023)

}

\usage{data(Zhang2023)}

\format{
  A list with 6 components: \code{GCTOF}, \code{HILICNEG}, \code{HILICPOS}, \code{metadata}.
  
  \describe{
    \item{\code{CHSNEG}}{A matrix whose rows are 68 mice and columns are 163 metabolomic variables obtained by LC-MS.}
    \item{\code{CHSPOS}}{A matrix whose rows are 68 mice and columns are 288 metabolomic variables obtained by LC-MS.}
    \item{\code{GCTOF}}{A matrix whose rows are 68 mice and columns are 108 metabolomic variables obtained by GC-MS.}
    \item{\code{HILICNEG}}{A matrix whose rows are 68 mice and columns are 44 metabolomic variables obtained by LC-MS.}
    \item{\code{HILICPOS}}{A matrix whose rows are 68 mice and columns are 133 metabolomic variables obtained by LC-MS.}
    \item{\code{metadata}}{A matrix whose rows are 68 mice and columns are the genotype group (Mutant/Wild) and the gender (Male/Female).}

  }
}

\details{

The 6 samples from 5 mutant groups (Dhfr, Gnpda1, Plk1, Sra1, Ulk3) and the 40 controls were retained, with the exception of 2 animals wM_035 (wild Male) et mM_102 (mutant Male) that had missing values in non-targeted metabolomics. 
Three HILIC-NEG (Ser_Asn, Théophylline, Val_Asp) and 2 HILIC-POS variables ("Ser-His", "Thr-Arg") were then removed due to missing or infinite values.

}

\source{

Zhang, Y.; Barupal, D.K.; Fan, S.; Gao, B.; Zhu, C.; Flenniken, A.M.; McKerlie, C.; Nutter, L.M.J.; Lloyd, K.C.K.; Fiehn, O. 
Sexual Dimorphism of the Mouse Plasma Metabolome Is Associated with Phenotypes of 30 Gene Knockout Lines. 
Metabolites 2023, 13, 947. https://doi.org/10.3390/metabo13080947

}

\examples{

data(Zhang2023)

X <- Zhang2023$GCTOF
head(X)

Y <- Zhang2023$metadata
head(Y)

}


\keyword{datasets}