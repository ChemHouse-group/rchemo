\name{Zhang2023}
\alias{Zhang2023}
\docType{data}
\encoding{latin1}

\title{Zhang2023}

\description{

A selection of metabolomics reduced datasets from the article of Zhang et al.(2023)

}

\usage{data(Zhang2023)}

\format{
  A list with 3 components: \code{GCTOF}, \code{HILICPOS}, \code{metadata}.
  
  \describe{
    \item{\code{GCTOF}}{A matrix whose rows are 68 mice and columns are 108 metabolomic variables obtained by GC-MS.}
    \item{\code{HILICPOS}}{A matrix whose rows are 68 mice and columns are 133 metabolomic variables obtained by LC-MS.}
    \item{\code{metadata}}{A matrix whose rows are 68 mice and columns are the genotype group and the gender.}

  }
}

\source{

Zhang, Y.; Barupal, D.K.; Fan, S.; Gao, B.; Zhu, C.; Flenniken, A.M.; McKerlie, C.; Nutter, L.M.J.; Lloyd, K.C.K.; Fiehn, O. 
Sexual Dimorphism of the Mouse Plasma Metabolome Is Associated with Phenotypes of 30 Gene Knockout Lines. 
Metabolites 2023, 13, 947. https://doi.org/10.3390/metabo13080947

}

\examples{

data(Zhang2023)

x <- Zhang2023$GCTOF
head(x)

y <- Zhang2023$metadata
head(y)

}


\keyword{datasets}