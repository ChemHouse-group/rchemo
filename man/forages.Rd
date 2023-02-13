\name{forages}
\alias{forages}
\docType{data}
\encoding{latin1}

\title{forages}

\description{
A NIRS dataset (pre-processed absorbance) describing the class membership of forages. Spectra were recorded from 1100 to 2498 nm at 2 nm intervals.
}

\usage{data(forages)}

\format{
  A list with 4 components: \code{Xtrain}, \code{ytrain}, \code{Xtest}, \code{ytest}.
  
  \describe{
  
For the reference (calibration) data:

    \item{\code{Xtrain}}{A matrix whose rows are the pre-processed NIR absorbance spectra (= log10(1 / Reflectance)).}

    \item{\code{ytrain}}{A vector of the response variable (class membership).}
    
For the test data:

    \item{\code{Xtest}}{A matrix whose rows are the pre-processed NIR absorbance spectra (= log10(1 / Reflectance)).}

    \item{\code{ytest}}{A vector of the response variable (class membership).}
    
  }
}

\examples{

data(forages)
str(forages)

}


\keyword{datasets}
