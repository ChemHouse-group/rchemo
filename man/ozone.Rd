\name{ozone}
\alias{ozone}
\docType{data}
\encoding{latin1}

\title{ozone}

\description{

Los Angeles ozone pollution data in 1976
(sources: Breiman & Friedman 1985,  Leisch &  Dimitriadou 2020).

}

\usage{data(ozone)}

\format{
  A list with 1 component: the matrix \code{X} with 366 observations, 13 variables. The variable to predict is V4.
  \describe{
    \item{\code{V1}}{Month: 1 = January, ..., 12 = December}
    \item{\code{V2}}{Day of month}
    \item{\code{V3}}{Day of week: 1 = Monday, ..., 7 = Sunday}
    \item{\code{V4}}{Daily maximum one-hour-average ozone reading}
    \item{\code{V5}}{500 millibar pressure height (m) measured at Vandenberg AFB}
    \item{\code{V6}}{Wind speed (mph) at Los Angeles International Airport (LAX)}
    \item{\code{V7}}{Humidity (\%) at LAX}
    \item{\code{V8}}{Temperature (degrees F) measured at Sandburg, CA}
    \item{\code{V9}}{Temperature (degrees F) measured at El Monte, CA}
    \item{\code{V10}}{Inversion base height (feet) at LAX}
    \item{\code{V11}}{Pressure gradient (mm Hg) from LAX to Daggett, CA}
    \item{\code{V12}}{Inversion base temperature (degrees F) at LAX}
    \item{\code{V13}}{Visibility (miles) measured at LAX}
  }
}

\source{

Breiman L., Friedman J.H. 1985. Estimating optimal transformations for multiple regression and correlation, JASA, 80, pp. 580-598.

Leisch, F. and Dimitriadou, E. (2010). mlbench: Machine Learning Benchmark Problems.
R package version 1.1-6. https://cran.r-project.org/.

}

\examples{

data(ozone)

z <- ozone$X
head(z)

plotxna(z)

}


\keyword{datasets}