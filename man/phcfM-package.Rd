\name{phcfM-package}
\alias{phcfM-package}
\alias{phcfM}
\docType{package}
\title{
  Modelling anthropogenic deforestation in tropical forests
}

\description{\code{phcfM} is an R package for modelling anthropogenic
deforestation. It was initially developed to obtain REDD baseline
scenarios of deforestation for the "programme holistique de conservation
des forets a Madagascar" (from which the package is named after). It
includes two main functions: (i) \code{demography()}, to model the
population growth with time using population census data and Gaussian
mixed linear regressions and (ii) \code{deforestation()}, to model the
deforestation process using land-cover change data and Binary logistic
regressions with variable time-intervals between land-cover
observations. Parameter inference is done in a hierarchical Bayesian
framework. Markov chains Monte Carlo (MCMC) are coded in C++ using the
Scythe statistical library to maximize computation efficiency.}

\details{
\tabular{ll}{
  Package: \tab phcfM\cr
  Type: \tab Package\cr
  Version: \tab 1.1\cr
  Date: \tab 2012-06-13\cr
  License: \tab GPL-3 \cr
  LazyLoad: \tab yes\cr
}
}

\author{
  Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
}

\keyword{tropical forest}
\keyword{deforestation}
\keyword{anthropogenic pressure}
