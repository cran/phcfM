\name{deforestation}

\alias{deforestation}

\title{Binary logistic regression with variable time-intervals between
  observations for modelling deforestation}

\description{The \code{deforestation()} function estimates the parameters of a
  Binary logistic regression model with variable time-intervals between
  observations in a hierarchical Bayesian framework. To estimate the
  posterior distribution of the parameters, an adaptive Metropolis
  algorithm is used. The user supplies data and priors and a sample
  from the posterior distribution is returned as an MCMC object, which
  can be subsequently analyzed with functions provided in the \code{coda}
  package.}
  
\usage{deforestation(formula, interval=1, data, burnin=1000, mcmc=1000,
   thin=1, verbose=1, seed=NA, tune=1, beta.start=NA, mubeta=0,
   Vbeta=1.0E6)}

\arguments{

  \item{formula}{A two-sided linear formula of the form 'y~x1+...+xp'
    describing the model, with the response on the left of a '~'
    operator and the p fixed terms, separated by '+' operators, on the
    right. Response variable y must be 0 or 1 (Bernoulli trials).}

  \item{interval}{A numeric scalar or a vector of length equal to the
    number of observations. \code{interval} specifies the time interval
    between land cover observations. Default to 1.}
  
  \item{data}{A data frame containing the variables of the model.}

  \item{burnin}{The number of burnin iterations for the sampler.}
    
  \item{mcmc}{The number of Gibbs iterations for the sampler. Total
    number of Gibbs iterations is equal to
    \code{burnin+mcmc}. \code{burnin+mcmc} must be divisible by 10 and
    superior or equal to 100 so that the progress bar can be displayed.}
    
  \item{thin}{The thinning interval used in the simulation. The number of
    Gibbs iterations must be divisible by this value.}

  \item{verbose}{A switch (0,1) which determines whether or not the
    progress of the sampler is printed to the screen. Default is 1: a
    progress bar is printed, indicating the step (in \%) reached by the
    Gibbs sampler.}

  \item{seed}{The seed for the random number generator. If NA, the
    Mersenne twister generator is used with default seed 12345; if an
    integer is passed it is used to seed the Mersenne twister.}

  \item{tune}{Metropolis tuning parameter. Must be a positive
    scalar. The tuning parameter is updated during the burning period
    to approach an acceptance rate of 0.44.}
    
  \item{beta.start}{The starting values for the \eqn{\beta}{beta}
    vector. This can either be a scalar or a p-length vector. The
    default value of NA will use the OLS \eqn{\beta}{beta} estimate of
    the corresponding Binary logistic regression. If this is a scalar,
    that value will serve as the starting value for all of the betas.}

  \item{mubeta}{The prior mean of \eqn{\beta}{beta}. This can either be
    a scalar or a p-length vector. If this takes a scalar value, then
    that value will serve as the prior mean for all of the betas. The
    default value of 0 will use a vector of zeros for an uninformative
    prior.}
  
  \item{Vbeta}{The prior variance of \eqn{\beta}{beta}. This can either
    be a scalar or a square p-dimension matrix. If this takes a scalar value,
    then that value times an identity matrix serves as the prior
    variance of beta. Default value of 1.0E6 will use a diagonal matrix
    with very large variance for an uninformative flat prior.}
}

\value{

  \item{mcmc}{An MCMC object that contains the posterior sample. This
    object can be summarized by functions provided by the coda package.}

  \item{deviance}{The posterior mean of the deviance \eqn{D}{D}, with
    \eqn{D=-2\log(\prod_i P(y_i|\theta_i))}{D=-2log(prod_i
      P(y_i|theta_i'))}.}

  \item{tune}{The optimized value for the tuning parameter. This value
  can be used for potential future runs.}
  
}

\details{
  
  The \code{deforestation()} function estimates the parameters of a
  logistic regression model with variable time-intervals between
  observations. The estimation is done in a hierarchical Bayesian
  framework using an adaptive Metropolis algorithm. Function is
  developped in C++ code using the Scythe statistical library (Pemstein
  et al. 2007) to maximize efficiency.
  
  The model takes the following form:
  
  \deqn{y_i \sim \mathcal{B}ernoulli(\theta'_i)}{y_i ~ Bernoulli(theta'_i)}

  With \eqn{\theta_i' = 1-(1-\theta_i)^{I_i}}{theta_i' =
  1-(1-theta_i)^I_i}, where \eqn{I_i} stands for the time-interval for
  observation \eqn{i}. Thus, \eqn{\theta_i}{theta_i} is a probability by
  unit of time (an annual rate for example if the unit of the
  time-interval is in year).

  Using the logit link function denoted \eqn{\phi}{phi}, we set:
  \deqn{\phi(\theta_i) = X_i \beta}{phi(theta_i) = X_i * beta}
  
  By default, we assume a multivariate Normal prior on
  \eqn{\beta}{beta}: \deqn{\beta \sim
  \mathcal{N}(\mu_{\beta},V_{\beta})}{beta ~ N(mu_beta,V_beta)}

  For the Metropolis algorithm, the proposal distribution is centered at
  the current value of \eqn{\beta}{beta} and has variance-covariance
  \eqn{V = T (V_{\beta}^{-1} + C^{-1})^{-1} T }{V = T (V_beta^(-1) +
  C^(-1))^(-1) T}, where \eqn{T}{T} is the diagonal positive definite
  matrix formed from the \code{tune}, \eqn{V_{\beta}}{V_beta} is the
  prior variance, and \eqn{C}{C} is the large sample variance-covariance
  matrix of the MLEs without considering the variation of the
  time-interval between observations. This last calculation is done via
  an initial call to \code{glm}.

}

\references{
  
  Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  
  \emph{Scythe Statistical Library 1.0.} \url{http://scythe.wustl.edu}.
   
  Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
  \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
  \url{http://www-fis.iarc.fr/coda/}.

  Ghislain Vieilledent, Clovis Grinand and Romuald
  Vaudry. 2013. Forecasting deforestation and carbon emissions in
  tropical developing countries facing demographic expansion: a case
  study in Madagascar. \emph{Ecology and Evolution}. DOI:
  10.1002/ece3.550
}

\author{
  Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
}

\seealso{
  \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}
}

\examples{

\dontrun{
#=====================================================================
# Logistic regression with variable time-interval between observations 
#=====================================================================

#= Library
library(phcfM)

#== Generating data

# Random seed
set.seed(1234)

# Constants
nobs <- 3000

# Covariates
X1 <- runif(n=nobs,min=-10,max=10)
X2 <- runif(n=nobs,min=-10,max=10)
X <- cbind(rep(1,nobs),X1,X2)
I <- runif(n=nobs,min=1,max=5) # Time-interval

# Target beta parameters
beta.target <- matrix(c(0.3,0.2,0.1),ncol=1)

# Response
theta <- vector()
theta_prim <- vector()
Y <- vector()
for (n in 1:nobs) {
  theta[n] <- inv.logit(X[n,]\%*\%beta.target)
  theta_prim[n] <- 1-(1-theta[n])^I[n]
  Y[n] <- rbinom(n=1,size=1,prob=theta_prim[n])
}

# Data-set
Data <- as.data.frame(cbind(Y,I,theta_prim,theta,X1,X2))
plot(Data$X1,Data$theta)
plot(Data$X2,Data$theta)

#== Call to deforestation()
model <- deforestation(formula=Y~X1+X2, interval=Data$I, data=Data, burnin=1000, mcmc=1000,
                       thin=1, verbose=1, seed=NA, tune=1, beta.start=NA, mubeta=0,
                       Vbeta=1.0E6)

#== MCMC analysis

# Graphics
plot(model$mcmc)

# Parameter estimates
str(model)
summary(model$mcmc)
model$deviance
model$tune

## We obtain good parameter estimates 
##
##                    Mean       SD  Naive SE Time-series SE
## beta.(Intercept) 0.3362 0.054534 0.0017245       0.006589
## beta.X1          0.2027 0.009345 0.0002955       0.001154
## beta.X2          0.1037 0.007269 0.0002299       0.000838

#== GLM resolution if time-interval is not taken into account

model.glm <- glm(Y~X1+X2,data=Data,family="binomial")
summary(model.glm)

## In this case, the parameter estimates are biased
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) 2.000761   0.072069   27.76   <2e-16 ***
## X1          0.247514   0.012060   20.52   <2e-16 ***
## X2          0.123137   0.009803   12.56   <2e-16 ***
}
}

\keyword{models}
\keyword{hierarchical models}
\keyword{glm}
\keyword{logit}
\keyword{Binary}
\keyword{MCMC}
\keyword{bayesian}