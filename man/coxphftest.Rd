% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxphftest.R
\name{coxphftest}
\alias{coxphftest}
\title{Penalized Likelihood Ratio Test in Cox Regression}
\usage{
coxphftest(
  formula,
  data,
  test = ~.,
  values,
  maxit = 50,
  maxhs = 5,
  epsilon = 1e-06,
  maxstep = 0.5,
  firth = TRUE,
  adapt = NULL,
  penalty = 0.5
)
}
\arguments{
\item{formula}{a formula object, with the response on the left of the  operator, and the model terms on the right. The response must be a survival object as returned by the 'Surv' function.}

\item{data}{a data.frame in which to interpret the variables named in the 'formula' argument.}

\item{test}{righthand formula of parameters to test (e.g. \code{~ B + D}). As default the null hypothesis that all parameters are 0 is tested.}

\item{values}{null hypothesis values, default values are 0. For testing the hypothesis H0: B1=1 and B4=2 and B5=0, specify \code{test= ~ B1 + B4 + B5} and \code{values=c(1, 2, 0)}.}

\item{maxit}{maximum number of iterations (default value is 50)}

\item{maxhs}{maximum number of step-halvings per iterations (default value is 5). The increments of the parameter vector in one Newton-Rhaphson iteration step are halved, unless the new likelihood is greater than the old one, maximally doing \code{maxhs} halvings.}

\item{epsilon}{specifies the maximum allowed change in penalized log likelihood todeclare convergence. Default value is 0.0001.}

\item{maxstep}{specifies the maximum change of (standardized) parameter values allowed in one iteration. Default value is 2.5.}

\item{firth}{use of Firth's penalized maximum likelihood (\code{firth=TRUE}, default) or the standard maximum likelihood method (\code{firth=FALSE}) for fitting the Cox model.}

\item{adapt}{optional: specifies a vector of 1s and 0s, where 0 means that the corresponding parameter is fixed at 0, while 1 enables 
parameter estimation for that parameter. The length of adapt must be equal to the number of parameters to be estimated.}

\item{penalty}{strength of Firth-type penalty. Defaults to 0.5.}
}
\value{
\item{testcov}{the names of the tested model terms}
\item{loglik}{the restricted and unrestricted maximized (penalized) log likelihood}
\item{df}{the number of degrees of freedom related to the test}
\item{prob}{the p-value}
\item{call}{the function call}
\item{method}{the estimation method (penalized ML or ML)}
}
\description{
Performs a penalized likelihood ratio test for hypotheses within a Cox regression analysis using Firth's penalized likelihood.
}
\details{
This function performs a penalized likelihood ratio test on some (or all) selected parameters.  
It can be used to test contrasts of parameters, or factors that are coded in dummy variables. 
The resulting object is of the class coxphftest and includes the information printed by the proper print method.
}
\examples{
library(survival)
testdata <- data.frame(list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
x1    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0),
x2    =c(0, 1, 1, 1, 0, 0, 1, 0, 1, 0),
x3    =c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0)))

summary( coxphf( formula=Surv(start, stop, event) ~ x1+x2+x3, data=testdata))

# testing H0: x1=0, x2=0

coxphftest( formula=Surv(start, stop, event) ~ x1+x2+x3, test=~x1+x2,  data=testdata)



# How to test total effect of a variable with time-dependent effect

# NOT RUN (works)
#fitt<- coxphf( formula=Surv(start, stop, event) ~ x1+x2+x3*stop, data=testdata, pl=FALSE)

#test <- coxphf(formula=Surv(start, stop, event) ~ x1+x2+x3*stop, data=testdata, adapt=c(1,1,0,0))

# PLR p-value for x3 + x3:stop
#pchisq((fitt$loglik[2]-test$loglik[2])*2, 2, lower.tail=FALSE)

#NOT RUN (does not work)
#test <- coxphf(formula=Surv(start, stop, event) ~ x1+x2+x3*stop, data=testdata, test=~x3+stop:x3)


}
\references{
Firth D (1993). Bias reduction of maximum likelihood estimates. \emph{Biometrika} 80:27--38.

Heinze G and Schemper M (2001). A Solution to the Problem of Monotone Likelihood in Cox Regression. \emph{Biometrics} 57(1):114--119. 

Heinze G (1999). Technical Report 10/1999: The application of Firth's procedure to Cox and logistic regression. Section of Clinical Biometrics, Department of Medical Computer Sciences, University of Vienna, Vienna.
}
