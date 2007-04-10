\name{decomposeSurv}
\alias{decomposeSurv}
\title{Decompose A Survival Formula}

\description{
Decompose a \code{Surv} formula which can have time-dependent effects and return the model matrix and the time function matrix. 
}

\usage{
decomposeSurv(formula, data, sort)
}

\arguments{
\item{formula}{A formula with \code{Surv} on the left-handed side and an arbitrary part on the right-handed side, including possibly time-dependent effects, e.g. \code{Surv(time, event) ~ x1 * x2 + I(x3^2) + x1 * log(time)}.}
\item{data}{The \code{data.frame} with the columns needed to extract the formula.}
\item{sort}{If TRUE, the returned model matrizes are sorted by stoptime and event.}
}

\details{
This function is mainly a helper function for \code{coxphf}, \code{coxphftest} and \code{coxphfplot}.
}

\value{
\item{NTDE}{The number of time dependent effects.}
\item{fac}{The factor matrix, see the value of \code{terms}.}
\item{resp}{N x 3 - response matrix, see \code{model.extract} with parameter \code{response}.}
\item{mm1}{Model matrix without time effects.}
\item{timedata}{Matrix with time functions as columns.}
\item{timeind}{Indicator of time-dependent effect: the i-th column of \code{timedata} interacts with the timeind[i]-th column of \code{mm1}.}
\item{covnames}{The names of the covariates.}
\item{ind}{An indicator vector. Elements with FALSE mean that these \code{covnames} are not part of the \code{formula}, but are needed to span a coherent model frame.}
}

\references{
Heinze, G. and Schemper, M. (2001). A Solution to the Problem of Monotone Likelihood in Cox Regression. Biometrics 57/1, 114-119.

Heinze, G. (1999). Technical Report 10/1999: The application of Firth's procedure to Cox and logistic regression. Section of Clinical Biometrics, Department of Medical Computer Sciences, University of Vienna, Vienna.

Ploner, M. and Heinze, G. (2001). Technical Report 1/2001: An SPLUS library to perform Cox regression without convergence problems. Section of Clinical Biometrics, Department of Medical Computer Sciences, University of Vienna, Vienna.

Heinze, G. and Ploner, M. (2002). SAS and SPLUS programs to perform Cox regression without convergence problems. Computer Methods and Programs in Biomedicine 67, 217-223.
}

\author{Meinhard Ploner, Georg Heinze}


\examples{
## prepares matrices for a model with time-dependent effects 
test3 <- data.frame(list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                         stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
                         event=c(1, 0, 0, 1, 0, 1, 1, 0, 0, 0),
                         x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) ))
                         
decomposeSurv(Surv(start,stop,event) ~ x + x * log(stop), data=test3)
}

\keyword{survival}
