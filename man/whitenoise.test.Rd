\name{whitenoise.test}
\alias{whitenoise.test}
\title{White Noise Test}
\description{
Run a series of white noise tests
}
\usage{
whitenoise.test(object, type = "deviance", k = 25, ...)
}
\arguments{
  \item{object}{a model fitted by \code{\link{fit.core}} or a vector of residuals}
  \item{type}{a quoted string indicating the type of residuals to extract. Default is \code{"deviance"}. See \code{\link{get.residuals}}}
  \item{k}{an integer indicating the number of lags to compute the partial autocorrelation coefficients}
  \item{\dots}{further options for \code{\link{ljungbox.test}}}
}
\details{
This function submits the extracted residuals to a series of white noise tests. 

For normality evaluation the robust tests Bownman-Shenton and Jarque-Bera tests are implemented. Besides, Kolmogorov-Smirnov and Shapiro tests are performed using \code{\link{ks.test}} and \code{\link{shapiro.test}} respectively.

Serial dependence of the residuals is evaluated using \code{\link{ljungbox.test}}.

Heteroscedasticity tests will be implemented soon.

Results are printed out on the console.
}
\value{
This function does not return a value.
}
\references{
Bowman, K.O. and Shenton L.R. (1975) Omnibus test contours for departures from normality based on b1 and b2. \emph{Biometrika} \bold{62}, 243--250.

Jarque, C. M. and Bera, A. K. (1987) A Test for Normality of Observations and Regression Residuals. \emph{International Statistical Review}, \bold{55(2)}, 163--172.

Spanos, A. (1999) \emph{Probability Theory and Statistical Inference}. Cambridge.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

Ares-Rio Program at \url{http://www.ims.uerj.br/ares-rio}
}

\seealso{\code{\link{get.residuals}},\code{\link{acf}}}
\examples{
library(ares)
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit.core(f)
r <- get.residuals(m)
whitenoise.test(r)
}
\keyword{ts}
\keyword{htest}
