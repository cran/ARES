\name{periodogram.test}
\alias{periodogram.test}
\title{Periodogram Test}
\description{
Uniformity test for the periodogram
}
\usage{
periodogram.test(object, plot = TRUE)
}
\arguments{
  \item{object}{an object created by \code{\link{periodogram}}}
  \item{plot}{a logical indicating whether the cumulative probability should be plotted}
}
\details{
A Kolmogorov-Smirnov (KS) test is performed to evaluate de uniformity of the \eqn{I\left(\omega\right)}. 

See \code{\link{ks.test}} for details on KS tests.
}
\value{
This function does not return a value.
}
\references{
Box, G., Jenkins, G., Reinsel, G. (1994) \emph{Time Series Analysis : Forecasting and Control}. 3 ed., Prentice Hall.

Diggle, P. J. (1990) \emph{Time Series: A Biostatistical Introduction}. Oxford University Press.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

Ares-Rio Program at \url{http://www.ims.uerj.br/ares-rio}
}

\seealso{\code{\link{fit.core}},\code{\link{get.residuals}},\code{\link{periodogram}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit.core(f)
p <-periodogram(m,test=FALSE)
periodogram.test(p)
}
\keyword{ts}
\keyword{regression}
