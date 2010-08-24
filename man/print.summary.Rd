\name{print.summary}
\alias{print.summary}
\title{Print Summary}
\description{
Print model information}
\usage{
print.summary(x, digits = getOption("digits"), ...)
}
\arguments{
  \item{x}{a model fitted by \code{\link{fit.core}}}
  \item{digits}{an integer indicating the decimal places. If not supplied, it is taken from \code{\link{options}}}
  \item{\dots}{further options for \code{\link{summary.glm}} or \code{\link{summary.gam}}}
}
\details{
This function prints out both \code{\link{summary.glm}} and \code{\link{summary.glm}} outputs added of extra relevant information.
}
\value{
The function invisibly returns a list with the following objects
  \item{summary.glm}{summary information of the linear part of the model. See \code{\link{summary.glm}}}
  \item{summary.gam}{summary information of the non-linear part of the model. See \code{\link{summary.gam}}}
  \item{dispersion}{estimate of the dispersion parameter}
  \item{pearson}{estimate of the generalized Pearson's statistic}
  \item{residuals.df}{residuals degrees of freedom}
  \item{deviance}{estimate of the deviance}
}
\references{
McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

Ares-Rio Program at \url{http://www.ims.uerj.br/ares-rio}
}

\seealso{\code{\link{gam}},\code{\link{glm}},\code{\link{summary.gam}},\code{\link{summary.glm}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit.core(f)
print.summary(m)
}
\keyword{print}
\keyword{regression}
