\name{pgps}
\alias{pgps}
\title{Generalized Pearson's Statistic}
\description{
Compute Generalized Pearson's statistic for a Poisson GAM model
}
\usage{
pgps(model)
}
\arguments{
  \item{model}{a model fitted by \code{\link{fit.core}}}
}
\details{
This function extracts necessary information from \code{model} and compute \eqn{\sum_{t=1}^{n}\left(y_{t}-E\left(y_{t}\right)\right)^{2}/E\left(y_{t}\right)}.
}
\value{
A scalar with the statistic.
}
\references{
McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

Ares-Rio Program at \url{http://www.ims.uerj.br/ares-rio}
}

\seealso{\code{\link{fit.core}},\code{\link{gam}},\code{\link{glm}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit.core(f)
pgps(m)
}
\keyword{regression}
