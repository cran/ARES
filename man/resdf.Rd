\name{resdf}
\alias{resdf}
\title{Residual Degrees of Freedom}
\description{
Extract residual degrees of freedom from a GAM model
}
\usage{
resdf(model)
}
\arguments{
  \item{model}{a model fitted by \code{\link{fit.core}}}
}
\details{
This function extracts necessary information from \code{model}. It include both linear and smooth terms.
}
\value{
A scalar with the residual degrees of freedom.
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
resdf(m)
}
\keyword{regression}
