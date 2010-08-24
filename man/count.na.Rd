\name{count.na}
\alias{count.na}
\title{Count Missing}
\description{
Count and return some statistics about missing values
}
\usage{
count.na(var)
}
\arguments{
  \item{var}{a vector}
}
\details{
The function take \code{var} and count the occurrence of missing data.
}
\value{
This function returns a list with
  \item{n.total}{number of observations of \code{var}}
  \item{na}{number of missing observations}
  \item{n.valid}{number of valid observations}
  \item{percent.na}{percentage of missing observations}
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

Ares-Rio Program at \url{http://www.ims.uerj.br/ares-rio}
}
\seealso{\code{\link{is.na}},\code{\link{na.action}}}
\examples{
data(admrio)
setup(admrio,"date")
count.na(pm10)
}
\keyword{ts}
\keyword{univar}
