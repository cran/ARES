\name{desc.data}
\alias{desc.data}
\title{Describe Dataset}
\description{
Report some information on a dataset 
}
\usage{
desc.data(dataset = NULL)
}
\arguments{
  \item{dataset}{an object of class \code{data.frame}}
}
\details{
This function outputs information on variables index, names, class, number of missing values, and labels. If \code{dataset} is omitted, the function tries to find \code{.ares.active.dataset}.
}
\value{
A data frame containing columns with the items described in Details.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

Ares-Rio Program at \url{http://www.ims.uerj.br/ares-rio}
}
\seealso{\code{\link{names}},\code{\link{class}}}
\examples{
data(admrio)
setup(admrio,"date")
## specific dataset

## .ares.active.data.set
d <- desc.data()
}
\keyword{data}
\keyword{datagen}

