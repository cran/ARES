\name{plot.cook}
\alias{plot.cook}
\title{Plot Cook}
\description{
Plot Cook's distance
}
\usage{
plot.cook(x, type = "deviance", line = 0.1, new = TRUE, ...)
}
\arguments{
  \item{x}{a model fitted by \code{\link{fit.core}}}
  \item{type}{a quoted string indicating the type of residuals to extract. Default is \code{"deviance"}. See \code{\link{get.residuals}}}
  \item{line}{a value indicating where to plot a horizontal line on the y-axis. Set \code{NULL} for none}
  \item{new}{if \code{TRUE} a new graph window is opened}
  \item{\dots}{further options for \code{\link{plot}}}
}
\details{
This function extracts the residuals of \code{x}, computes de Cook's distance, and generates a plot for each observation. It is a very useful tool for influence analysis.
}
\value{
This function does not return a value.
}
\references{
McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
 
Cook, R. D., Weisberg, S. (1982) \emph{Residuals and Influence in Regression}. Chapman and Hall.

Atkinson, A.C. (1985) \emph{Plots, Transformations and Regression}. Oxford University Press.

Pierce, D. A., Schafer, D. W. (1986) Residuals in generalized linear models. \emph{Journal of the American Statistical Association}, \bold{81(396)},977--986.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

Ares-Rio Program at \url{http://www.ims.uerj.br/ares-rio}
}

\seealso{\code{\link{fit.core}},\code{\link{get.residuals}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit.core(f)
plot.cook(m)
}
\keyword{hplot}
\keyword{regression}
