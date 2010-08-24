\name{save.plot}
\alias{save.plot}
\title{Save Plot}
\description{
Conveniently save the active plot window to some known graph file format
}
\usage{
save.plot(file, width = 520, height = 480)
}
\arguments{
  \item{file}{quoted string with the file name}
  \item{width}{width in pixels}
  \item{height}{height in pixels}
}
\details{
This function tries to guess the filter to use from the file name extension. Supported graph formats are Windows meta-file (wmf), portable bitmap (png), jpeg (jpg), portable document format (pdf), encapsulated post script (eps), and scalable vector graphics (svg). If \code{off} is set to \code{TRUE}, the active window device is closed after the graph is saved.

Alternatively, one can use the quite new \R native \code{\link{savePlot}} that essentially the same job.

This function has been rewritten and GhostScript is not necessary any more.
}
\value{ 
This function does not return a value.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

Ares-Rio Program at \url{http://www.ims.uerj.br/ares-rio}
}
\seealso{\code{\link{dev.copy}},\code{\link{dev.copy2eps}}}
\examples{
data(admrio)
setup(admrio,"date")
plot.event(resp5)
save.plot("plot_resp5.png")
}
\keyword{device}
