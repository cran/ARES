\name{aresPlot}
\alias{aresPlot}
\title{ Plotting an allelic richness accumulation curve }
\description{
  Plots an allelic richness accumulation curve, on the basis of the output
  from \code{aresCalc}.
}

\synopsis{
aresPlot(output_aresCalc, T = NULL)
}
\usage{
aresPlot(output_aresCalc, T = "mytitle")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output_aresCalc}{the output matrix from aresCalc, a dataframe with 
  three columns}
  \item{T}{the plot title}
}

\value{
The function returns a plot of the estimated allelic richness with 95\% 
confidence bounds.
}

\references{ 
van Loon EE, Cleary DFR, Fauvelot C (2006) 
ARES: software to compare allelic richness between uneven samples.
Molecular Ecology Notes, (in review)
}

\author{ Emiel van Loon }

\examples{

# allelic richness of the butterfly Drupadia theda (Felder),  
# at a study site B1 in the year 1998 compared to that in 2000

data(butterfly_borneo)

b198 <- aresCalc(bb198,bootsize=4,maxsize=60)
aresPlot(b198, T="1998, B1")

b100 <- aresCalc(bb100,bootsize=4,maxsize=60)
aresPlot(b100, T="2000, B1")

# note!!! these examples use a very small bootstrap size (for speed of display).
#  A more realistic size is bootsize=200.

}
\keyword{ models }
