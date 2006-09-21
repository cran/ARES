\name{aresCalc}
\alias{aresCalc}
\alias{aftbino}
\alias{aftbino.mbe}
\alias{batbino}
\alias{dmtbino}
\alias{dtbino}
\alias{emtbino}
\alias{fmtbino}
\alias{gtbino}
\alias{mpair}
\alias{mtbino}
\alias{ntbino}
\alias{omtbino}
\alias{removetsp}
\alias{rpldtbino}
\alias{sltbino}
\alias{vtbino}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Allelic richness estimation, with extrapolation beyond the sample size }
\description{
  Generates an allelic richness accumulation curve. This curve shows the expected 
  number of unique alleles in a population when taking a sample of individuals,
  ranging from one till \code{maxsize} individuals. 
  It takes a binary data matrix as input, showing the presence of alleles per 
  individual, and gives an accumulation curve (mean with 95\% confidence bounds) 
  back. The total number of individuals present in the input data can be smaller 
  then \code{maxsize}.
}

\usage{
aresCalc(rawdata, bootsize = NULL, maxsize = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawdata}{ a binary data matrix with columns representing different 
  individuals and rows representing different alleles }
  \item{bootsize}{ the number of bootstrap runs to be made, for calculating confidence bounds (typically between 100 and 1000) }
  \item{maxsize}{ the number of individuals to where the accumulation curve has to be calculated}
}

\value{
The function returns a table of 3 columns and \code{maxsize} rows. The first column
holds the mean (named \code{esti}), the second the lower 97.5 \% conf. bound 
(named \code{lwr}) and the third the upper 97.5\% conf. bound (named \code{upr})
The output can be visualised by the function \code{aresPlot}.
}

\references{ 
Colwell RK, Mao CX, Chang J (2004) 
Interpolating, extrapolating and comparing incidence-based species accumulation curves. 
Ecology, 85, 2717-2727.

Mao CX, Colwell RK, Chang J (2006) Estimating the species accumulation curve using mixtures. 
Biometrics, 61, 433-441.

van Loon EE, Cleary DFR, Fauvelot C (2006) 
ARES: software to compare allelic richness between uneven samples.
Molecular Ecology Notes, (in review)
}

\author{ Emiel van Loon, based on code by Chang Xuan Mao }

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

