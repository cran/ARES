\name{ARES-package}
\alias{ARES-package}
\alias{ARES}
\docType{package}
\title{
Allelic richness estimation, with extrapolation beyond the sample size
}
\description{
  Allows to easily generates and plot an allelic richness accumulation curve. 
  This curve shows the expected number of unique alleles in a population when 
  taking a sample of individuals, ranging from one to \code{maxsize} individuals. 
  It takes a binary data matrix as input, showing the presence of alleles per 
  individual, and gives an accumulation curve (mean with 95\% confidence bounds) 
  back. The total number of individuals present in the input data can be smaller 
  then \code{maxsize}.
  The package contains two functions that can be called by users: \code{aresCalc} and \code{aresPlot}

}
\details{
\tabular{ll}{
Package: \tab ARES\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-09-21\cr
License: \tab GPL 2 or newer\cr
}
}


\author{
Emiel van Loon, based on code by Chang Xuan Mao
Maintainer: Emiel van Loon <vanloon@uva.nl>
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

\keyword{ package }

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
