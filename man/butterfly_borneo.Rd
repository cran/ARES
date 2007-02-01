\encoding{latin1}
\name{butterfly_borneo}
\alias{butterfly_borneo}
\alias{bb100}
\alias{bb198}
\alias{bb200}
\alias{bb204}
\alias{bb297}
\alias{bb298}
\alias{bb299}
\alias{bi100}
\alias{bi200}
\alias{bi298}
\alias{bi300}
\docType{data}
\title{ Microsatellite data from the tropical butterfly Drupadia theda (Felder)}
\description{
This data set contains a number of dataframes: bb100, bb198, bb200, bb204, bb297,
bb298, bb299, bi100, bi200, bi298, bi300.
Each row in the dataframes refers to an allele, and each column to an individual 
butterfly that has been sampled.
Each dataframe has 74 rows and a varying number of columns. 
All individuals belong to the species Drupadia theda (Felder), 
a rainforest lycaenid butterfly typical of pristine and moderately 
disturbed forest. 
The naming of the dataframes is as follows. The second and third charcters 
(b1, b2, i1, i2 and i3) refer to different locations where the data were
colletcted, and the last two characters (97,98,99,00,04) refer to the year
when the data were collected.
The data were originally collected to investigate the 
species-genetic diversity correlation (as introduced by Vellend and Geber, 2005) in 
disturbed landscapes. Butterflies were collected in five differentially 
disturbed landscapes (i1, i2, i3, b1 and b2) in East Kalimantan, Indonesian 
Borneo. Also a time series of samples were taken in one landscape (b2) 
before and after the 1997/98 ENSO (El Niño Southern Oscillation) event. 
The species was sufficiently abundant in all environments to enable reliable 
sample size for robust estimates of allelic richness. However, sample sizes 
differed among landscapes (sample sizes were much lower in disturbed than in 
pristine forests). A more elaborate description of landscapes and data is 
available in Cleary et al. (2004, 2006).
Note that the same data is also available in a genepop-format
\(butterfly_borneo.dat\). The genepop format can be read into R with the function
 \code{read.genepop}.

}

\usage{
data(butterfly_borneo)
}
\format{
  The data frames (bb100, bb198, bb200, bb204, bb297, bb298, bb299, bi100, bi200, bi298, bi300)
  contain logical data (TRUE or FALSE). Each data frame has 74 rows, and a varying number
  of columns. 
}

\references{
Cleary DFR, Mooers AØ, Eichhorn KAO, van Tol J, de Jong R, Menken SBJ (2004) Diversity and community composition of butterflies and odonates in an ENSO-induced fire affected habitat mosaic: a case study from East Kalimantan, Indonesia. Oikos, 105, 426-446.

Cleary DFR, Fauvelot C, Genner MJ, Menken SBJ, Mooers AØ (2006) Parallel responses of species and genetic diversity to El Niño Southern Oscillation-induced environmental destruction. Ecology Letters, 9, 304-310.

Vellend M, Geber MA (2005) Connections between species diversity and genetic diversity. Ecology Letters 8, 767-781.
}

\examples{
data(butterfly_borneo)
}
\keyword{datasets}
