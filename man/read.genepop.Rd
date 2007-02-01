\name{read.genepop}
\alias{read.genepop}
\title{ Imports a genepop file }
\description{
  Imports a genepop formatted file.
}

\synopsis{
read.genepop( filename = NULL  )
}

\usage{
read.genepop( filename = "butterfly_borneo" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ path (absolute or relative) & filename of genepop file to be read. }
}

\value{
The function returns an object that holds information regarding the genepop file.
The output can be used by the function \code{aresCalc}.
}

\references{

Raymond M, Rousset F (1995) GENEPOP (version 1.2): population genetics software for exact tests and ecumenicism. J. Heredity, 86, 248-249.

\url{http://wbiomed.curtin.edu.au/genepop/help_input.html}

}

\author{ Scott Davis, based on code by Emiel van Loon }

\examples{

# import of genepop file butterfly_borneo
path <- system.file( package = 'ARES' );
file <- paste( list( path, '/data/butterfly_borneo' ), collapse="" );

butterfly_data <- read.genepop( filename=file );

# select fourth population from the imported genepop file
I298 <- butterfly_data[[4]];

# the data is accessed through the attribute 'output'
output <- aresCalc(I298@output, bootsize=4, maxsize=100);
aresPlot(output, T ='I2, 1998');

# note!!! this example uses a very small bootstrap size (for speed of display).
#  A more realistic size is bootsize=200.


}
\keyword{ models }

