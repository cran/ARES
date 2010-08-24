\name{import.data}
\alias{import.data}
\title{Import Data}
\description{
Import datasets in different formats
}
\usage{
import.data(file, text.format = "csv", ...)
}
\arguments{
  \item{file}{quoted string with the file name}
  \item{text.format}{rules for text files}
  \item{\dots}{further options for \code{\link{read.table}}}
}
\details{
The function will try to select the proper filter based on the file extension. Valid file types are some of those supported by package \pkg{foreign}. There are rules implemented for Stata up to version 9 (*.dta), SPSS up to version 11 (confirmed) (*.sav), xBase files up to version 4 (*.dbf), Epi-Info files version 6 (*.rec). R binary files \code{*.rda} can also be read directly by this function. Excel files (*.xls) are no longer supported. Please, save the files in the \code{CSV} format.

Several text formats can be read by setting \code{text.format} or passing options through \code{\dots}. If none of the above proprietary formats extension is part of the file name, the function will assume an ASCII text filter disregarding the extension. Available options for \code{text.format} are \code{"csv"} for comma separated values, \code{"tab"} for columns separated by a tab character and \code{"spc"} for space separated columns and, also, its variants for comma as the decimal separator, \code{"csv2"}, \code{"tab2"} and \code{"spc2"}.

If \code{file} is set to \code{"clipboard"} the function will try to import data previously copied into the system clipboard using text data rules set in \code{text.format}. It can render unexpected results though.
}
\value{
An object of class \code{data.frame}.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

Ares-Rio Program at \url{http://www.ims.uerj.br/ares-rio}
}
\seealso{\code{\link{attach}},\code{\link{search}}}
\examples{
data(admrio)
export.data(admrio,"admrio.dta")
new.data=import.data("admrio.dta")
}
\keyword{data}
\keyword{datagen}
