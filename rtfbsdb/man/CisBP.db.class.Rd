\name{CisBP.db-class}
\Rdversion{1.1}
\docType{class}
\alias{CisBP.db-class}
\alias{tfbs.createFromCisBP,CisBP.db-method}
\alias{CisBP.group,CisBP.db-method}

\title{Class \code{"CisBP.db"}}
\description{
The motif library from CisBP web site.\cr
Link: \code{\link{http://cisbp.ccbr.utoronto.ca/}}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{\link{CisBP.extdata}}, \code{\link{CisBP.zipload}}, \code{\link{CisBP.download}}.
}
\section{Slots}{
  \describe{
    \item{\code{species}:}{String indicating the species name defined in the CisBP dataset.}
    \item{\code{zip.file}:}{String indicating the filename of temporary data file.}
    \item{\code{zip.url}:}{String indicating the download source.}
    \item{\code{file.tfinfo}:}{String indicating the TF filename, defulat is TF_Information.txt.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{tfbs.db}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\link{tfbs.createFromCisBP}}{Build a tfbs object by querying the meta file of CisBP dataset and subsetting the results.}
    \item{\link{CisBP.group}}{Get the statistical summary by grouping the field in the CisBP dataset.}
	 }
}
\references{
Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero, A., Drewe, P., ... & Hughes, T. R. (2014). Determination and inference of eukaryotic transcription factor sequence specificity. Cell, 158(6), 1431-1443.
}

\seealso{
}
\examples{
showClass("CisBP.db")
}
\keyword{classes}
