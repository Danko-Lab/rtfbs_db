\name{CisBP.download}
\alias{CisBP.download}
\title{ Download CisBP dataset.}
\description{
Download TF data file from CisBP dataset and store it to temporary folder
}
\usage{
CisBP.download(species = "Homo_sapiens", 
      url = "http://cisbp.ccbr.utoronto.ca/bulk_archive.php")
}
\arguments{
  \item{species}{ String, indicating the species name in the CisBP dataset }
  \item{url}{ String, the URL of bulk dowbnloads from CisBP dataset, \cr default is \code{http://cisbp.ccbr.utoronto.ca/bulk_archive.php}}
}
\details{
The dowload function has been confirmed in the web site of \emph{cisbp.ccbr.utoronto.ca} o June, 2015.
}
\value{
A CisBP object (class name: \code{"\linkS4class{CisBP.db}"}) is returned with four items:
  \item{species}{String indicating the species name}
  \item{zip.file}{String indicating the filename of temporary data file.}
  \item{zip.url}{String indicating the download source }
  \item{file.tfinfo}{String indicating the TF filename, default is TF_Information.txt.}
}
\references{
Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero, A., Drewe, P., ... & Hughes, T. R. (2014). Determination and inference of eukaryotic transcription factor sequence specificity. Cell, 158(6), 1431-1443.
}

\seealso{
See Also as \code{\link{CisBP.zipload}}, \code{\link{CisBP.extdata}}.
}
\examples{

#download human dataset
db1 <- CisBP.download("Homo_sapiens");

#download mouse dataset
db2 <- CisBP.download("Mus_musculus");

}
\keyword{CisBP object}

