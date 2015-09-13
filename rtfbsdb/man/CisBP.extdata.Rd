\name{CisBP.extdata}
\alias{CisBP.extdata}
\title{
Load internal CisBP dataset.
}
\description{
Build a CisBP object from the internal zip file stored in this package
}
\usage{
CisBP.extdata(species)
}
\arguments{
  \item{species}{ String, only valid for human and mouse species, i.e. Homo_sapiens, Mus_musculus, or Drosophila_melanogaster }
}
\details{
The CisBP data for Homo_sapiens and Mus_musculus are delivered by this package. When you use the newest dataset, you should download it from the website by \code{\link{CisBP.download}}.
}
\value{
A CisBP object (class name: \code{"\linkS4class{CisBP.db}"}) is returned with four items:
  \item{species}{String indicating the species name defined in the CisBP dataset.}
  \item{zip.file}{String indicating the filename of temporary data file.}
  \item{zip.url}{String indicating the download source }
  \item{file.tfinfo}{String indicating the TF filename, default is TF_Information.txt.}
}

\seealso{
See Also as \code{\link{CisBP.zipload}}, \code{\link{CisBP.download}}.
}
\examples{

#reading human data from extension data file in the package
db.human <- CisBP.extdata("Homo_sapiens")

#reading Drosophila_melanogaster from extension data file in the package
db.dm3 <- CisBP.extdata("dm3")

}
\keyword{CisBP object}
