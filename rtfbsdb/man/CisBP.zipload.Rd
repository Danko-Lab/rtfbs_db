\name{CisBP.zipload}
\alias{CisBP.zipload}
\title{
Load the zipped CisBP file.
}
\description{
Build a CisBP object from the zipped CisBP file.
}
\usage{
CisBP.zipload(zip.file, species = "Homo_sapiens")
}
\arguments{
  \item{zip.file}{ String, indicating the zipped file data }
  \item{species}{ String, indicating the species name in the CisBP database }
}
\details{
The zip data canbe downloaded from the web site, please check \code{\link{CisBP.download}}.
}

\value{
A CisBP object (class name: \code{"\linkS4class{CisBP.db}"}) is returned with four items:
  \item{species}{String indicating the species name}
  \item{zip.file}{String indicating the filename of temporary data file.}
  \item{zip.url}{String indicating the download source }
  \item{file.tfinfo}{String indicating the TF filename, default is TF_Information.txt.}
}
\seealso{
See Also as \code{\link{CisBP.extdata}}, \code{\link{CisBP.download}}.
}
\examples{

# Download the dataset 
db2 <- CisBP.download("Mus_musculus");

# Loading the zip file, the db2 and db3 have same TF data. 
# Here is an example to show how to use CisBP.zipload. 
# We dont nee to download it by CisBP.download and then load it by CisBP.zipload
db3 <- CisBP.zipload(db2@zip.file, species="Mus_musculus");

}
\keyword{CisBP object}
