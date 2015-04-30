\name{tfbs.dirs}
\alias{tfbs.dirs}
\title{
Create a tfbs object from the folders.
}
\description{
Create a tfbs object from all the PWM files found in the supplied folders.
}
\usage{
tfbs.dirs(..., 
      species = "Homo_sapiens", 
      args.read.motif = NULL, 
      pattern = glob2rx("*.pwm"), 
      recursive = FALSE)
}
\arguments{
  \item{\dots}{ Multiple strings, one or more folders can be used in this function. }
  \item{species}{ String, including the species name.}
  \item{args.read.motif}{ List, including \emph{pseudocount}, \emph{force_even} or other parameters used in \code{read.table} function.}
  \item{pattern}{String, a character vector specifying regular expression and wlidcards. }
  \item{recursive}{Logical, indicating the loading recursively descends into subfolders or not, default: FALSE.}
}
\details{

Two parameters in the list of \code{args.read.motif} can be used:\cr
pseudocount: log value for zero value in PWM matrix, default is -7.\cr
force_even: whether the PWM matrix with odd size needs to be even.\cr

}
\value{
A tfbs object collecting all the PWM files in the specified folders. For the details of tfbs object, please see \code{\link{tfbs}}
}
\seealso{
The structure of tfbs object is described in \code{"\linkS4class{tfbs}"}
}
\examples{

fs.dir <- system.file("extdata","", package="rtfbsdb")
tfs <- tfbs.dirs( fs.dir, 
      args.read.motif = list(pseudocount=-7, header=TRUE, sep="\t" , row.names=1) );
str(tfs);

}

\keyword{tfbs object}
