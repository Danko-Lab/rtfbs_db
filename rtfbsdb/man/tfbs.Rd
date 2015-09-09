\name{tfbs}
\alias{tfbs}
\title{
Create a tfbs object from the supplied PWM files.
}
\description{
Create a tfbs object from the supplied PWM files.
}
\usage{
tfbs(filenames, 
      names, 
      species="Homo_sapiens", 
      tf_info = NULL,
      tf_missing = NULL, ...)
}
\arguments{
  \item{filenames}{ Vector of PWM files }
  \item{names}{Vector of unique gene symbols. }
  \item{species}{String indicating species name }
  \item{tf_info}{Data frame including meta information copied from CisBP data file for all existing motifs., Default: NULL}
  \item{tf_missing}{Data frame including meta information copied from CisBP data file for missing motifs., Default: NULL}
  \item{\dots}{Parameters,such as pseudocount, force_even, and the parameters used in read.table function. }
}
\details{
Load the PWM files to build a \code{"\linkS4class{tfbs}"} object.
}
\value{
A tfbs object (class: \code{"\linkS4class{tfbs}"}) including all PWM matrics.The all attributes are as follows:

  \item{TFID}{Vector of non-unique ID for TF.}
  \item{species}{String indicating the species name}
  \item{ntfs}{Number of motifs in matrix.}
  \item{pwm}{A list including PWM matics.}
  \item{filename}{Vector of PWM filename.}
  \item{mgisymbols}{Unique gene symbols for TF.}
  \item{tf_info}{Data frame, including extra information for all existing PWMs, it maybe different with motif dataset, default:NULL.}
  \item{tf_missing}{Data frame, including extra information for missing PWMs, it maybe different with motif dataset, default:NULL.}
  \item{distancematrix}{Distance matrix between motifs returned by \code{\link{tfbs.clusterMotifs}}, default:NULL.}
  \item{expressionlevel}{Data frame indicatig the result of expression level returned by \code{\link{tfbs.getExpression}}, default:NULL. }
  \item{cluster}{Matrix with 2 columns returned by \code{\link{tfbs.clusterMotifs}}, 1st column is the index of motifs and 2nd column is the group number of clustering, default:NULL. }

The tfbs object can be created by the function of \code{\link{tfbs}}, \code{\link{tfbs.dirs}}, \code{\link{tfbs.createFromCisBP}}.
}
\seealso{
\code{\link{tfbs}}, \code{\link{tfbs.dirs}}, \code{\link{tfbs.createFromCisBP}}
}
\examples{

# M3590_1.01 PAX5 ENSG00000196092	
# M3590_1.01 PAX5 ENSG00000196092	
fs1 <- system.file("extdata","M3590_1.01.pwm", package="rtfbsdb")
fs2 <- system.file("extdata","M3591_1.01.pwm", package="rtfbsdb")

cat(fs1, "\n");

tfs <- tfbs( c( fs1, fs2 ), names=c("M3590_1.01","M3591_1.01"), 
      header=TRUE, sep="\t" , row.names=1 );
str(tfs);

}
\keyword{tfbs object}

