\name{tfbs.importMotifs}
\alias{tfbs.importMotifs}
\title{ Import licensed motifs to tfbs object }
\description{ Import licensed motifs to tfbs object }
\usage{
tfbs.importMotifs(tfbs, motif_ids, file.pwms)
}
\arguments{
  \item{tfbs}{A tfbs object (\code{"\linkS4class{tfbs}"}) returned by \code{\link{tfbs.createFromCisBP}}, \code{\link{tfbs}}, \code{\link{tfbs.dirs}}. }
  \item{motif_ids}{Vector of motif IDs, motif IDs are in accordance with the TF information which can be exported from TF_information,txt by \code{CisBP.getTFinformation}. }
  \item{file.pwms}{Vector of file names corresponding to the motif IDs specified in 'motif_ids'.}
}

\details{
The motif IDs will be checked according to the TF information in the Cis-BP database.
}

\value{
A new tfbs object (\code{"\linkS4class{tfbs}"}) merged with licensed motifs.
}

\seealso{
\code{\link{tfbs.createFromCisBP}}
}
\examples{

library(rtfbsdb);

db  <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, family_name="AP-2");
tfs;

motif_ids <- c( "M2938_1.02", "M2940_1.02", "M2940_2.02", "M4056_1.02" );

path <- system.file("extdata", package="rtfbsdb");
file_pwms <- paste(path, c("fake_M2938_1.02.pwm", "fake_M2940_1.02.pwm", "M2940_2.02.pwm", "M4056_1.02.pwm"), sep="/");

tfs <- tfbs.importMotifs(tfs, motif_ids, file_pwms );
show(tfs);

}
\keyword{ tfbs object }
