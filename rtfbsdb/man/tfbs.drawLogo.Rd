\name{tfbs.drawLogo}
\alias{tfbs.drawLogo}
\title{
Draw single motif logo
}
\description{
Draw the logo for a single TF motif.
}
\usage{
tfbs.drawLogo(tfbs, index)
}
\arguments{
  \item{tfbs}{A tfbs object(\code{"\linkS4class{tfbs}"})}
  \item{index}{Vector of number, indicating the motif index.}
}
\details{
}
\value{
No return values.
}

\seealso{
See Also as \code{"\linkS4class{tfbs}"}
}
\examples{
 
 db <- CisBP.extdata("Homo_sapiens");
 tfs <- tfbs.createFromCisBP(db, family_name="AP-2");
 pdf("test-logos.pdf");
 tfbs.drawLogo(tfs, 1:10 );
 dev.off();
 
 #unlink("test-logos.pdf");
 
}
\keyword{Logo}
\keyword{tfbs object}
