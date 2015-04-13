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
  \item{tfbs}{A tfbs object}
  \item{index}{Number, indicating the motif index.}
}
\details{
}
\value{
No return values.
}

\seealso{
See Also as \code{\link{tfbs}}
}
\examples{
 
 db <- CisBP.extdata("Homo_sapiens");
 tfs <- CisBP.find(db, family_name="AP-2");
 pdf("test-logos.pdf");
 for(i in 1:10) tfbs.drawLogo(tfs, i );
 dev.off();
 
 #unlink("test-logos.pdf");
 
}
\keyword{Logo}
\keyword{tfbs object}
