\name{tfbs.drawLogosForClusters}
\alias{tfbs.drawLogosForClusters}
\title{
Draw the motif logos by clustering.
}
\description{
Draw the motif logos by one cluster per page.
}
\usage{
tfbs.drawLogosForClusters(tfbs, file.pdf, nrow.per.page=6 )
}
\arguments{
  \item{tfbs}{A tfbs object(\code{"\linkS4class{tfbs}"}).}
  \item{file.pdf }{ String indicating a PDF eilname.}
  \item{nrow.per.page }{ Number indicating row count in each page.}
}
\details{
It is different with \code{\link{tfbs.drawLogo}} which is capable of printing out motif logos in group. This group is calculated by the \code{\link{tfbs.clusterMotifs}}, not is classfied by any group filed.
}
\value{
No return value.
}
\seealso{
See Also as \code{\link{tfbs.clusterMotifs}}
}
\examples{

# Load the internal CisBP data set
db <- CisBP.extdata("Homo_sapiens");

# Create a tfbs object by querying the meta file of CisBP dataset. 
tfs <- tfbs.createFromCisBP(db, motif_type="ChIP-seq", tf.information.type=1 );

# Cluster the motifs using the "cors" method
tfs <- tfbs.clusterMotifs(tfs, method="apcluster", 
      pdf.heatmap = "test-heatmap1.pdf" );
show(tfs@cluster);

# draw motif logos on one group per page.
tfbs.drawLogosForClusters(tfs, "test-cluster1.pdf")

# draw motif logos on one group per page.
tfbs.drawLogosForClusters(tfs, "test-cluster2.pdf", nrow.per.page=5 )


}
\keyword{tfbs object}
\keyword{Clustering}
\keyword{Logo}
