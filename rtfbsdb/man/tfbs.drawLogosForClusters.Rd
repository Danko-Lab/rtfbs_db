\name{tfbs.drawLogosForClusters}
\alias{tfbs.drawLogosForClusters}
\title{
Draw the motif logos by clustering.
}
\description{
Draw the motif logos by one group per page.
}
\usage{
tfbs.drawLogosForClusters(tfbs, cluster.mat, pdf.logos)
}
\arguments{
  \item{tfbs}{A tfbs object.}
  \item{cluster.mat}{ A matrix with 2 columns returned by \code{\link{tfbs.clusterMotifs}}, 1st column is the index of motifs and 2nd column is the group number of clustering.}
  \item{pdf.logos}{ String indicating a PDF eilname.}
}
\details{
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
tfs <- CisBP.find(db, motif_type="ChIP-seq", motif_info_type=1 );

# Calculate the distance matrix 
tfs <- tfbs.getDistanceMatrix( tfs, ncores=1);

# Cluster the motifs using the "cors" method
cluster1 <- tfbs.clusterMotifs(tfs, pdf.heatmap = "test-heatmap1.pdf", method="cors" );
show(cluster1);

# draw motif logos on one group per page.
tfbs.drawLogosForClusters(tfs, cluster1, "test-cluster1.pdf")


}
\keyword{tfbs object}
\keyword{Clustering}
\keyword{Logo}
