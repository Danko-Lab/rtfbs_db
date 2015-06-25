\name{tfbs.drawLogosForClusters}
\alias{tfbs.drawLogosForClusters}
\title{
Draw the motif logos by clustering.
}
\description{
Draw the motif logos by one cluster per page.
}
\usage{
tfbs.drawLogosForClusters(tfbs, cluster.mat, file.pdf )
}
\arguments{
  \item{tfbs}{A tfbs object(\code{"\linkS4class{tfbs}"}).}
  \item{cluster.mat}{ A matrix with 2 columns returned by \code{\link{tfbs.clusterMotifs}}, 1st column is the index of motifs and 2nd column is the group number of clustering.}
  \item{file.pdf }{ String indicating a PDF eilname.}
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
