\name{tfbs.clusterMotifs}
\alias{tfbs.clusterMotifs}
\title{
Clustering the specified motifs and drawing the heatmap.
}
\description{
Clustering the specified motifs and drawing the heatmap.
}
\usage{
tfbs.clusterMotifs(tfbs, 
      subset = NA, 
      pdf.heatmap = NA, 
      method = NA, 
      group.k = NA)
}
\arguments{
  \item{tfbs}{A tfbs object (\code{"\linkS4class{tfbs}"})  returned by \code{\link{tfbs.createFromCisBP}}, \code{\link{tfbs.dirs}} or other functions.}
  \item{subset}{Vector, the indexes of partial motifs if not all motifs are clustered.}
  \item{pdf.heatmap}{String, a PDF filename for heatmap.}
  \item{method}{String, availabe values are "agnes" and "cors".}
  \item{group.k}{Integer, if the method of agnes is used to do clustering, the parameter of k is optional to use as preset group number.}
}
\details{
This result of clustering will be used in the 
}
\value{
A matrix with 2 columns is returned, 1st column is the index of motifs and 2nd column is the group number of clustering.
}

\seealso{
See Also as \code{\link{tfbs.selectByGeneExp}}  and \code{\link{tfbs.selectByRandom}}
}
\examples{

# Load the internal CisBP data set
db <- CisBP.extdata("Homo_sapiens");

# Create a tfbs object by querying the meta file of CisBP dataset. 
tfs <- tfbs.createFromCisBP(db, motif_type="ChIP-seq", tf.information.type=1 );

# Calculate the distance matrix 
tfs <- tfbs.getDistanceMatrix( tfs, ncores=1 );

# Cluster the motifs using the "cors" method
cluster1 <- tfbs.clusterMotifs(tfs, pdf.heatmap = "test-heatmap1.pdf", method="cors" );
show(cluster1);

# draw motif logos on one group per page.
tfbs.drawLogosForClusters(tfs, cluster1, "test-cluster1.pdf")

# Cluster the motifs using the "agnes" function
cluster2 <- tfbs.clusterMotifs(tfs, pdf.heatmap = "test-heatmap2.pdf", method="agnes" );
show(cluster2);

}

\keyword{Clustering}
\keyword{tfbs object}
