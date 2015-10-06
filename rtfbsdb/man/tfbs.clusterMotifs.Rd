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
      method = c("agnes", "apcluster"), 
      pdf.heatmap = NA, 
      group.k = NA,
      apcluster.q = 0.99,
      ncores = 1,
      BG = log(c(0.25, 0.25, 0.25, 0.25)) )
}
\arguments{
  \item{tfbs}{A tfbs object (\code{"\linkS4class{tfbs}"})  returned by \code{\link{tfbs.createFromCisBP}}, \code{\link{tfbs.dirs}} or other functions.}
  \item{pdf.heatmap}{String, a PDF filename for heatmap.}
  \item{method}{String, availabe values are "agnes" and "apcluster".}
  \item{group.k}{Integer, if the method of agnes is used to do clustering, the parameter of k is optional to use as preset group number.}
  \item{apcluster.q}{Numeric value between 0 and 1, if the method of apcluster is used to do clustering, the parameter of q is optional to use as preset group number.}
  \item{ncores}{Number, the number of cores to use simultaneously.}
  \item{BG}{The log value of probabilities for nucleotide A, C, G and T as Backgroud computing.}
}
\details{
This result of clustering will be used in the \code{\link{tfbs.drawLogosForClusters}}, \code{\link{tfbs.selectByGeneExp}}, \code{\link{tfbs.enrichmentTest}}.

tfbs@cluster will be updated by the clustering matrix which 1st column is the index of motifs and 2nd column is the group number of clustering.

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

# Cluster the motifs using the "agnes" method
tfs <- tfbs.clusterMotifs(tfs, method="agnes", pdf.heatmap = "test-heatmap-agnes.pdf" );
show(tfs@cluster);

# Cluster the motifs using the "apcluster" method
tfs <- tfbs.clusterMotifs(tfs, method="apcluster", pdf.heatmap = "test-heatmap-apcluster.pdf" );
show(tfs@cluster);

# draw motif logos on one group per page.
tfbs.drawLogosForClusters(tfs, file.pdf = "test-cluster-logos.pdf");

}

\keyword{Clustering}
\keyword{tfbs object}
