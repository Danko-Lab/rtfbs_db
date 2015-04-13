\name{tfbs.selectByRandom}
\alias{tfbs.selectByRandom}
\title{
Random motif selection
}
\description{
Select the motifs randomly from each group of clustering.
}
\usage{
tfbs.selectByRandom(tfbs, cluster.mat)
}
\arguments{
  \item{tfbs}{A tfbs object.}
  \item{cluster.mat}{ A matrix with 2 columns returned by \code{\link{tfbs.clusterMotifs}}, 1st column is the index of motifs and 2nd column is the group number of clustering.}
}
\details{
The indexes of selected motifs can be used in the function of \code{\link{tfbs.compareTFsite}} or \code{\link{tfbs.scanTFsite}}. 
}
\value{
A vector of motif indexes is returned.
}
\seealso{
See Also as \code{\link{tfbs.selectByGeneExp}}, \code{\link{tfbs.getExpression}}
}
\examples{

db <- CisBP.extdata("Homo_sapiens");

tfs <- CisBP.find(db, family_name="AP-2");

tfs <- tfbs.getDistanceMatrix(tfs, ncores=1);
 
cluster1 <- tfbs.clusterMotifs(tfs, pdf.heatmap="test-AP2-heatmap.pdf" );

usemotif <- tfbs.selectByRandom(tfs, cluster1);

show(usemotif);

}
\keyword{tfbs object}
\keyword{Selection}
\keyword{Clustering}

