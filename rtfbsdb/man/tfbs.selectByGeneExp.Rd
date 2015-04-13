\name{tfbs.selectByGeneExp}
\alias{tfbs.selectByGeneExp}
\title{
Motif selection by gene expression level.
}
\description{
Select the motifs with minimum p-value from each group of clustering.
}
\usage{
tfbs.selectByGeneExp(tfbs, cluster.mat)
}
\arguments{
  \item{tfbs}{A tfbs object with the data frame of gene expression level.}
  \item{cluster.mat}{ A matrix with 2 columns returned by \code{\link{tfbs.clusterMotifs}}, 1st column is the index of motifs and 2nd column is the group number of clustering.}
}
\details{
The function of \code{\link{tfbs.getExpression}} should be successfully called and te results of gene expression should be returned before this function is called. The indexes of selected motifs will be used in the function of \code{\link{tfbs.compareTFsite}} or \code{\link{tfbs.scanTFsite}}. 
}
\value{
A vector of motif indexes is returned.
}
\seealso{
See Also as \code{\link{tfbs.selectByRandom}}, \code{\link{tfbs.getExpression}}
}
\examples{

db <- CisBP.extdata("Homo_sapiens");

tfs <- CisBP.find(db, family_name="AP-2");

if(0)
{
tfs <- tfbs.getExpression(tfs, file.bigwig.plus, file.bigwig.minus, file.hg19);

cluster1 <- tfbs.clusterMotifs(tfs, pdf.heatmap="test-AP2-heatmap.pdf" );

usemotif <- tfbs.selectByGeneExp(tf, cluster1);
}


}
\keyword{tfbs object}
\keyword{Clustering}
\keyword{Selection}
