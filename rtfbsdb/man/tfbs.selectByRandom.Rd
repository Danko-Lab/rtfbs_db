\name{tfbs.selectByRandom}
\alias{tfbs.selectByRandom}
\title{
Random motif selection
}
\description{
Select the motifs randomly from each group of clustering.
}
\usage{
tfbs.selectByRandom(tfbs)
}
\arguments{
  \item{tfbs}{A tfbs object(\code{"\linkS4class{tfbs}"}).}
}
\details{
The indexes of selected motifs can be used in the function of \code{\link{tfbs.enrichmentTest}} or \code{\link{tfbs.scanTFsite}}. 
}
\value{
New tfbs object with the selected indices is returned(tfbs@cluster). The 3rd column of tfbs@cluster is added or updated as the select flag(1:selected, 0:unselected).
}
\seealso{
See Also as \code{\link{tfbs.selectByGeneExp}}, \code{\link{tfbs.getExpression}}
}
\examples{

db <- CisBP.extdata("Homo_sapiens");

tfs <- tfbs.createFromCisBP(db, family_name="AP-2");

tfs <- tfbs.clusterMotifs(tfs, pdf.heatmap="test-AP2-heatmap.pdf" );

tfs <- tfbs.selectByRandom(tfs );

show(tfs@cluster);

}
\keyword{tfbs object}
\keyword{Selection}
\keyword{Clustering}

