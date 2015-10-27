\name{tfbs.selectByGeneExp}
\alias{tfbs.selectByGeneExp}
\title{
Motif selection by gene expression level.
}
\description{
Select the motifs with minimum p-value from each group of clustering.
}
\usage{
tfbs.selectByGeneExp(tfbs)
}
\arguments{
  \item{tfbs}{A tfbs object (\code{"\linkS4class{tfbs}"}) with the data frame of gene expression level.}
}
\details{
The function of \code{\link{tfbs.getExpression}} should be successfully called and the results of gene expression should be returned before this function is called. The indexes of selected motifs will be used in the function of \code{\link{tfbs.enrichmentTest}} or \code{\link{tfbs.scanTFsite}}. 

}
\value{
New tfbs object with the selected indices is returned(tfbs@cluster). The 3rd column of tfbs@cluster is added or updated as the select flag(1:selected, 0:unselected).
}
\seealso{
See Also as \code{\link{tfbs.selectByRandom}}, \code{\link{tfbs.getExpression}}
}
\examples{

db <- CisBP.extdata("Homo_sapiens");

tfs <- tfbs.createFromCisBP(db, family_name="AP-2");


file.bigwig.minus <- system.file("extdata",
      "GSM1480327_K562_PROseq_chr19_minus.bw", package="rtfbsdb")
file.bigwig.plus <- system.file("extdata",
      "GSM1480327_K562_PROseq_chr19_plus.bw", package="rtfbsdb")
hg19.twobit <- system.file("extdata","hg19.chr19.2bit", package="rtfbsdb")
gencode.gtf <- system.file("extdata",
      "gencode.v21.annotation.chr19.gtf.gz", package="rtfbsdb")

tfs <- tfbs.getExpression(tfs, hg19.twobit, gencode.gtf, 
      file.bigwig.plus, file.bigwig.minus, ncore=1); 

tfs <- tfbs.clusterMotifs(tfs, pdf.heatmap="test-AP2-heatmap.pdf" );

tfs <- tfbs.selectByGeneExp( tfs );

show(tfs@cluster);

}
\keyword{tfbs object}
\keyword{Clustering}
\keyword{Selection}
