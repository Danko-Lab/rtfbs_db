\name{tfbs.compareTFsite}
\alias{tfbs.compareTFsite}
\title{
Comparative TFBS search with the BED ranges
}
\description{

}
\usage{
tfbs.compareTFsite(tfbs, 
      file.twoBit, 
      positive.bed, 
      negative.bed, 
      file.prefix, 
      usemotifs=NA,
      background.correction=FALSE, 
      fdr=0.1, 
      threshold=NA, 
      background.order=2, 
      background.length=100000, 
      ncores=3)
}
\arguments{
  \item{tfbs}{A tfbs object, see also \code{"\linkS4class{tfbs}"} }
  \item{file.twoBit}{String, the file name of genome data( hg19.2bit, mm10.2bit) }
  \item{positive.bed}{Data frame, bed-formatted down-regulatory ranges}
  \item{negative.bed}{Data frame, bed-formatted up-regulatory ranges}
  \item{file.prefix}{String, the prefix for outputted file }
  \item{usemotifs}{Number, the index of which motifs are used to compare. All motifs are used if NA. }
  \item{background.correction}{Logical value, if the difference between positive and negative TREs is significant,the resampling will be applied to the correction for the negative TREs. Default:FALSE. }
  \item{fdr}{Numeric value, False Discovery Rate (FDR) of possible binding sites, only binding sites with FDR less than this value can be selected.}
  \item{threshold}{Numeric value, only sites with scores above this threshold are returned (default = 6)}
  \item{background.order}{ Number, order of Markov model to build background}
  \item{background.length}{ Number, length of the sequence to simulate background}
  \item{ncores}{Number, comupting nodes in parallel environment.}
}
\details{
The background difference between positive.bed and negative.bed is checked before the comparson. The p-value of Wilcox test and a vioplot figure 
show this difference and help the user to determine whether the correction is necessary. If the difference is very significant, please set 
\emph{background.correction} to do background correction by resampling the TREs from negative bed data based on the frequency of TREs in 
negative bed data.

}
\value{
A data frame with the following columns: 
  \item{tf.name}{TF name.}
  \item{Npos}{TF site count found in positive ranges.}
  \item{Nneg}{TF site count found in negative ranges.}
  \item{es.ratio}{Enrichment score.}
  \item{assoc.pvalue}{p-value calculated by fisher test.}
  \item{pv.bonferoni}{p-value corrected by Bonferoni.}
  \item{starch}{Binary filename of detected TF sites.}
  \item{motif.id}{Motif ID.}
}
\seealso{

}
\examples{

file.dREG.H.change.bed <- "/home/zw355/src/rtfbs_db/rtfbsdb/test/dREG.H.change.bed"
file.dREG.all.bed    <- "/home/zw355/src/rtfbs_db/rtfbsdb/test/dREG.all.bed" 
file.twoBit     <- "/local/storage/data/hg19/hg19.2bit"

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, family_name="AP-2");

dREG_H_change_bed <- read.table(file.dREG.H.change.bed, header=FALSE);
dREG_all_bed <- read.table(file.dREG.all.bed, header=FALSE);

t <- tfbs.compareTFsite( tfs, 
      file.twoBit, 
      dREG_H_change_bed, 
      dREG_all_bed, 
      background.correction=TRUE,
      file.prefix="comp.db",
      ncores = 1); #ncores=3

}
\keyword{tfbs object}
\keyword{Comparing}

