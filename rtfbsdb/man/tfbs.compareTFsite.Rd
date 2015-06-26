\name{tfbs.compareTFsite}
\alias{tfbs.compareTFsite}
\title{
Comparative TS sites between positive and negative TRE loci
}
\description{
Comparative TS sites between positive and negative TRE loci for all motifs.
}
\usage{
tfbs.compareTFsite(tfbs, 
      file.twoBit, 
      positive.bed, 
      negative.bed, 
      file.prefix=NA, 
      usemotifs=NA,
      ncores=3,
      negative.correction=FALSE, 
      fdr=0.1, 
      threshold=NA, 
      gc.groups=4, 
      background.order=2, 
      background.length=100000, 
      pv.adj = p.adjust.methods)
}
\arguments{
  \item{tfbs}{A tfbs object, see also \code{"\linkS4class{tfbs}"} }
  \item{file.twoBit}{String, the file name of genome data( e.g. hg19.2bit, mm10.2bit) }
  \item{positive.bed}{Data frame, bed-formatted TRE loci.}
  \item{negative.bed}{Data frame, bed-formatted background loci.}
  \item{file.prefix}{String, the prefix for outputted BED file, no bed files output if NA }
  \item{usemotifs}{Vector of index, the index of which motifs are used to compare. All motifs are used if NA. }
  \item{ncores}{Number, comupting nodes in parallel environment.(default=3)}
  \item{negative.correction}{Logical value, if the difference between positive and negative TREs is significant,the resampling will be applied to the correction for the negative TREs. (default=FALSE) }
  \item{fdr}{Numeric value between 0 and 1, False Discovery Rate (FDR) of possible binding sites in \code{rtfbs} package, only binding sites with FDR less than this value can be selected.(default=0.1)}
  \item{threshold}{Numeric value, only sites with scores above this threshold are returned, not be used if NA. (default = NA)}
  \item{gc.groups}{Numeric value,indicating number of quantiles to group sequences into in \code{rtfbs} package. (default = 4)}
  \item{background.order}{Number, order of Markov model to build background.(default=2).}
  \item{background.length}{Number, length of the sequence to simulate background.(default=100000).}
  \item{pv.adj}{String, P-values correct method for \code{p.adjust} function. The available values are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr" or "none".(default="bonferroni").}
}
\details{
The difference of GC contents between positive.bed and negative.bed is checked before the comparson. The p-value of Wilcox test and a vioplot figure 
show this difference and help the user to determine whether the correction is necessary. If the difference is very significant, please set 
\emph{background.correction} to do background correction by resampling the TREs from negative bed data based on the frequency of TREs in 
negative bed data.

}
\value{
A object with the class name of "tfbs.comparson" will be resturned in this comparson function. It includes one list of parameters \code{parm} and one data frame of results \code{result}.\cr


\code{result} is a data frame with the following columns: 
  \item{motif.id}{Motif ID.}
  \item{tf.name}{TF name.}
  \item{Npos}{TF site count found in positive ranges.}
  \item{Nneg}{TF site count found in negative ranges.}
  \item{es.ratio}{Enrichment score.}
  \item{pvalue}{p-value calculated by fisher test.}
  \item{pv.adj}{p-value corrected by the multiple correction.}
  \item{starch}{Binary filename of detected TF sites.}
  
The \code{result} can be outputted to a report by the function \code{\link{tfbs.reportComparson}}. 

}
\seealso{
\code{\link{print.tfbs.comparson}}, \code{\link{summary.tfbs.comparson}}, \code{\link{tfbs.reportComparson}}.
}

\examples{

library(rtfbsdb);

file.twoBit     <- "/local/storage/data/hg19/hg19.2bit"

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, family_name="AP-2");


#make two dummy BED data frame for positive loci and negative loci
pos.bed <- data.frame(chr="chr1", 
		start=round(runif(100,1000000, 2000000)), 
		stop=0, 
		name="", 
		score=0, 
		strand=".");
pos.bed$stop <- pos.bed$start + 3000;

neg.bed <- data.frame(chr="chr1", 
		start=round(runif(200, 800000, 1800000)), 
		stop=0, 
		name="", 
		score=0, 
		strand=".");
neg.bed$stop <- neg.bed$start + round(runif(200, 1000, 3000));

t1 <- tfbs.compareTFsite( tfs, 
      file.twoBit, 
      pos.bed, 
      neg.bed, 
      negative.correction=TRUE,
      ncores = 1); #ncores=3

#Show a brief result
t1;

#Show the comparson results of all motifs 
show(t1$result);

#Output the result to one pdf report. 
tfbs.reportComparson(tfs, t1, file.pdf="test-tfbs-comp.pdf", sig.only=FALSE);


}
\keyword{tfbs object}
\keyword{Comparing}

