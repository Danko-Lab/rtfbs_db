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
      use.cluster=NA,
      ncores=3,
      gc.correction=TRUE, 
      gc.correction.pdf=NA, 
      fdr=NA, 
      threshold=6, 
      gc.groups=1, 
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
  \item{use.cluster}{ Clustering matrix with 2 columns, 1st column is the index of motifs and 2nd column is the group number of clustering. It can be obtained from \code{\link{tfbs.clusterMotifs}}. If no clustering matrix, all motifs are used to do the comparson. see \emph{details} }
  \item{ncores}{Number, comupting nodes in parallel environment.(default=3)}
  \item{gc.correction}{Logical value, if the difference between positive and negative TREs is significant,the resampling will be applied to the correction for the negative TREs. (default=TRUE) }
  \item{gc.correction.pdf}{String, indicating th pdf file name if the GC correction is checked. (default=NA) }
  \item{fdr}{Numeric value between 0 and 1, False Discovery Rate (FDR) of possible binding sites in \code{rtfbs} package, only binding sites with FDR less than this value can be selected.(default=NA)}
  \item{threshold}{Numeric value, only sites with scores above this threshold are returned, not be used if NA. (default = 6)}
  \item{gc.groups}{Numeric value,indicating number of quantiles to group sequences into in \code{rtfbs} package. (default = 1)}
  \item{background.order}{Number, order of Markov model to build background.(default=2).}
  \item{background.length}{Number, length of the sequence to simulate background.(default=100000).}
  \item{pv.adj}{String, P-values correct method for \code{p.adjust} function. The available values are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr" or "none". (default="bonferroni").}
}
\details{
The difference of GC contents between positive.bed and negative.bed is checked before the comparson. The p-value of Wilcoxon-Mann-Whitney test 
shows this difference and helps the user to determine whether the GC correction is necessary.If the difference is very significant, please set 
\code{gc.correction} to do GC content correction by resampling the TREs from negative bed data based on the frequency of TREs in negative bed data. 
Use the parameter of \code{gc.correction.pdf} to output vioplot figurs in a pdf file if you want to check the visualized difference.\cr

The clustering matix indicates which motifs in the 1st column are slected to do comparson and which clustering group in the 2nd columns are applied to adjust p-values for multiple comparsons. 
The function applys the p-values adjust for each clustering group. If no clustering information, all motifs in the tfbs object will be selected and adjusted as one group, which is the most conservative method.

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
      gc.correction=TRUE,
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

