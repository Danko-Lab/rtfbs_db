\name{tfbs.enrichmentTest}
\alias{tfbs.enrichmentTest}
\title{
Comparative TS sites between positive and negative TRE loci
}
\description{
Comparative TS sites between positive and negative TRE loci for all motifs.
}
\usage{
tfbs.enrichmentTest(tfbs, 
      file.genome, 
      positive.bed, 
      negative.bed=NA, 
      file.prefix=NA, 
      use.cluster=FALSE,
      ncores=1,
      gc.correction=TRUE, 
      gc.correction.pdf=NA, 
      gc.min.sample = 500, 
      gc.robust.rep=NA, 
      threshold = 6, 
      threshold.type = c("score", "fdr"),
      gc.groups=1, 
      background.order=2, 
      background.length=100000, 
      pv.adj = p.adjust.methods)
}
\arguments{
  \item{tfbs}{A tfbs object, see also \code{"\linkS4class{tfbs}"} }
  \item{file.genome}{String, the file name of genome data, 2bit or FastA format( e.g. hg19,fasta, hg19.2bit, mm10.2bit) }
  \item{positive.bed}{Data frame, bed-formatted TRE loci.}
  \item{negative.bed}{Data frame, bed-formatted background loci. If not specified, the genomic loci adjacent to positive one are randomly extracted as the negative bed. }
  \item{file.prefix}{String, the prefix for outputted BED file, no bed files output if NA }
  \item{use.cluster}{ Clustering matrix with 2 columns, 1st column is the index of motifs and 2nd column is the group number of clustering. It can be obtained from \code{\link{tfbs.clusterMotifs}}. If no clustering matrix, all motifs are used to do the comparson. see \emph{details} }
  \item{ncores}{Number, comupting nodes in parallel environment.(default=1)}
  \item{gc.correction}{Logical value, if the difference between positive and negative TREs is significant,the resampling will be applied to the correction for the negative TREs. (default=TRUE) }
  \item{gc.correction.pdf}{String, indicating the pdf file name if the GC correction is checked. (default=NA) }
  \item{gc.min.sample}{Number, indicating minimum sample size when resamping the negative genomic loci for the GC correct. (default=500) }
  \item{gc.robust.rep}{Number, indicating whether resampling background set multiple times is applied to get the median of binding sites. (default=NA) }
  \item{threshold}{Numeric value, if 'score' is specified in \code{threshold.type}, only binding sites with scores above this threshold are returned, 
  if 'fdr' is specified in \code{threshold.type}, only binding sites with FDR (False Discovery Rate) less than this value can be selected. Default value is 6 for 'score' and 0.1 for 'fdr'.} 
  \item{threshold.type}{String value, two options are available. only sites with scores above this threshold are returned, not be used if NA. (default = 'score')}
  \item{gc.groups}{Numeric value,indicating number of quantiles to group sequences into in \code{rtfbs} package. (default = 1)}
  \item{background.order}{Number, order of Markov model to build background.(default=2).}
  \item{background.length}{Number, length of the sequence to simulate background.(default=100000).}
  \item{pv.adj}{String, P-values correct method for \code{p.adjust} function. The available values are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr" or "none". (default="bonferroni").}
}
\details{
(1)The difference of GC contents between positive.bed and negative.bed is checked before the comparson. The p-value of Wilcoxon-Mann-Whitney test 
shows this difference and helps the user to determine whether the GC correction is necessary.If the difference is very significant, please set 
\code{gc.correction} to do GC content correction by resampling the TREs from negative bed data based on the frequency of TREs in negative bed data. 
Use the parameter of \code{gc.correction.pdf} to output vioplot figurs in a pdf file if you want to check the visualized difference.\cr

(2)The clustering matix indicates which motifs in the 1st column are slected to do comparson and which clustering group in the 2nd columns are applied to adjust p-values for multiple comparsons. 
The function applys the p-values adjust for each clustering group. If no clustering information, all motifs in the tfbs object will be selected and adjusted as one group, which is the most conservative method.


(3)This function heavily relys on some Unix commands to operate bed data and gencode data. Please make sure
the following commands work normally in R terminal.\cr

\code{starch}, \code{sort-bed}, \code{bedtools},\code{twoBitToFa}.\cr

The function \code{\link{Sys.which}} and \code{\link{system}} help you to locate these commands and test its availability.\cr

}
\value{
A object with the class name of "tfbs.enrichment" will be resturned in this comparson function. It includes one list of parameters \code{parm} and one data frame of results \code{result}.\cr


\code{result} is a data frame with the following columns: 
  \item{motif.id}{Motif ID.}
  \item{tf.name}{TF name.}
  \item{Npos}{TF site count found in positive ranges.}
  \item{expected}{TF site count found in negative ranges.}
  \item{fe.ratio}{Ratio of fold enrichment.}
  \item{pvalue}{p-value calculated by fisher test.}
  \item{pv.adj}{p-value corrected by the multiple correction.}
  \item{starch}{Binary filename of detected TF sites.}
  
The \code{result} can be outputted to a report by the function \code{\link{tfbs.reportEnrichment}}. 

}
\seealso{
\code{\link{print.tfbs.enrichment}}, \code{\link{summary.tfbs.enrichment}}, \code{\link{tfbs.reportEnrichment}}.
}

\examples{

library(rtfbsdb);

file.twoBit <- system.file("extdata","hg19.chr19.2bit", package="rtfbsdb")

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, family_name="AP-2");


#make two dummy BED data frame for positive loci and negative loci
pos.bed <- data.frame(chr="chr19", 
      start=round(runif(1000,1000000, 2000000)), 
      stop=0, 
      name="", 
      score=0, 
      strand=".");
pos.bed$stop <- pos.bed$start + round(runif(1000, 20, 30));

neg.bed <- data.frame(chr="chr19", 
      start=round(runif(8000, 800000, 1800000)), 
      stop=0, 
      name="", 
      score=0, 
      strand=".");
neg.bed$stop <- neg.bed$start + round(runif(8000, 20, 30));

t1 <- tfbs.enrichmentTest( tfs, 
      file.twoBit, 
      pos.bed, 
      neg.bed, 
      gc.correction=TRUE,
      gc.min.sample = 1000,
      ncores = 1); #ncores=3

#Show a brief result
t1;

#Show the comparson results of all motifs 
show(t1$result);

summary(t1);

#Output the result to one pdf report. 
tfbs.reportEnrichment(tfs, t1, file.pdf="test-tfbs-enrich-all.pdf", sig.only=FALSE);

file.ELF1 <- system.file("extdata","Chipseq-k562-chr19-ELF1.bed", package="rtfbsdb")
pos.bed<- read.table(file.ELF1)

tfs <- tfbs.createFromCisBP(db, family_name="Ets");

t2 <- tfbs.enrichmentTest( tfs, 
      file.twoBit, 
      pos.bed, 
      neg.bed, 
      gc.correction=TRUE,
      gc.robust.rep=5,
      ncores = 1); #ncores=3

show(t2)

#Output the result to one pdf report. 
tfbs.reportEnrichment(tfs, t2, file.pdf="test-tfbs-enrich-both.pdf", 
      sig.only=TRUE, enrichment.type="both");

#plot QQ-like plot for the pvalues of all motifs
tfbs.plotEnrichment(tfs, t2, file.pdf="plot-tfbs-enrich-both.pdf", 
      enrichment.type="both", 
      plot.title="QQ plot", 
      top.motif.labels=3, bottom.motif.labels=3, plot.type="polar", color.scheme=1);


t3 <- tfbs.enrichmentTest( tfs, 
      file.twoBit, 
      pos.bed, 
      gc.correction=TRUE,
      gc.robust.rep=5,
      ncores = 1); #ncores=3

show(t3)

#Output the result to one pdf report. 
tfbs.reportEnrichment(tfs, t3, file.pdf="test-elf1-enrich-depleted.pdf", 
      sig.only=TRUE, enrichment.type="depleted");


#plot QQ-like plot for the pvalues of all motifs
tfbs.plotEnrichment(tfs, t3, file.pdf="plot-elf1-enrich-enriched.pdf", 
      enrichment.type="enriched", 
      plot.title="QQ plot", 
      top.motif.labels=3, plot.type="nonpolar", color.scheme=2);


#Use FastA file to call enrichment test
file.fastfa = system.file("extdata","dna.fasta", package="rtfbsdb")

if(0)
{
t4 <- tfbs.enrichmentTest( tfs, 
      file.fastfa, 
      pos.bed, 
      neg.bed, 
      gc.correction=TRUE,
      gc.min.sample = 1000,
      ncores = 1); #ncores=3

t4;
}


}
\keyword{tfbs object}
\keyword{Enrichment}

