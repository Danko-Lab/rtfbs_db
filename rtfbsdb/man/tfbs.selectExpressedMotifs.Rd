\name{tfbs.selectExpressedMotifs}
\alias{tfbs.selectExpressedMotifs}
\title{
Select expressed Motifs for GRO-seq, PRO-seq and RNA-seq data
}
\description{ 
Select expressed Motifs for GRO-seq, PRO-seq and RNA-seq data
}
\usage{
tfbs.selectExpressedMotifs(tfbs, 
      file.twoBit, 
      file.gencode.gtf, 
      file.bigwig.plus=NA, 
      file.bigwig.minus=NA, 
      file.bam=NA, 
      seq.datatype= c("GRO-seq", "PRO-seq", "RNA-seq"),
      pvalue.threshold = 0.05, 
      lowest.reads.RPKM = NA,
      include.DBID.missing=TRUE, 
      use.strand = FALSE,
      ncores = 1 )       

}
\arguments{
  \item{tfbs}{A tfbs object (\code{"\linkS4class{tfbs}"}) returned by \code{\link{tfbs.createFromCisBP}}, \code{\link{tfbs}}, \code{\link{tfbs.dirs}}. }
  \item{file.bigwig.plus}{String, indicating bigwig file for strand plus(+) if \code{seq.datatype} is GRO-seq or PRO-seq.}
  \item{file.bigwig.minus}{String, indicating bigwig file for strand minus(-) if \code{seq.datatype} is GRO-seq or PRO-seq.}
  \item{file.bam}{String, indicating BAM file for rna reads if \code{seq.datatype} is RNA-seq.}
  \item{file.twoBit}{String, indicating the binary data of sequence. (e.g. hg19.2bit, mm10.2bit)}
  \item{file.gencode.gtf}{String,indicating Gencode GTF file downloaded from the Gencode web site.}
  \item{seq.datatype}{String,indicating which kind of seq data is applied to this function, three values are available: GRO-seq, PRO-seq and RNA-seq. Default: GRO-seq}
  \item{pvalue.threshold}{Numeric, indicating p-value criteria for expressed motifs. Default: 0.05 }
  \item{lowest.reads.RPKM}{Numeric, implying the motifs with lower reads than this threshold will be removed from expressed list. Default: NA}
  \item{include.DBID.missing}{Logical, indicating whether the TFs without association with GENCODE through the DBID are selected. }
  \item{use.strand}{Logical, indicating whether same strandedness is required when getting the number of mapped reads from BAM files for RNA-seq. }
  \item{ncores}{Number, comuputing nodes in parallel environment for gencode data converting.}
}
\details{

1) If \code{seq.datatype} is GRO-seq or PRO-seq and the bigwig files are provided, the gene expression values are calculated 
through querying the TREs region from the GENCODE database( for human, gencode.v19.annotation.gtf, for mouse: gencode.vM3.annotation.gtf)\cr
and querying the reads count in the plus and minus bigWig files.\cr

If \code{seq.datatype} is RNA-seq and the BAM file is provided, read counts for each TRE regions will be queried from the BAM file.\cr

2) If the expressed TFs only is used in the tfbs object, the TFs with p-values corrected by Bonfrroni less than 0.05 will be selected.\cr

The following part explains how to calculate the gene expression.\cr

For each motif, the occurance ranges can be queried by the gene ID  
After the searching, one range obtained from the merge of the multiple ranges will be 
used to detect the reads count in the specified  bigwig files(including plus and minus). 
The probability of each motif can be calcuated by the reads count and lambda.\cr

The lambda is determined by the following formulation:\cr

For GRO-seq and PRO-seq data:\cr
\code{r.lambda = 0.04 * sum(reads_in_all_chromosomes)/10751533/1000.}\cr

For RNA-seq data:\cr

\code{r.lambda = mode( reads_in_1000_bp_windows_cross_all_gene_deserts )/1000.}\cr

3) This function heavily relys on some Unix commands to operate bed data and gencode data. Please make sure
the following commands work normally in R terminal.\cr

\code{awk}, \code{zcat}, \code{sort-bed}, \code{twoBitInfo}, \code{bedtools}.\cr

The function \code{\link{Sys.which}} and \code{\link{system}} help you to locate these commands and test its availability.\cr

4) This function will be failed to get the reads count if the BAM file is not indexed. Please use the command samtools to make the index file for the BAM file\cr

\code{samtools index your_bam_file}\cr


}
\value{
A new tfbs object (\code{"\linkS4class{tfbs}"}) with the matrix of gene expression level.

}
\examples{

library(rtfbsdb);

# Load the internal CisBP data set
db.human <- CisBP.extdata("Homo_sapiens");

# Create a tfbs object by querying the meta file of CisBP dataset. 
tfs <- tfbs.createFromCisBP(db.human, motif_type="ChIP-seq", 
      tf.information.type=1 );

file.bigwig.minus <- system.file("extdata",
      "GSM1480327_K562_PROseq_chr19_minus.bw", package="rtfbsdb")
file.bigwig.plus <- system.file("extdata",
      "GSM1480327_K562_PROseq_chr19_plus.bw", package="rtfbsdb")
hg19.twobit <- system.file("extdata","hg19.chr19.2bit", package="rtfbsdb")
gencode.gtf <- system.file("extdata",
      "gencode.v19.annotation.chr19.gtf.gz", package="rtfbsdb")

tfs1 <- tfbs.selectExpressedMotifs(tfs, 
      hg19.twobit, 
      gencode.gtf, 
      file.bigwig.plus, 
      file.bigwig.minus, 
      seq.datatype = "PRO-seq",
      pvalue.threshold=0.001, 
      include.DBID.missing=TRUE, 
      ncore=1); 

show(tfs1)

file.bam <- "/local/storage/projects/NHP/AllData/bams/H3_U.fastq.gz.sort.bam"

tfs2 <- tfbs.selectExpressedMotifs(tfs, 
      hg19.twobit, 
      gencode.gtf, 
      file.bam = file.bam, 
      seq.datatype = "RNA-seq",
      pvalue.threshold=0.01, 
      include.DBID.missing=TRUE, 
      ncore=1); 

show(tfs2)

}
\keyword{tfbs object}