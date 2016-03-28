\name{tfbs.getExpression}
\alias{tfbs.getExpression}
\title{
Estimate gene expression of target TF.
}
\description{
Gets expression level of target TF.\cr
USE extra_info$DBID to find gene information encoded by GENCODE\cr
}
\usage{
tfbs.getExpression(tfbs, 
      file.twoBit, 
      file.gencode.gtf,
      file.bigwig.plus=NA, 
      file.bigwig.minus=NA, 
      file.bam=NA, 
      seq.datatype = c("GRO-seq", "PRO-seq", "RNA-seq"),
      ncores =1 )
}
\arguments{
  \item{tfbs}{A tfbs object(\code{"\linkS4class{tfbs}"}).}
  \item{file.twoBit}{String, indicating the binary data of sequence. (e.g. hg19.2bit, mm10.2bit)}
  \item{file.gencode.gtf}{Gencode RDATA file encoded by ths package.}
  \item{file.bigwig.plus}{String, indicating bigwig file for strand plus(+) if \code{seq.datatype} is GRO-seq or PRO-seq.}
  \item{file.bigwig.minus}{String, indicating bigwig file for strand minus(-) if \code{seq.datatype} is GRO-seq or PRO-seq.}
  \item{file.bam}{String, indicating BAM file for rna reads if \code{seq.datatype} is RNA-seq.}
  \item{seq.datatype}{String,indicating which kind of seq data is applied to this function, three values are available: GRO-seq, PRO-seq and RNA-seq. (Default=GRO-seq)}
  \item{ncores}{Number, comupting nodes in parallel environment.}
}

\details{

1) For each motif, the occurance ranges can be queried by the gene ID in the GENCODE 
database( for human, gencode.v19.annotation.gtf, for mouse: gencode.vM3.annotation.gtf). 
After the searching, one range obtianed from the merge of the multiple ranges will be 
used to detect the reads count in the specified  bigwig files(including plus and minus). 
The probability of each motif can be calcuated by the reads count and lambda.\cr

The lambda is determined by the following formulation:\cr

\code{r.lambda = 0.04 * sum(reads_in_all_chromosomes)/10751533/1000.}\cr

2) This function heavily relys on some Unix commands to operate bed data and gencode data. Please make sure
the following commands work normally in R terminal.\cr

\code{awk}, \code{zcat}, \code{sort-bed}, \code{twoBitInfo}, \code{bedtools}, \code{samtools}.\cr

The function \code{\link{Sys.which}} and \code{\link{system}} help you to locate these commands and test its availability.\cr

}

\value{
A tbfs object with new expression data frame including the follwing columns:
\item{Motif_ID}{ Motif_ID from CisBP dataset or other data source.}
\item{DBID}{ DBID from CisBP dataset or other data source.}
\item{chr}{String, chromosome name.}
\item{start}{Integer, start postion in which gene ID can be detected.}
\item{end}{Integer, end postion in which gene ID can be detected.}
\item{strand}{String, + or -, indicating the strand direction.}
\item{bed_length}{Integer, the length of range which gene ID can be detected.}
\item{reads}{The reads number queried by BigWig function from the bigwig files( plus and minus) }
\item{lambda}{The lambda parameter in poison distribution. }
\item{reads.RPKM}{The RPKM value of \code{reads} column.}
\item{lambda.RPKM}{The RPKM value of \code{lambda} column.}
\item{prob}{The probability calculated based on Poisson distribution.}
}

\seealso{
See Also as \code{"\linkS4class{tfbs}"}
}
\examples{

# Load the internal CisBP data set
db.human <- CisBP.extdata("Homo_sapiens");

# Create a tfbs object by querying the meta file of CisBP dataset. 
tfs <- tfbs.createFromCisBP(db.human, motif_type="ChIP-seq", 
      tf.information.type=1 );

file.bigwig.minus <- system.file("extdata","GSM1480327_K562_PROseq_chr19_minus.bw", 
      package="rtfbsdb")
file.bigwig.plus <- system.file("extdata","GSM1480327_K562_PROseq_chr19_plus.bw", 
      package="rtfbsdb")
hg19.twobit <- system.file("extdata","hg19.chr19.2bit", package="rtfbsdb")
gencode.gtf <- system.file("extdata","gencode.v19.annotation.chr19.gtf.gz", 
      package="rtfbsdb")

tfs <- tfbs.getExpression(tfs, hg19.twobit, gencode.gtf, 
      file.bigwig.plus, file.bigwig.minus, ncore=1); 

}
\keyword{tfbs object}
\keyword{Gene expression}
