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
      file.bigwig.plus, 
      file.bigwig.minus, 
      file.bam=NA, 
      file.twoBit=NA, 
      file.gencode.gtf=NA, 
      seq.datatype=NA,
      pvalue.threshold = 0.05, 
      include.DBID.missing=TRUE, 
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
  \item{pvalue.threshold}{Numeric, indicating . }
  \item{include.DBID.missing}{Logical, indicating whether the TFs without association with GENCODE through the DBID are selected. }
  \item{ncores}{Number, comuputing nodes in parallel environment for gencode data converting.}
}
\details{
}
\value{
}
\examples{
}
\keyword{tfbs object}