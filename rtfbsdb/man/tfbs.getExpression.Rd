\name{tfbs.getExpression}
\alias{tfbs.getExpression}
\title{
Estimate gene expression of target TF.
}
\description{
Gets expression level of target TF.\cr
USE extra_info$DBID to find gene information encoded by GENCODE V21\cr
}
\usage{
tfbs.getExpression(tfbs, 
      file.bigwig.plus, file.bigwig.minus, 
      file.twoBit=NA, 
      gencode.ext.rdata=NA)
}
\arguments{
  \item{tfbs}{A tfbs object outputted by CisBP.extdata(), CisBP.dowmload(), CisBP.zipload().}
  \item{file.bigwig.plus}{String, indicating bigwig file for strand plus(+).}
  \item{file.bigwig.minus}{String, indicating bigwig file for strand minus(-).}
  \item{file.twoBit}{String, indicating the binary data of sequence.}
  \item{gencode.ext.rdata}{Gencode RDATA file encoded by ths package.}
}

\details{
For each motif, the occurance ranges can be queried by the gene ID in the GENCODE 
database( for human, gencode.v21.annotation.gtf, for mouse: gencode.vM3.annotation.gtf). 
After the searching, one range obtianed from the merge of the multiple ranges will be 
used to detect the reads count in the specified  bigwig files(including plus and minus). 
The probability of each motif can be calcuated by the reads count and lambda.\cr

The lambda is determined by the following formulation:\cr

\code{r.lambda = 0.04 * sum(reads_in_all_chromosomes)/10751533/1000.}\cr

The dataset of GENECODE v21 (human) and vM3 (mouse) have been compiled into RDATA file 
and attached in this package.\cr

The \code{gencode_transcript_ext} object can be accessed after the following command is 
executed successfully.\cr

\code{load( system.file("extdata", "gencode_human21_transcript_ext.rdata", 
      package="rtfbsdb"), environment() );}\cr
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
\item{lambda}{The lammbda parameter in poison distribution. }
\item{prob}{The probability calculated based on poison distribution.}
}

\seealso{
See Also as \code{\link{tfbs}}, ~~~
}
\examples{

# Load the internal CisBP data set
db.human <- CisBP.extdata("Homo_sapiens");

# Create a tfbs object by querying the meta file of CisBP dataset. 
tfs <- CisBP.find(db.human, motif_type="ChIP-seq", motif_info_type=1 );

if(0)
{
	file.bigwig.plus <- "/work/03350/tg826494/hg19/GSM1480327_K562_PROseq_plus.bw";
	file.bigwig.minus <- "/work/03350/tg826494/hg19/GSM1480327_K562_PROseq_minus.bw";

	tfs <- tfbs.getExpression(tfs, file.bigwig.plus, file.bigwig.minus );
}

}
\keyword{tfbs object}
\keyword{Gene expression}
