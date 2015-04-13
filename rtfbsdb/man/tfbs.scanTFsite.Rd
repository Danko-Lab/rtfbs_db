\name{tfbs.scanTFsite}
\alias{tfbs.scanTFsite}
\title{
Find TF sites from genome data within the BED ranges
}
\description{
Find TF sites from genome data within the BED ranges. Please notice that this package \cr
does not provided genome data such as hg19.2bit, mm10.2bit. 
}
\usage{
tfbs.scanTFsite(tfbs, 
      file.twoBit, 
      dnase.peaks.bed, 
      file.prefix="scan.db", 
      usemotifs=NA,
      ncores=3, 
      return.type=c("matches", "posteriors", "maxposterior", "writedb"), 
      threshold=6, ...)
}

\arguments{
  \item{tfbs}{A tfbs object returned by \code{CisBP.find}, \code{tfbs}, \code{tfbs.dirs}. }
  \item{file.twoBit}{String, the file name of genome data( hg19.2bit or mm10.2bit) }
  \item{dnase.peaks.bed}{Data frame, bed-formatted peak information}
  \item{file.prefix}{String, the prefix for outputted file, only used when the return.type is \emph{writedb} }
  \item{usemotifs}{Vector indtciting indexes of motif to be used in scanning.}
  \item{ncores}{Number, comupting nodes in parallel environment.}
  \item{return.type}{String, four available values explained in th details }
  \item{threshold}{Numeric value, only sites with scores above this threshold are returned (default = 6)}
  \item{\dots}{Any parameters used in the function of score.ms.}
}
\details{

(1) Four options are availabe for the function of tfbs.scanTFsite as follows.\cr
\itemize{
\item{matches}{returns all matching motifs.}
\item{writedb}{writes a bed file with matches.  Assuems that sort-bed and starch tools are availiable in $PATH}
\item{posteriors}{returns the posteriors at each position.}
\item{maxposterior}{returns the max(posterior) in each dnase-1 peak.}
}

(2) In order to make the binary file with the parameter of writedb, make sure that starchcat 
and sort-bed command (in  BEDOPS) can be accessed from R environment. If not, please put the 
folder in $PATH.\cr

}

\value{
The option of \emph{matches} returns a list including the result of every motif, which result is BED style data frame with the following columns.\cr
\item{chrom}{chrmosome}
\item{chromStart}{start position}
\item{chromEnd}{chrmosomeend position}
\item{name}{}
\item{score}{}
\item{strand}{strand}

The option of \emph{writedb} will return a binary BED filename in which store all bed ranges.\cr

The option of \emph{posteriors} will return a list for each motif returned by score.ms function. Scores represent the motif ’match score’,
or the product of the probability of observing each base under the motif or background models.  Scores are returned under the motif model for all positions in the sequence, on both forward and reverse strands, and under the background
model. \cr

The option of \emph{maxposterior} will return a probability matrix which the row inditcates the range of dnase peak and the column indicates the motif.\cr

}


\seealso{
}
\examples{

library(rtfbsdb);

file.dREG.H.change.bed <- "/home/zw355/src/rtfbs_db/rtfbsdb/test/dREG.H.change.bed"
file.twoBit     <- "/local/storage/data/hg19/hg19.2bit"

db <- CisBP.extdata("Homo_sapiens");
tfs <- CisBP.find(db, family_name="AP-2");

dREG_H_change_bed <- read.table(file.dREG.H.change.bed, header=FALSE);

t0 <- tfbs.scanTFsite( tfs, 
	file.twoBit, 
	dREG_H_change_bed, 
	file.prefix="test.db",
	ncores = 1);

str(t0[[1]])

t1 <- tfbs.scanTFsite( tfs, 
	file.twoBit, 
	dREG_H_change_bed, 
	file.prefix="test.db",
	return.type="writedb",
	ncores = 1);

t1

t2 <- tfbs.scanTFsite( tfs, 
	file.twoBit, 
	dREG_H_change_bed, 
	file.prefix="test.db",
	return.type="posteriors",
	ncores = 1);

str(t2[[1]]);

t3 <- tfbs.scanTFsite( tfs, 
	file.twoBit, 
	dREG_H_change_bed, 
	file.prefix="test.db",
	return.type="maxposterior",
	ncores = 1);	

str(t3);

}

\keyword{tfbs object}
\keyword{Scanning}
