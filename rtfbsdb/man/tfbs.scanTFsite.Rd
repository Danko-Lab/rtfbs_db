\name{tfbs.scanTFsite}
\alias{tfbs.scanTFsite}
\title{
Find TF sites from genome data within the BED loci
}
\description{
Find TF sites from genome data within the BED loci. Please notice that this package \cr
does not provided genome data such as hg19.2bit, mm10.2bit. 
}
\usage{
tfbs.scanTFsite(tfbs, 
      file.twoBit, 
      tre.bed, 
      return.type=c("matches", "posteriors", "maxposterior", "writedb"), 
      file.prefix=NA, 
      usemotifs = NA,
      ncores = 3, 
      fdr = NA, 
      threshold = 6, 
      gc.groups = NA, 
      background.order = 2,
      background.length = 100000)
}

\arguments{
  \item{tfbs}{A tfbs object (\code{"\linkS4class{tfbs}"}) returned by \code{\link{tfbs.createFromCisBP}}, \code{\link{tfbs}}, \code{\link{tfbs.dirs}}. }
  \item{file.twoBit}{String, the file name of genome data( e.g. hg19.2bit or mm10.2bit) }
  \item{tre.bed}{Data frame, bed-formatted loci information with 6 columns}
  \item{return.type}{String, four available values explained in th details(default = "matches")}
  \item{file.prefix}{String, the prefix for outputted file, only used when the return.type is \emph{writedb}}
  \item{usemotifs}{Vector indicating indexes of motif to be used in scanning.}
  \item{ncores}{Number, computing nodes in parallel environment (default = 3).}
  \item{fdr}{Numeric value between 0 and 1, False Discovery Rate (FDR) of possible binding sites in \code{rtfbs} package, only binding sites with FDR less than this value can be selected. If fdr value is assigned, the \code{threshold} will be ignored.}
  \item{threshold}{Numeric value, only sites with scores above this threshold are returned in \code{rtfbs} package (default = 6).}
  \item{gc.groups}{Numeric value,indicating number of quantiles to group sequences into in \code{rtfbs} package (default = 1).}
  \item{background.order}{Numeric value,indicating the order of Markov model to build in \code{rtfbs} package (default = 2).}
  \item{background.length}{Numeric value, indicating length of the sequence to simulate in \code{rtfbs} package (default = 100000)}
}
\details{

(1) Four options are availabe for the function of tfbs.scanTFsite as follows.\cr
\itemize{
\item{matches:}{ returns all matching TF sites for all motifs.}
\item{writedb:}{ writes a bed file with matches sites. Assumes that sort-bed and starch tools are availiable in $PATH}
\item{posteriors:}{ returns the posteriors at each position in bed-formatted loci.}
\item{maxposterior:}{ returns the max(posterior) in each position in bed-formatted loci.}
}

(2) In order to make the binary file with the parameter of writedb, make sure that starchcat 
and sort-bed command (in  BEDOPS) can be accessed from R environment. If not, please put the 
folder in $PATH.\cr

}

\value{
A list object will be returned with the class name of \code{tfbs.finding}. The object wraps four sub-list as follows:\cr

1) \code{parm}:Calling parameters(\code{fdr}, \code{threshold}), \code{gc.groups}...).\cr
2) \code{bed}: Calling bed-formatted loci(\code{tre.bed}).\cr
3) \code{summary}: A data frame including summrized information about matched TF sites for all motifs.\cr
4) \code{result}: Scanning results which data type is depend on the parameter of \code{return.type}.

The option of \emph{matches} returns a list including the result of every motif, which result is BED style data frame with the following columns.\cr
\item{chrom}{chromosome}
\item{chromStart}{start position}
\item{chromEnd}{chromosome end position}
\item{name}{}
\item{score}{ The score is given by the log likelihood ratio against the Marklov model(backgound).}
\item{strand}{strand}

The option of \emph{writedb} will return a binary BED filename in which store all bed ranges.\cr

The option of \emph{posteriors} will return a list for each motif returned by score.ms function. Scores represent the motif 'match score',
or the product of the probability of observing each base under the motif or background models. Scores are returned under the motif model for all positions in the sequence, on both forward and reverse strands, and under the background
model. \cr

The option of \emph{maxposterior} will return a probability matrix which the row indicates the target loci and the column indicates the motif.\cr

}

\seealso{
\code{\link{print.tfbs.finding}}, \code{\link{summary.tfbs.finding}}, \code{\link{tfbs.reportFinding}}. 
}
\examples{

library(rtfbsdb);
file.twoBit     <- "/local/storage/data/hg19/hg19.2bit"

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, family_name="AP-2");

tre.bed <- data.frame(chr="chr1", 
			start=round(runif(10,1000000, 2000000)), 
			stop=0, 
			name="", 
			score=0, 
			strand=".");
tre.bed$stop <- tre.bed$start + 3000;

t1 <- tfbs.scanTFsite( tfs, 
	file.twoBit, 
	tre.bed, 
	file.prefix="test.db",
	ncores = 1);

#show a brief information about the result
t1

#show the summary information in the result
show(t1$summary);

#show the matched TF sites for first motif
show(t1$result[[1]]);

#Output a PDF report for all motifs.
tfbs.reportFinding(tfs, t1, file.pdf="Test Results");

t2 <- tfbs.scanTFsite( tfs, 
	file.twoBit, 
	tre.bed, 
	file.prefix="test.db",
	return.type="writedb",
	ncores = 1);

t2

t3 <- tfbs.scanTFsite( tfs, 
	file.twoBit, 
	tre.bed, 
	return.type="posteriors",
	ncores = 1);

t3

t4 <- tfbs.scanTFsite( tfs, 
	file.twoBit, 
	tre.bed, 
	return.type="maxposterior",
	ncores = 1);	

t4;

t4$result;

}

\keyword{tfbs object}
\keyword{Scanning}
