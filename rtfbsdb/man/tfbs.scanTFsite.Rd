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
      gen.bed, 
      return.type=c("matches", "posteriors", "maxposterior", "writedb"), 
      file.prefix=NA, 
      usemotifs = NA,
      ncores = 1, 
      threshold = 6, 
      threshold.type = c("score", "fdr"), 
      gc.groups = NA, 
      background.order = 2,
      background.length = 100000)
}

\arguments{
  \item{tfbs}{A tfbs object (\code{"\linkS4class{tfbs}"}) returned by \code{\link{tfbs.createFromCisBP}}, \code{\link{tfbs}}, \code{\link{tfbs.dirs}}. }
  \item{file.twoBit}{String, the file name of genome data( e.g. hg19.2bit or mm10.2bit) }
  \item{gen.bed}{Data frame, bed-formatted loci information with 6 columns}
  \item{return.type}{String, four available values explained in th details(default = "matches")}
  \item{file.prefix}{String, the prefix for outputted file, only used when the return.type is \emph{writedb}}
  \item{usemotifs}{Vector indicating indexes of motif to be used in scanning.}
  \item{ncores}{Number, computing nodes in parallel environment (default = 1).}
  \item{threshold}{Numeric value, if 'score' is specified in \code{threshold.type}, only binding sites with scores above this threshold are returned, 
  if 'fdr' is specified in \code{threshold.type}, only binding sites with FDR (False Discovery Rate) less than this value can be selected. Default value is 6 for 'score' and 0.1 for 'fdr'.} 
  \item{threshold.type}{String value, two options are available. only sites with scores above this threshold are returned, not be used if NA. (default = 'score')}
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
2) \code{bed}: Calling bed-formatted loci(\code{gen.bed}).\cr
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

file.twoBit <- system.file("extdata","hg19.chr19.2bit", package="rtfbsdb")

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, family_name="Ets");

gen.bed <- data.frame(chr="chr19", 
      start=round(runif(10,1000000, 2000000)), 
      stop=0, 
      name="", 
      score=0, 
      strand=".");
gen.bed$stop <- gen.bed$start + 3000;

t1 <- tfbs.scanTFsite( tfs, 
      file.twoBit, 
      gen.bed, 
      file.prefix="test.db",
      ncores = 1);

#show a brief information about the result
t1

#show the summary information in the result
show(t1$summary);

#show the matched TF sites for first motif
show(t1$result[[1]]);

#Output a PDF report for all motifs.
tfbs.reportFinding(tfs, t1, file.pdf="test-rtfbs-scan.pdf", report.title="ELF1");

file.ELF1 <- system.file("extdata","Chipseq-k562-chr19-ELF1.bed", package="rtfbsdb")
gen.bed<- read.table(file.ELF1)

t2 <- tfbs.scanTFsite( tfs, 
      file.twoBit, 
      gen.bed, 
      file.prefix="test.db",
      return.type="writedb",
      ncores = 1);

t2

t3 <- tfbs.scanTFsite( tfs, 
      file.twoBit, 
      gen.bed, 
      return.type="posteriors",
      ncores = 1);

t3

t4 <- tfbs.scanTFsite( tfs, 
      file.twoBit, 
      gen.bed, 
      return.type="maxposterior",
      ncores = 1);	

t4;

dim(t4$result);

}

\keyword{tfbs object}
\keyword{Scanning}
