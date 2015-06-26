\name{tfbs.reportComparson}
\alias{tfbs.reportComparson}
\title{
Output report for comparson results.
}
\description{
Output comparson results to a PDF report which includes motif names, counts of TF site, p-value, enrichment ratio and motif logos. 
}
\usage{
tfbs.reportComparson(tfbs, r.comp, 
      file.pdf = NA, 
      report.size = "letter", 
      report.title = "", 
      sig.only = TRUE, 
      pv.cutoff = NA, 
      pv.adj = NA)
}
\arguments{
  \item{tfbs}{A tfbs object, see also \code{"\linkS4class{tfbs}"} }
  \item{r.comp}{A result object from the function of \code{\link{tfbs.compareTFsite}}}
  \item{file.pdf}{String, the file name of PDF report.}
  \item{report.size}{String, the page size ( default="letter")}
  \item{report.title}{String, the report title. }
  \item{sig.only}{String, indicating whether only significant motifs are outputted or not.( default=TRUE).}
  \item{pv.cutoff}{Numeric value,indicating whether the different cutoff of p-value is applied to select the significant motifs. }
  \item{pv.adj}{String,indicating whether the different correction metod of p-value is applied to select the significant motifs. }
}
\details{
The table with 7 columns is outputted into a PDF report within letter size. \cr
Two color bars are used to display p-values and enrichment ratios. Motif logos are shown visually in each row.
}
\value{
No return values.
}
\seealso{
\code{\link{tfbs.compareTFsite}}, \code{\link{summary.tfbs.comparson}}.
}
\examples{
# see examples in tfbs.compareTFsite
}

\keyword{ tfbs object }
\keyword{ Comparing }
