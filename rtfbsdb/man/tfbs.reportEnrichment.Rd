\name{tfbs.reportEnrichment}
\alias{tfbs.reportEnrichment}
\title{
Output report for enrichment results.
}
\description{
Output enrichment results to a PDF report which includes motif names, counts of TF site, p-value, enrichment ratio and motif logos. 
}
\usage{
tfbs.reportEnrichment(tfbs, r.comp, 
      file.pdf = NA, 
      report.size = "letter", 
      report.title = "", 
      enrichment.type = c ("both", "enriched", "depleted"),
      sig.only = TRUE, 
      pv.threshold = 0.05, 
      pv.adj = NA)
}
\arguments{
  \item{tfbs}{A tfbs object, see also \code{"\linkS4class{tfbs}"} }
  \item{r.comp}{A result object from the function of \code{\link{tfbs.enrichmentTest}}}
  \item{file.pdf}{String, the file name of PDF report.}
  \item{report.size}{String, the page size ( default="letter")}
  \item{report.title}{String, the report title. }
  \item{enrichment.type}{String, three values are available for significant motifs to be printed out.( default="both").}
  \item{sig.only}{String, indicating whether only significant motifs are outputted or not.( default=TRUE).}
  \item{pv.threshold}{Numeric value,indicating whether the different threshold of p-value is applied to select the significant motifs. }
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
\code{\link{tfbs.enrichmentTest}}, \code{\link{summary.tfbs.enrichment}}.
}
\examples{
# see examples in tfbs.enrichmentTest
}

\keyword{ tfbs object }
\keyword{ Enrichment }
