\name{tfbs.plotEnrichment}
\alias{tfbs.plotEnrichment}
\title{
Output report for enrichment results.
}
\description{
Output enrichment results to a PDF report which includes motif names, counts of TF site, p-value, enrichment ratio and motif logos. 
}
\usage{
tfbs.plotEnrichment(tfbs, r.comp, 
      file.pdf, 
      plot.title = "", 
      top.motif.labels = 5,
      bottom.motif.labels = 5, 
      xlab="Order", 
      ylab="-log10(p-value)", 
      y.max=NULL,
      enrichment.type = c ("both", "enriched", "depleted"), 
      plot.type=c("nonpolar", "polar"),
      color.scheme=2)
}
\arguments{
  \item{tfbs}{A tfbs object, see also \code{"\linkS4class{tfbs}"} }
  \item{r.comp}{A result object from the function of \code{\link{tfbs.enrichmentTest}}}
  \item{file.pdf}{String, the file name of PDF report.}
  \item{plot.title}{String, the plot title. }
  \item{top.motif.labels}{String or number, the labels for top motifs.}
  \item{bottom.motif.labels}{String or number, the labels for top motifs in the depleted side.}
  \item{xlab}{String, the label for X axis. }
  \item{ylab}{String, the label for Y axis. }
  \item{y.max}{Number, the maximum value for Y axis. The points larger than this value will be drawn at the border. }
  \item{enrichment.type}{String, three values are available for significant motifs to be drawn out.(default="both").}
  \item{plot.type}{String, indicating whether the depleted motifs are shown at the negative side.(default="nonpolar") }
  \item{color.scheme}{Number, 1 or 2 are available, color scheme for the points based on enrichment ratios. (default=2) }
}
\details{
The function draws a QQ-like plot which all motifs are ordered by its p-values along X axis and -log10(pvalue) are drawn along Y axis. 
If plot type is 'ploar', depleted motifs are seperated from enriched motifs and drawn at the negative side of Y axis.
The motif color is determined by the enrichment ratio.\cr

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
