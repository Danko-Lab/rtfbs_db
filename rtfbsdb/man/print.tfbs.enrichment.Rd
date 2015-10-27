\name{print.tfbs.enrichment}
\alias{print.tfbs.enrichment}
\title{
Print the brief enrichment results
}
\description{
Print the brief enrichment results.
}
\usage{
\method{print}{tfbs.enrichment}(x, \dots, pv.threshold=0.05, pv.adj=NA )
}
\arguments{
  \item{x}{The result obtained by \code{\link{tfbs.enrichmentTest}}.}
  \item{\dots}{Additional arguments affecting the print produced.}
  \item{pv.threshold}{Numeric value,indicating whether the different cutoff of p-value is applied to select the significant motifs. }
  \item{pv.adj}{String, P-values correct method for \code{p.adjust} function. The available values are "holm", \cr
  "hochberg", "hommel", "bonferroni", "BH", "BY","fdr" or "none".\cr
  (default="bonferroni").}
}
\details{
This command shows the calling parameters and significant motifs from the result object.\cr
The significant motifs are selected by the corrected p-value cutoff(0.05) and at most 20 significant motifs are listed.
The adjust method of p-value is defined in the calling function.
}
\value{
No return values.
}

\seealso{
See also as \code{\link{tfbs.enrichmentTest}}.
}
\examples{
#See example in tfbs.enrichmentTest
}
\keyword{print}
\keyword{Enrichment}
