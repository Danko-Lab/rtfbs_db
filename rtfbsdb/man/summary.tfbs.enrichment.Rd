\name{summary.tfbs.enrichment}
\alias{summary.tfbs.enrichment}
\title{
Summarize the enrichment result
}
\description{
Return the significant motifs based on the adjust p-values using multiple comparisons.
}
\usage{
\method{summary}{tfbs.enrichment}(object, pv.threshold = 0.05, pv.adj = NA, \dots)
}

\arguments{
  \item{object}{The result obtained by \code{\link{tfbs.enrichmentTest}}.}
  \item{pv.threshold}{The p-value threshold for significant motifs.}
  \item{pv.adj}{P-values adjust method for \code{p.adjust} function. The available values are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr" or "none".}
  \item{\dots}{Additional arguments affecting the summary produced.}
}
\details{
A data frame with 6 colums is returned.
}
\value{
The results is a data frame including 6 columns,
\item{motif.id }{Motif ID}
\item{tf.name }{TF Name}
\item{Npos }{Read count in positive loci.}
\item{expected }{Read count in negative loci.}
\item{fe.ratio}{The ratio of read counts between positive loci and negative loci. }
\item{starch}{Cpmporessed Bed filename}
\item{pvalue}{p-value}
\item{pv.adj}{adjusted p-value by multiple comparson method.}
}
\seealso{
See also as \code{\link{tfbs.enrichmentTest}}.
}
\keyword{summary}
\keyword{Enrichment}
