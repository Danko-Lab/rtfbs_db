\name{summary.tfbs.comparson}
\alias{summary.tfbs.comparson}
\title{
Summarize the comparson result
}
\description{
Return the significant motifs based on the adjust p-values using multiple comparisons.
}
\usage{
\method{summary}{tfbs.comparson}(object, pv.cutoff = 0.05, pv.adj = NA, \dots)
}

\arguments{
  \item{object}{The result obtained by \code{\link{tfbs.compareTFsite}}.}
  \item{pv.cutoff}{The p-value cutoff for significant motifs.}
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
\item{Nneg }{Read count in negative loci.}
\item{pv.adj}{p-value}
\item{es.ratio}{The ratio of read counts between positive loci and negative loci. }
}
\seealso{
See also as \code{\link{tfbs.compareTFsite}}.
}
\keyword{summary}
\keyword{comparson}
