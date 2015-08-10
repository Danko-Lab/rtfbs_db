\name{print.tfbs.comparson}
\alias{print.tfbs.comparson}
\title{
Print the brief comparson results
}
\description{
Print the brief comparson results.
}
\usage{
\method{print}{tfbs.comparson}(x, \dots, pv.cutoff=0.05, pv.adj=NA )
}
\arguments{
  \item{x}{The result obtained by \code{\link{tfbs.compareTFsite}}.}
  \item{\dots}{Additional arguments affecting the print produced.}
  \item{pv.cutoff}{Numeric value,indicating whether the different cutoff of p-value is applied to select the significant motifs. }
  \item{pv.adj}{String, P-values correct method for \code{p.adjust} function. The available values are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr" or "none".(default="bonferroni").}
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
See also as \code{\link{tfbs.compareTFsite}}.
}
\examples{
#See example in tfbs.compareTFsite
}
\keyword{print}
\keyword{Comparing}
