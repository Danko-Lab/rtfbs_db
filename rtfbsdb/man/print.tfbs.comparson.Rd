\name{print.tfbs.comparson}
\alias{print.tfbs.comparson}
\title{
Print the brief comparson results
}
\description{
Print the brief comparson results.
}
\usage{
\method{print}{tfbs.comparson}(x, \dots)
}
\arguments{
  \item{x}{The result obtained by \code{\link{tfbs.compareTFsite}}.}
  \item{\dots}{Additional arguments affecting the print produced.}
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
