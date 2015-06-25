\name{summary.tfbs.finding}
\alias{summary.tfbs.finding}
\title{
Summarize scanning results.
}
\description{
Return a data frame with summarized TF sites for every motif if the calling parameter is "matches".
}
\usage{
\method{summary}{tfbs.finding}(object, \dots)
}
\arguments{
  \item{object}{The result obtained by \code{\link{tfbs.scanTFsite}}. }
  \item{\dots}{Additional arguments affecting the summary produced.}
}
\details{
\code{summary} in class of \code{tfbs.finding} is returned.
}
\value{
This function will return a data frame with summarized TF sites for every motif if the calling parameter is "matches", otherwise,  NULL will be returned.
}
\seealso{
See Also as \code{\link{tfbs.scanTFsite}}
}
\keyword{ summary }
\keyword{Scanning}
