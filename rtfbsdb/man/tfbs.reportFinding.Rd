\name{tfbs.reportFinding}
\alias{tfbs.reportFinding}
\title{
Make report for scanning results.
}
\description{
Output a PDF report includes motif names, counts of TF site and motif logos. 
}
\usage{
tfbs.reportFinding(tfbs, 
      r.scan, 
      file.pdf = NA, 
      report.size = "letter", 
      report.title = "")
}
\arguments{
  \item{tfbs}{A tfbs object, see also \code{"\linkS4class{tfbs}"} }
  \item{r.scan}{A result object from the function of \code{\link{tfbs.scanTFsite}}}
  \item{file.pdf}{String, the file name of PDF report.}
  \item{report.size}{String, the page size ( default="letter")}
  \item{report.title}{String, the report title. }

}
\details{
The table with 4 columns is outputted into a PDF report within letter size. \cr
Motif logos are shown visually in each row.
}
\value{
No return values.
}

\seealso{
\code{\link{tfbs.scanTFsite}}, \code{\link{print.tfbs.finding}}
}

\examples{
#See example in tfbs.scanTFsite 
}
\keyword{ tfbs object }
\keyword{Scanning}
