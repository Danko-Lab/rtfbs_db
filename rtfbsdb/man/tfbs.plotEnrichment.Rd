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
      enrichment.type = c ("both", "enriched", "depleted"), 
      plot.type=c("nonpolar", "polar"), 
      options=list() )
}
\arguments{
  \item{tfbs}{A tfbs object, see also \code{"\linkS4class{tfbs}"} }
  \item{r.comp}{A result object from the function of \code{\link{tfbs.enrichmentTest}}}
  \item{file.pdf}{String, the file name of PDF report.}
  \item{enrichment.type}{String, three values are available for significant motifs to be drawn out.(default="both").}
  \item{plot.type}{String, indicating whether the depleted motifs are shown at the negative side.(default="nonpolar") }
  \item{options}{List, containing the pre-defined parameters to control the plot, see details.}
}
\details{
The function draws a QQ-like plot which all motifs are ordered by its p-values along X axis and -log10(pvalue) are drawn along Y axis. 
If plot type is 'ploar', depleted motifs are seperated from enriched motifs and drawn at the negative side of Y axis.
The motif color is determined by the enrichment ratio.\cr

The parameter \code{options} has the following options:\cr
\preformatted{
options=list(
	title  = "",
	width = 7,
	height = 7,
	y.max  = NULL,   
	xlab   = "Order",
	ylab   = "-log10(p-value)",
	top.motif.labels = 5,
	bottom.motif.labels = 5,
	color.scheme = 2,
	abline = NULL,
	zoom.tick = 1,
	zoom.label = 1,
	zoom.motif.logo = 1,
	zoom.legend.label=1,
	zoom.motif.label = 1 );
}

\itemize{
\item{\code{title}:                 }{ String, the plot title. }
\item{\code{width}:                 }{ Number, page width, default is 7 inches. }
\item{\code{height}:                }{ Number, page height, default is 7 inches. }
\item{\code{y.max}:                 }{ Number, the maximum value for Y axis. The points larger than \cr  this value will be drawn at the border.}
\item{\code{xlab}:                  }{ String, the label for X axis. }
\item{\code{ylab}:                  }{ String, the label for Y axis. }
\item{\code{top.motif.labels}:      }{ String or number, the labels for top motifs. }
\item{\code{bottom.motif.labels}:   }{ String or number, the labels for top motifs in the depleted side.}
\item{\code{color.scheme}:          }{ Number, 1 or 2 are available, color scheme for the points based on enrichment ratios. (default=2).}
\item{\code{abline}:                }{ Number, indicating to draw horizontal lines for significant p-values.}
\item{\code{zoom.tick}:             }{ Number, the ratio to zoom in or out the ticks on the X and Y axis. (default=1).}
\item{\code{zoom.label}:            }{ Number, the ratio to zoom in or out the labels on the X and Y axis. (default=1).}
\item{\code{zoom.motif.logo}:       }{ Number, the ratio to zoom in or out the motif logos located  at the top or bottom. (default=1)}
\item{\code{zoom.motif.label}:      }{ Number, the ratio to zoom in or out the motif labels located  at the top or bottom. (default=1)}
\item{\code{zoom.legend.label}:     }{ Number, the ratio to zoom in or out the legend label. (default=1)}
}

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
