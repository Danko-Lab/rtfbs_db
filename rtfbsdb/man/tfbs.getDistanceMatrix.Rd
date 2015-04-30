\name{tfbs.getDistanceMatrix}
\alias{tfbs.getDistanceMatrix}
\title{
Calcuate a distance matrix with Pearson's R values.
}
\description{
Compare any two motifs and return a matrix with Pearson's R values.
}
\usage{
tfbs.getDistanceMatrix(tfbs, ncores = 3, BG = log(c(0.25, 0.25, 0.25, 0.25)))
}
\arguments{
  \item{tfbs}{A tfbs object(\code{"\linkS4class{tfbs}"}).}
  \item{ncores}{Number, the number of cores to use simultaneously.}
  \item{BG}{The log value of probabilities for nucleotide A, C, G and T as Backgroud computing.}
}
\details{
}

\value{
A tfbs object with new distance matrix (@distancematrix).
}
\seealso{
}
\examples{

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, family_name="AP-2");
tfs0 <- tfbs.getDistanceMatrix(tfs, ncores=1);

}
\keyword{tfbs object}
\keyword{Distance matrix}
