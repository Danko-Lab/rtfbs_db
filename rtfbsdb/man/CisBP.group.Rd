\name{CisBP.group}
\alias{CisBP.group}
\title{
Summarize the motif number.
}
\description{
Get the statistical summary by grouping the field in the CisBP dataset.
}
\usage{
CisBP.group(cisbp.db, 
      group_by=c("tf_name", "tf_species", "tf_status", "family_name", 
            "motif_type", "msource_id"), 
      motif_info_type=1)
}
\arguments{
  \item{cisbp.db}{ A CisBP object including the TF_Information.txt.}
  \item{group_by}{ String, indicating which field will be used to group values. Available values are tf_name, tf_species, tf_status, family_name, motif_type and msource_id.}
  \item{motif_info_type}{ Number, indicating which TF meta file will be used. Available values are 1 for TF_Information.txt, 2 for TF_Information_all_motifs.txt and 3 for F_Information_all_motifs_plus.txt.}
}
\details{
Three TF information files in CisBP dataset.\cr\cr
1: TF_Information.txt : (direct motifs) or (no direct but inferred motifs with 90\%)\cr
2: TF_Information_all_motifs.txt: (direct motifs) and (inferred motifs above the threshold)\cr
3: F_Information_all_motifs_plus.txt: All motifs\cr
}
\value{
  A data frame returned includes two columns
  \item{group_by}{Values of grouping field}
  \item{number}{Counts of group value  }
}
\seealso{
See Also as \code{\link{CisBP.find}}
}
\examples{

# Load the internal CisBP dataset
db_human <- CisBP.extdata("Homo_sapiens");

# Group the motif count by the column of family_name in TF_Information.txt
gr1 <- CisBP.group(db_human, group_by="family_name", motif_info_type=1 ); 

# Group the motif count by the column of tf_status in TF_Information.txt
gr2 <- CisBP.group(db_human, group_by="tf_status", motif_info_type=1 ); 

# Group the motif count by the column of tf_status in TF_Information_all_motifs.txt
gr3 <- CisBP.group(db_human, group_by="tf_status", motif_info_type=2); 

# Group the motif count by the column of tf_status in F_Information_all_motifs_plus.txt
gr4 <- CisBP.group(db_human, group_by="tf_status", motif_info_type=3); 

}

\keyword{CisBP object}
