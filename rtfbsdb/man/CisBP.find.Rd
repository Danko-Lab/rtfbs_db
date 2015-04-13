\name{CisBP.find}
\alias{CisBP.find}
\title{
Query the CisBP dataset.
}
\description{
Build a tfbs object by querying the meta file of CisBP dataset and subseting the results.
}
\usage{
CisBP.find(cisbp.db, 
      tf_name = NULL, 
      tf_status = NULL, 
      family_name = NULL, 
      motif_type = NULL, 
      msource_id = NULL, 
      motif_info_type = 1)
}
\arguments{
  \item{cisbp.db}{A CisBP object, including the file of TF_Information.txt.}
  \item{tf_name}{String, indicting the TF_name field will be used to select motifs.}
  \item{tf_status}{String, indicting the TF_Status field will be used to select motifs. }
  \item{family_name}{String, indicting the Family_Name field will be used to select motifs. }
  \item{motif_type}{String, indicting the Motif_Type field will be used to select motifs.  }
  \item{msource_id}{String, indicting the MSource_Identifier field will be used to select motifs. }
  \item{motif_info_type}{Number,indicating which TF meta file will be used. Available values are 1 \cr
  for TF_Information.txt, 2 for TF_Information_all_motifs.txt and 3 for \cr
  F_Information_all_motifs_plus.txt.}
}
\details{
}
\value{
A tfbs object is returned with PWM matrices, see Also as \code{\link{tfbs}}
}
\seealso{
See Also as \code{\link{tfbs}}
}
\examples{

# Load the internal CisBP dataset
db_human <- CisBP.extdata("Homo_sapiens");

# Load all motifs and return a tfbs object.
tfs0 <- CisBP.find(db_human); 

# Query the motifs by the conditins and return a tfbs object 
tfs1 <- CisBP.find(db_human, family_name="Homeodomain", tf_status="D",  
      motif_type="ChIP-seq", msource_id= "MS01_1.01", motif_info_type=1 ); 

# Query the motifs by the conditins and return a tfbs object 
tfs2 <- CisBP.find(db_human, family_name="Homeodomain", tf_status="D",  motif_info_type=1 ); 

# Query the motifs by the conditins and return a tfbs object 
tfs3 <- CisBP.find(db_human, motif_type="ChIP-seq", motif_info_type=1 ); 

# Query the motifs by the conditins and return a tfbs object 
tfs4 <- CisBP.find(db_human, motif_info_type=2); 

}
\keyword{CisBP object}
\keyword{tfbs object}