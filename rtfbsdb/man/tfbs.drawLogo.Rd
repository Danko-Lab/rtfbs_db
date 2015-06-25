\name{tfbs.drawLogo}
\alias{tfbs.drawLogo}
\title{
Draw single motif logo.
}
\description{
Draw the motif logos in two models, 1 logo within a page or 1 group within one page.
}
\usage{
tfbs.drawLogo(tfbs, file.pdf = NULL, index = NULL, tf_id = NULL,
                 motif_id = NULL, tf_name = NULL, family_name = NULL,
                 tf_status = NULL, groupby = NULL)
}
\arguments{
  \item{tfbs}{A tfbs object(\code{"\linkS4class{tfbs}"})}
  \item{file.pdf}{String, the file name of PDF report.}
  \item{index}{Vector of number, indicating the motif index.}
  \item{tf_id}{Vector of string, indicating the TF_ID string, TF_ID is one motif attribute in TF_Information.txt. (Default=NULL).}
  \item{motif_id}{Vector of string, indicating the Motif_ID string, Motif_ID is one motif attribute in TF_Information.txt. (Default=NULL).}
  \item{tf_name}{Vector of string, indicating the TF_Name string, TF_Name is one motif attribute in TF_Information.txt. (Default=NULL).}
  \item{family_name}{Vector of string, indicating Family_Name string, Family_Name is one motif attribute in TF_Information.txt. (Default=NULL).}
  \item{tf_status}{String, indicating the TF_status value, TF_status is one motif attribute in TF_Information.txt. (Default=NULL).}
  \item{groupby}{String, indicating the group field is applied to print the motif, each group is printed in one page, the available values are NA, "Family_Name", "TF_Name", "TF_Status" or "Motif_Type". (Default=NA). }
}
\details{

Multiple selection is provided for outputting logos. The selected motifs by each criteria will be combined into one set.\cr 

Draw the motif logos in two models: \cr
(1) 1 logo within a page 
(2) 1 group within one page.The motif logos are splitted if motif count is greater than 10.


}
\value{
No return values.
}

\seealso{
See Also as \code{"\linkS4class{tfbs}"}
}
\examples{
 
 db <- CisBP.extdata("Homo_sapiens");

 tfs <- tfbs.createFromCisBP(db);
 
 motif_id   <- c( "M5604_1.01", "M5441_1.01", "M5162_1.01", "M5352_1.01");
 tf_id      <- c( "T093250_1.01", "T093251_1.01","T093252_1.01","T093253_1.01");
 family_name<- c( "p53", "Homeodomain", "Paired box", "Pipsqueak");


 #Draw 10 motif logos from first one. 
 tfbs.drawLogo(tfs, file.pdf="tfbs.drawLogo1.pdf", index=c(1:10) );
 
 #Draw logos for specified Motif_ID, or TF_ID, or TF_Name, or Family_Name 
 tfbs.drawLogo(tfs, file.pdf="tfbs.drawLogo2.pdf", 
 			motif_id=motif_id, 
 			tf_id=tf_id, 
 			tf_name="AP-2", 
 			family_name=family_name, 
 			groupby="TF_Status");
 
 #Draw logos for specified TF_Status
 tfbs.drawLogo(tfs, file.pdf="tfbs.drawLogo3.pdf", tf_status="D", groupby="TF_Status");

 #unlink("tfbs.drawLogo1.pdf");
 #unlink("tfbs.drawLogo2.pdf");
 #unlink("tfbs.drawLogo3.pdf");
 
}
\keyword{Logo}
\keyword{tfbs object}
