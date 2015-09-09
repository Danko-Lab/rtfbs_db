\name{CisBP.getTFinformation}
\alias{CisBP.getTFinformation}
\title{Get TF information with PWM status }
\description{ Get TF information with PWM status }
\usage{
CisBP.getTFinformation(cisbp.db, tf.information.type = NA)
}
\arguments{
  \item{cisbp.db}{ A CisBP object (\code{"\linkS4class{CisBP.db}"}) including the TF_Information.txt.}
  \item{tf.information.type}{Number, indicating which TF meta file will be used. Available values are 1 for \cr
  TF_Information.txt, 2 for TF_Information_all_motifs.txt and 3 \cr
  for F_Information_all_motifs_plus.txt.}
}
\details{
Three TF information files in CisBP dataset.\cr\cr
1: TF_Information.txt : (direct motifs) or (no direct but inferred motifs with 90\%)\cr
2: TF_Information_all_motifs.txt: (direct motifs) and (inferred motifs above the threshold)\cr
3: F_Information_all_motifs_plus.txt: All motifs\cr

The following parts are copied from RAEDME.txt in zipped CisBP data file.\cr

TF_Information.txt, TF_Information_all_motifs.txt, TF_Information_all_motifs_plus.txt - These files contain information on the TFs. \cr\cr

'TF_Information.txt' contains, for each TF, all directly determined motifs (see below). If a TF does not have a directly determined motif, this file will also include its best inferred motif.  'Best' is defined as the motif(s) obtained from the most similar TF (based on the %ID in the amino acids of its TF) that has a directly determined motif.\cr\cr

'TF_Information_all_motifs.txt' is a superset of 'TF_Information.txt'.  It also includes any motif that can be inferred for a given TF, given the TF family-specific threshold.  For example, if a TF has a directly determined motif, and two TFs with motifs with 90% and 80% identical DBDs to it, TF_Information.txt will only include the direct motif, but 
TF_Information_all_motifs.txt will include all three motifs.  Likewise, if a TF does not have a direct motif, but has two TFs with 90% and 80% identical DBDs to it, TF_Information.txt will only include the motif from the 90% indentical TF, while TF_Information_all_motifs.txt would include both.\cr\cr

'TF_Information_all_motifs_plus.txt' is a superset of the other two files.  It contains all motifs for a given TF, which includes all direct motifs, and all inferred motifs above the threshold.\cr

}
\value{
  A data frame returned with the status indicating PWM data is existing or not 
  \item{TF_ID}{Internal CisBP ID for the TF.  Each gene has a unique TF_ID}
  \item{Family_ID}{Internal CisBP ID for the TF family.  A family is the unique set of DNA binding domains (DBDs) present in the protein.}
  \item{TSource_ID}{Internal CisBP ID for the source of the TF (i.e. where its genome sequence was obtained).}
  \item{Motif_ID}{Internal CisBP ID for the associated motif.}
  \item{MSource_ID}{Internal CisBP ID for the source of the motif (i.e. which database or study it came from)}
  \item{DBID}{External ID of the RBP (e.g., Ensembl ID)}
  \item{TF_Name}{Name of the TF}
  \item{TF_Species}{Species of the TF}
  \item{TF_Status}{Motif status of the TF. 'D' stands for directly determined motif.  'I' indicates that the motif is inferred from another TF, based on DBD similarity (see Weirauch et al. 2013 for details). 'N' means no motif is available.}
  \item{Family_Name}{Name of the TF's family}
  \item{DBDs}{The unique set of DBDs (Pfam names) present in the TF}
  \item{DBD_Count}{Number of unique DBDs in the TF}
  \item{Cutoff}{Cutoff used to infer motifs for the TF family}
  \item{DBID}{Motif ID from the associated database or study}
  \item{Motif_Type}{Experimental assay used to determine the motif}
  \item{MSource_Identifier}{ID for the source of the motif (i.e., its project name)}
  \item{MSource_Type}{Internal CisBP ID for the motif category}
  \item{MSource_Author}{First author for the source of the motif}
  \item{MSource_Year}{Year of publication of the motif source}
  \item{PMID}{Pubmed ID of the motif source}
  \item{MSource_Version}{Version of the source (i.e. database build)}
  \item{TFSource_Name}{Source of the TF (i.e. where did the genome build come from?)}
  \item{TFSource_URL}{URL of the TF source}
  \item{TFSource_Year}{Year the genome data was downloaded}
  \item{TFSource_Month}{Month the genome data was downloaded}
  \item{TFSource_Day}{Day the genome data was downloaded}
  \item{\emph{motif_existing}}{Status indicating PWM data is existing or not}
}
\seealso{
See Also as \code{\link{CisBP.group}}, \code{\link{CisBP.extdata}}, \code{\link{CisBP.zipload}}, \code{\link{CisBP.download}}
}
\examples{

# Load the internal CisBP dataset
db_human <- CisBP.extdata("Homo_sapiens");

df.tfinfo <- CisBP.getTFinformation( db_human, tf.information.type = 2)
show(head(df.tfinfo));

}
\keyword{ CisBP object }
