\name{tfbs.createFromCisBP}
\alias{tfbs.createFromCisBP}
\title{
Create TF object by querying the CisBP dataset.
}
\description{
Build a tfbs object by querying the meta file of CisBP dataset and subsetting the results.
}
\usage{
tfbs.createFromCisBP(cisbp.db, 
      tf_name = NULL, 
      tf_status = NULL, 
      family_name = NULL, 
      motif_type = NULL, 
      msource_id = NULL, 
      tf.information.type = 1)       
}
\arguments{
  \item{cisbp.db}{A CisBP object(\code{"\linkS4class{CisBP.db}"}), including the file of TF_Information.txt.}
  \item{tf_name}{String, indicating the TF_name field will be used to select motifs.}
  \item{tf_status}{String, indicating the TF_Status field will be used to select motifs. }
  \item{family_name}{String, indicating the Family_Name field will be used to select motifs. }
  \item{motif_type}{String, indicating the Motif_Type field will be used to select motifs.  }
  \item{msource_id}{String, indicating the MSource_Identifier field will be used to select motifs. }
  \item{tf.information.type}{Number,indicating which TF meta file will be used. Available values are 1 \cr
  for TF_Information.txt, 2 for TF_Information_all_motifs.txt and 3 for \cr
  TF_Information_all_motifs_plus.txt.}
}
\details{
The function includes three steps to build a tfbs object:\cr
1) Searching the TF information and PWM files in the CisBP dataset according to the criteria specified by 
the parameters of \emph{tf_name}, \emph{tf_status}, \emph{family_name}, \emph{motif_type} and \emph{msource_id}.\cr
2) If \code{seq.datatype} is GRO-seq or PRO-seq and the bigwig files are provided, the gene expression values are calculated 
through querying the TREs region from the GENCODE database( for human, gencode.v21.annotation.gtf, for mouse: gencode.vM3.annotation.gtf)\cr
and querying the reads count in the plus and minus bigWig files.\cr

If \code{seq.datatype} is RNA-seq and the BAM file is provided, read counts for each TRE regions will be queried from the BAM file.\cr

3) If the expressed TFs only is used in the tfbs object, the TFs with p-values corrected by Bonfrroni less than 0.05 will be selected.\cr

The following part explains how to calculate the gene expression.\cr

For each motif, the occurance ranges can be queried by the gene ID  
After the searching, one range obtained from the merge of the multiple ranges will be 
used to detect the reads count in the specified  bigwig files(including plus and minus). 
The probability of each motif can be calcuated by the reads count and lambda.\cr

The lambda is determined by the following formulation:\cr

\code{r.lambda = 0.04 * sum(reads_in_all_chromosomes)/10751533/1000.}\cr

The dataset of GENECODE v21 (human) and vM3 (mouse) have been compiled into RDATA file 
and attached in this package.\cr

The \code{gencode_transcript_ext} object can be accessed after the following command is 
executed successfully.\cr

\code{load( system.file("extdata", "gencode_v21_transcript_ext.rdata",
      package="rtfbsdb"), environment() );}\cr

}
\value{
A tfbs object is returned with PWM matrices, see Also as \code{"\linkS4class{tfbs}"}
}
\seealso{
See Also as \code{\link{tfbs}}
}
\examples{

# Load the internal CisBP dataset
db_human <- CisBP.extdata("Homo_sapiens");

# Load all motifs and return a tfbs object.
tfs0 <- tfbs.createFromCisBP(db_human); 

# Query the motifs by the conditions and return a tfbs object 
tfs1 <- tfbs.createFromCisBP(db_human, family_name="Homeodomain", tf_status="D",  
      motif_type="ChIP-seq", msource_id= "MS01_1.01", tf.information.type=1 ); 

# Query the motifs by the conditions and return a tfbs object 
tfs2 <- tfbs.createFromCisBP(db_human, family_name="Homeodomain", tf_status="D" ); 

# Query the motifs by the conditions and return a tfbs object 
tfs3 <- tfbs.createFromCisBP(db_human, motif_type="ChIP-seq"  ); 

# Query the motifs by the conditions and return a tfbs object 
tfs4 <- tfbs.createFromCisBP(db_human, tf.information.type=2); 

}
\keyword{CisBP object}
\keyword{tfbs object}