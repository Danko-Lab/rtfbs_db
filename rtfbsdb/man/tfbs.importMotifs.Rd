\name{tfbs.importMotifs}
\alias{tfbs.importMotifs}
\title{ Import motifs to tfbs object }
\description{ Import licensed motifs or the source other than Cis-BP to tfbs object }
\usage{
tfbs.importMotifs(tfbs, format, filenames, motif_ids=NULL, PPM.format=TRUE,
      skip.lines=0, pseudocount= -7, force_even= FALSE, ...)
}
\arguments{
  \item{tfbs}{A tfbs object (\code{"\linkS4class{tfbs}"}) returned by \code{\link{tfbs.createFromCisBP}}, \code{\link{tfbs}}, \code{\link{tfbs.dirs}}. }
  \item{format}{String value indicating predefined format or string vector containing some tags to define a customize format. Predefined formats include \strong{pwm.matrix}, \strong{jaspar}, \strong{meme}, \strong{mscan}, \strong{transfac} and \strong{HOCOMOCO}. }
  \item{filenames}{Vector of file names.}
  \item{motif_ids}{Vector of motif ID only for the format \strong{pwm.matrix} . Theses motif IDs are assigned to the motif data loaded from 'pwm.matrix' files. If these motifs are missing or licensed motifs, these IDs have to be in accordance with the TF information which can be exported from 'TF_information.txt' by the function \code{CisBP.getTFinformation}. }
  \item{PPM.format}{Boolean value indicating whether the matrix data is format by PPM(position probability matrix) or PFM(position frequency matrix), otherwise, the original data is PWM format. Default: FALSE. }
  \item{skip.lines}{Number indicating specified non-empty lines are skipped from the head of data file. This parameter is not used for the  format 'pwm.matrix'. Default: 0.}
  \item{pseudocount}{Number value indicating to replace -Inf in PWM. Default: -7. }
  \item{force_even}{Boolean value indicating whether to append a row to make the row number even. Default: FALSE. }
  \item{\dots}{Optional paramaters used in the function \code{read.table} for the 'pwm.matrix' format. }
}

\details{

Two goals to import the motifs from the sources other than Cis-BP. \cr
1) Fill the missing motifs mainly licensed in the Cis-BP database.\cr
2) Make use of the different source.\cr

Five predefined formats are available, \strong{'pwm.matrix'}, \strong{'jaspar'}, \strong{'meme'}, \strong{'mscan'}, \strong{'transfac'} and \strong{'HOCOMOCO'}.\cr

\strong{'pwm.matrix'} \cr
Single text file containing 5 colums, Position, A, C, G and T as shown below.\cr
\preformatted{
Pos	A	C	G	T
1	0.2	0.0	0.4	0.4
2	0.0	0.4	0.4	0.2
3	0.2	0.2	0.2	0.4
4	0.6	0.0	0.4	0.0
5	0.0	0.4	0.0	0.6
6	0.0	0.4	0.4	0.2
7	0.0	0.0	1.0	0.0
8	0.2	0.0	0.0	0.8
9	0.0	0.6	0.0	0.4
10	0.6	0.4	0.0	0.0
11	0.0	1.0	0.0	0.0
}

\strong{'jaspar'}\cr
The package use a predifined template to load the 'jaspar' format as show below.\cr
\preformatted{
> Mycn
A [ 0 29 0 2 0 0 ]
C [31 0 30 1 3 0 ]
G [ 0 0 0 28 0 31]
T [ 0 2 1 0 28 0 ]
}

\strong{'meme'}\cr
The 'meme' format reference: \url{http://meme-suite.org/}\cr
The package use a predifined template to load the 'meme' format as show below.\cr

\preformatted{
MOTIF crp alternative name
	letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009 
 0.000000  0.176471  0.000000  0.823529 
 0.000000  0.058824  0.647059  0.294118 
 0.000000  0.058824  0.000000  0.941176 
 0.176471  0.000000  0.764706  0.058824 
 0.823529  0.058824  0.000000  0.117647 
 0.294118  0.176471  0.176471  0.352941 
 0.294118  0.352941  0.235294  0.117647 
 0.117647  0.235294  0.352941  0.294118 
 0.529412  0.000000  0.176471  0.294118 
 0.058824  0.235294  0.588235  0.117647 
 0.176471  0.235294  0.294118  0.294118 
 0.000000  0.058824  0.117647  0.823529 
 0.058824  0.882353  0.000000  0.058824 
 0.764706  0.000000  0.176471  0.058824 
 0.058824  0.882353  0.000000  0.058824 
 0.823529  0.058824  0.058824  0.058824 
 0.176471  0.411765  0.058824  0.352941 
 0.411765  0.000000  0.000000  0.588235 
 0.352941  0.058824  0.000000  0.588235 
}

\strong{'mscan'}\cr
The 'mscan' format reference: \url{http://www.cisreg.ca/cgi-bin/mscan/MSCAN}\cr
The package use a predifined template to load the 'mscan' format as show below.\cr

\preformatted{
>mef2
10  0  0  0 22  0  6  2  3  4 22 10
0  2 12  0  0  0  0  0  0  0  0  0
9 20  2  0  0  0  0  0  0  0  0 10
3  0  8 22  0 22 16 20 19 18  0  2
>myf
7  9  4  0 16  7  0  6  0  0  6  0
8  0  2 15  0  0 15  0  0 10  0  0
1  7 10  1  0  9  1  0 16  6  0 16
0  0  0  0  0  0  0 10  0  0 10  0
}

\strong{'transfac'}\cr
The transfac format reference: \url{http://www.cisreg.ca/cgi-bin/mscan/MSCAN}\cr
The package use a predifined template to load the 'transfac' format as show below.\cr

\preformatted{
AC  MA0001.1
XX
ID  AGL3
XX
DE  MA0001.1 AGL3; from JASPAR
PO       A     C     G     T
1        0    94     1     2
2        3    75     0    19
3       79     4     3    11
4       40     3     4    50
5       66     1     1    29
6       48     2     0    47
7       65     5     5    22
8       11     2     3    81
9       65     3    28     1
10       0     3    88     6
XX
CC  program: jaspar
XX
//
}

The package implemented a simple parser to load motifs from the different sources. This parser basically reads the data file word by word according to the format rules defined in advance for diffierent source. \cr
The format rules are comprised of some fixed vocabularies and tags defined by the parser or the user. So there are two kinds of tags to describe a motif format, predefined tags and user-defined tags. All tags start with the dollar($) symbol and meet the requirement of program identifier, such as $Express_Pvalue,$D1.\cr
The predefined tags include the following names, which define the motif information and control the parser's cursor.\cr

\tabular{rlll}{
    [1] \tab '>' '[' ']'	\tab control  \tab Start-stop characters for a line or a motif \cr
    [2] \tab $SKIP n		\tab control  \tab Ignore or skip the rest part or next \emph{n} lines.\cr
    [3] \tab $REPEAT    	\tab control  \tab Repeat use the current format until it can't be matched. \cr
    [4] \tab $EOM    		\tab control  \tab End of Motif.\cr
    [5] \tab $LOM   		\tab control  \tab Line of Motif.\cr
    [6] \tab $Motif_ID     	\tab variable \tab Motif_ID required in the 'tfbs' object.\cr
    [7] \tab $TF_Name     	\tab variable \tab TF_Name required in the 'tfbs' object.\cr
    [8] \tab $A+     		\tab variable \tab Multiple \emph{A} nucleobase frequencies in one row\cr
    [9] \tab $C+     		\tab variable \tab Multiple \emph{C} nucleobase frequencies in one row\cr
    [10] \tab $G+     		\tab variable \tab Multiple \emph{G} nucleobase frequencies in one row\cr
    [11] \tab $T+     		\tab variable \tab Multiple \emph{T} nucleobase frequencies in one row\cr
    [12] \tab $A     		\tab variable \tab Single \emph{A} nucleobase frequence in one column\cr
    [13] \tab $C     		\tab variable \tab Single \emph{C} nucleobase frequence in one column\cr
    [14] \tab $G     		\tab variable \tab Single \emph{G} nucleobase frequence in one column\cr
    [15] \tab $T     		\tab variable \tab Single \emph{T} nucleobase frequence in one column\cr
  }

The user-defined tags reprsent the variables which values can be collected to the 'tf_info' slot in the tfbs object (\code{"\linkS4class{tfbs}"}). e.g. $Description, $anyword.

The following section shows 5 predefined format:\cr

\strong{'jaspar'}
\preformatted{
>$Motif_ID $TF_Name
A [ $A+ ]
C [ $C+ ]
G [ $G+ ]
T [ $T+ ]
}

\strong{'meme'}
\preformatted{
MOTIF $Motif_ID $TF_Name $SKIP 
letter-probability matrix: $SKIP
$A  $C  $G  $T $REPEAT
URL $SKIP
}

\strong{'mscan'}
\preformatted{
>$Motif_ID $TF_Name $SKIP
$A+
$C+    
$G+     
$T+
}

\strong{'transfac'}
\preformatted{
AC  $Motif_ID
XX  $SKIP
ID  $TF_Name
XX  $SKIP
DE  $Description
PO       A     C     G     T
$LOM    $A    $C    $G     $T $REPEAT
XX $SKIP
CC $SKIP $REPEAT
// $EOM
}


\strong{'HOCOMOCO'}
\preformatted{
>$Motif_ID,
$A $C $G $T $REPEAT
}


}

\value{
A new tfbs object (\code{"\linkS4class{tfbs}"}) merged with licensed motifs.
}

\seealso{
\code{\link{tfbs.createFromCisBP}}
}
\examples{

require(rtfbsdb);

db  <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, family_name="AP-2");
tfs;


#It is test code
#motif_ids <- c( "M2938_1.02", "M2940_1.02", "M2940_2.02", "M4056_1.02" );
#path <- system.file("extdata", package="rtfbsdb");
#file_pwms <- paste(path, c( 
#      "fake_M2938_1.02.pwm", 
#      "fake_M2940_1.02.pwm",
#      "M3591_1.01.pwm", 
#      "M3590_1.01.pwm",
#      ), sep="/");
#
#
#tfs <- tfbs.importMotifs(tfs, 'pwm.matrix', file_pwms, motif_ids, header=T );
#show(tfs);


## import 2 motifs to fill the licensed motifs in Cis-BP and 1 new motif fromother source
motif_ids <- c( "M2938_1.02", "M3591_1.01", "M3590_1.01" );
path <- system.file("extdata", package="rtfbsdb");
file_pwms <- paste(path, c( 
      "fake_M2938_1.02.pwm",
      "M3591_1.01.pwm", 
      "M3590_1.01.pwm"), sep="/");

tfs <- tfbs.importMotifs(tfs, 'pwm.matrix', file_pwms, motif_ids, header=TRUE );
show(tfs);


## Data is copied from http://rsat01.biologie.ens.fr/rsa-tools/help.convert-matrix.html
##
data.transfac <- system.file("extdata", "pwm.example.transfac.txt", package="rtfbsdb");
tfs.transfac <- tfbs.importMotifs(tfs, "transfac", data.transfac, skip.lines=0);
show(tfs.transfac);

## Data from 
## http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt
##
data.jaspar <- system.file("extdata", "pwm.example.jaspar.2015.txt", package="rtfbsdb");
tfs.jaspar <- tfbs.importMotifs(tfs, "jaspar", data.jaspar, skip.lines=0);
show(tfs.jaspar);

## Data from 
## http://jaspar.genereg.net/html/DOWNLOAD/ARCHIVE/JASPAR2010/JASPAR_CORE/pfms/pfms_all.txt
##
data.mscan <- system.file("extdata", "pwm.example.mscan.jaspar2010.txt", package="rtfbsdb");
tbs <- tfbs();
tfs.mscan <- tfbs.importMotifs(tbs, "mscan", data.mscan, skip.lines=0);
show(tfs.mscan);

## Data from http://meme-suite.org/doc/download.html?man_type=web
##
data.meme <- system.file("extdata", "pwm.example.meme.Homo_sapiens.txt", package="rtfbsdb");
tbs <- tfbs();
tfs.meme <- tfbs.importMotifs(tbs, "meme", data.meme, skip.lines=5);
show(tfs.meme);

## Data from http://www.nature.com/nature/journal/v527/n7578/full/nature15518.html
##
format.style <- c("Base $Motif_ID $TF_Name $Experiment 
$Ligand_sequbatch $Seed $Multinomial $Cycle $Is_matrix $Comment 
$Genomic_pvalue $Enrichment_pvalue $Category $SKIP",
		"A $A+",
		"C $C+",
		"G $G+",
		"T $T+" );
data.file <- system.file("extdata", "pwm.example.nature15518.s1.txt", package="rtfbsdb");

tbs <- tfbs();
tfs.nature15518 <- tfbs.importMotifs(tbs, format.style, data.file, skip.lines=19)
show(tfs.nature15518);


## HOCOMOCO_10
## Data from http://hocomoco.autosome.ru/downloads
## The matrix is PWM format.

library(rtfbsdb)
format.style <- c(">$Motif_ID",
		"$A $C $G $T $REPEAT");
data.file <- system.file("extdata", "HOCOMOCO_10_pwms_HUMAN_mono.txt", package="rtfbsdb");

tbs <- tfbs();
tfs.HOCOMOCO.human <- tfbs.importMotifs(tbs, format.style, data.file, skip.lines=0, PPM.format=FALSE)
show(tfs.HOCOMOCO.human);

}
\references{
1.\url{http://rsat01.biologie.ens.fr/rsa-tools/help.convert-matrix.html}\cr
2.\url{http://www.cisreg.ca/cgi-bin/mscan/MSCAN}\cr
3.\url{http://meme-suite.org/}\cr
4.\url{http://hocomoco.autosome.ru/downloads}\cr
}
\keyword{ tfbs object }
