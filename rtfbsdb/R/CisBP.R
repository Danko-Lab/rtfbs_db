## Create a genral TFBS database class (tfbs.db) 
## 
## Two methods fortfbs.db
## tfbs.group: get the statistical summary( SQL/sum ) 
##             according to a field
## CisBP.create : get the subset by the specified parameter.
##

setClass("tfbs.db", 
  representation(
    species = "character"   ## 
    )
  )

setGeneric("tfbs.createFromCisBP", 
    def=function(cisbp.db, tf_name = NULL, tf_status = NULL,
                 family_name = NULL, motif_type = NULL, msource_id =NULL, tf.information.type = 1,
                 expressed.only=TRUE, include.DBID.Missing=TRUE, seq.datatype=NA,
    		 	 file.bigwig.plus=NA, file.bigwig.minus=NA, file.bam=NA, file.twoBit=NA, file.gencode.gtf=NA, ncores = 1) {
	  standardGeneric("tfbs.createFromCisBP")
	})

setGeneric("CisBP.group", 
    def=function(cisbp.db, group.by = c("tf_name", "tf_species",
                 "tf_status", "family_name", "motif_type",
                 "msource_id"), tf.information.type = 1) {
	  standardGeneric("CisBP.group")
	})

## Create a CisBP database class (extending tfbs.db)
##  
## Download Link: http://cisbp.ccbr.utoronto.ca/bulk_archive.php
## paper: Determination and inference of eukaryotic transcription 
##        factor sequence specificity, Cell, 2014

setClass("CisBP.db", 
  representation(
    zip.url     = "character",  ## http url for data download
    zip.file    = "character",  ## zip file including pwm files and TF_information.txt
    file.tfinfo = "character"   ## "TF_Information.txt" 
    ), contains = "tfbs.db"
  )

#' download zip file and store it to temporary folder. 
#'
#' URL: http://cisbp.ccbr.utoronto.ca/bulk.php 
#' @param species: the values are available on the web page
#' @param URL: the URL of bulk dowbnloads.
#'
#' @return: tfbs.db object;

CisBP.download <- function( species="Homo_sapiens", url="http://cisbp.ccbr.utoronto.ca/bulk_archive.php" ) 
{
  if(! (require(RCurl) && require (stringr) ) )
    stop("Package RCurl and stringr are necessary to download files.");
  
  #test the HTTP header 
  #if(is.null(url)) url <- "http://127.0.0.1:8888/";
  if( is.null(url) ) url <- "http://cisbp.ccbr.utoronto.ca/bulk_archive.php";
 
  if(!RCurl::url.exists(url))
     stop( paste("URL is invalid(", url, ")." ));

  cat("  CisBP Download\n");
  cat("* Species =", species, "\n");
  cat("* URL =", url, "\n");
  
  f.encode <- function(x)
  {
    if(x!="Download+Species+Archive%21")
     	return( RCurl::curlEscape(x) )
    else	
		return(x);
  }

  # to append specific headers in HTTP request(Accept).

  h.curl = RCurl::getCurlHandle();
  hidden <- RCurl::curlSetOpt( .opts = list(httpheader = c(Accept="text/html,application/xhtml+xml,application/xml"), verbose = TRUE), curl=h.curl)

  cat("  Posting the form data to the web page and get the link of zip file ...\n");
  
  # HTTP post menthod is called to send the form parameters(selSpec, ...)
  t <- RCurl::postForm( url, curl=h.curl, 
	selSpec   = species,
	"Spec%5B%5D"  = "Logos",
	"Spec%5B%5D"  = "TF_Information",
	"Spec%5B%5D"  = "PWMs",
	submit    = "Download+Species+Archive%21",
	style     = "POST", 
	.contentEncodeFun =f.encode);
	
  

  # parser the return page and get the link of zip file
  zip.path <- stringr::str_extract(as.character(t), "tmp/[^ \f\n\r\t\v]+.zip\">Download")
  zip.path <- substring( zip.path, 1, nchar(zip.path)-10 );
	
  zip.url <- paste("http://cisbp.ccbr.utoronto.ca/", zip.path, sep="");
  cat("* Zip file =", zip.url, "\n");
  
  cat("  Downloading the zip file ...\n");
  
  zip.file = tempfile();
  if( download.file( zip.url, zip.file ) != 0 )
    stop( paste("Failed to download file from url=", zip.url, sep="") );

  cat("  Detecting TF_Information.txt ...\n");
  tmp.dir <- tempdir();
  r.file <- unzip( zip.file, c("TF_Information.txt"), exdir=tmp.dir );
  if(length(r.file)<1)
  {
    cat("! TF_Information.txt can not be found in the zip file(", zip.file, ")\n");
  	return(NULL);
  }

  #unlink(r.file);
  
  cat("* TF_Information =", r.file, "\n");

  new("CisBP.db", 
	species = species, 
	file.tfinfo = paste(tmp.dir, "TF_Information.txt", sep="/"),
	zip.file= zip.file,
	zip.url= zip.url);  
}

#' Construct tfbs.DB object from zipped CisBP file.
#'
#' @param zip.file: the zip file downloaded from the CisBP.
#' @param species: the species for note
#'
#' @return: tfbs.db object;

CisBP.zipload <- function( zip.file, species="Homo_sapiens" ) 
{
  tmp.dir <- tempdir();
  r.file <- unzip( zip.file, c("TF_Information.txt"), exdir=tmp.dir );
  if(length(r.file)<1)
  {
    cat("! TF_Information.txt can not be found in the zip file(", zip.file, ")\n");
  	return(NULL);
  }

  #unlink(r.file);

  new("CisBP.db", 
	species = species, 
	file.tfinfo = paste(tmp.dir, "TF_Information.txt", sep="/"),
	zip.file= zip.file,
	zip.url= ""); 
}

#' Construct tfbs.DB object from inner zip file stored in this package.
#'
#' @param species: for human, Homo_sapiens_2015_04_09.zip is used.
#'                 for mouse, Mus_musculus_2015_04_09.zip
#'
#' @return: tfbs.db object;

CisBP.extdata<-function( species )
{
  zip.file <- "";	
  if ( species=="Homo_sapiens" ||  species=="Human" ||  species=="human") 
  {
  	zip.file <- system.file("extdata","Homo_sapiens_2015_04_09.zip", package="rtfbsdb");
  	species <- "Homo_sapiens";
  }
  else if ( species=="Mus_musculus" || species=="Mouse" || species=="mouse") 
  {
  	zip.file <- system.file("extdata","Mus_musculus_2015_04_09.zip", package="rtfbsdb");
  	species <- "Mus_musculus";
  }
  else
    stop( paste("No zip file for ", species, "."));
    
  tmp.dir <- tempdir();
  r.file <- unzip( zip.file, c("TF_Information.txt"), exdir=tmp.dir );
  if(length(r.file)<1)
  {
    cat("! TF_Information.txt can not be found in the zip file(", zip.file, ")\n");
  	return(NULL);
  }
  
  unlink(r.file);
  
  new("CisBP.db", 
	species = species, 
	file.tfinfo = paste(tmp.dir, "TF_Information.txt", sep="/"),
	zip.file= zip.file,
	zip.url= "extdata");  
}

#' Select the motif table for CisBP.create and CisBP.group function. 
#'
#' Three TF information files in zip file
#'
#' 1: TF_Information.txt : (direct motifs) or (no direct but inferred motifs with 90%)
#' 2: TF_Information_all_motifs.txt: (direct motifs) and (inferred motifs above the threshold)
#' 3: F_Information_all_motifs_plus.txt: All motifs
#' 
#' @param cisbp.db: CisBP.db object
#' @param tf.information.type: 1,2 or 3 indicate the index of above files.
#'
#' @return: temporary file;

CisBP.active_motif_info<-function( cisbp.db, tf.information.type=1 )
{
    motif_infos <- c( "TF_Information.txt", "TF_Information_all_motifs.txt", "TF_Information_all_motifs_plus.txt");
    
    if(tf.information.type<1 || tf.information.type>3)
    {
       cat("! The meta data can be stored in TF_Information.txt(1), TF_Information_all_motifs.txt(2) or TF_Information_all_motifs_plus.txt(3).");
       return(cisbp.db@file.tfinfo);
    }
    
    file.tfinfo <- motif_infos[ tf.information.type ];
    
    tmp.dir <- tempdir();
    r.file <- unzip(cisbp.db@zip.file, c(file.tfinfo), exdir=tmp.dir );
    if(length(r.file)<1)
    {
    
      cat("! The meta file (", file.tfinfo ,") can not be found in the zip file.");
      return( cisbp.db@file.tfinfo );
    }
    
    return( paste(tmp.dir, file.tfinfo, sep="/") );
}

#' Get the statistical summary by grouping the fields in the motif table
#'
#' @param cisbp.db: CisBP object
#' @param group.by: available values are tf_name, tf_species, tf_status, family_name, motif_type and msource_id.
#' @param tf.information.type: 1,2 or 3 indicate which motif file will be used.
#'
#' @return: data.frame;

setMethod("CisBP.group", signature(cisbp.db="CisBP.db"),
    function(cisbp.db, group.by=c("tf_name", "tf_species", "tf_status", "family_name", "motif_type", "msource_id"), 
    		tf.information.type=1 )
    {
      group_by <- match.arg(group.by);
      
      cisbp.db@file.tfinfo <- CisBP.active_motif_info( cisbp.db, tf.information.type );
      tb.motifs <- read.csv(cisbp.db@file.tfinfo, header=T, sep="\t");
      
      col.idx <- which(toupper(group_by) == toupper(colnames(tb.motifs)) )
      gr <- aggregate( tb.motifs["TF_ID"], by=as.data.frame(tb.motifs[,col.idx]), FUN=length)    
      gr;	    
})


#' Find the subset by querying the motif table
#'
#' @param cisbp.db: cisbp.db object
#' @param tf_name: string, the query value for tf_name 
#' @param tf_status: string, the query value for tf_status 
#' @param family_name: string, the query value for family_name  
#' @param motif_type: string, the query value for motif_type  
#' @param msource_id: string, the query value for msource_id  
#' @param tf.information.type: 1,2 or 3 indicate which motif file will be used.
#'
#' @return: NULL or tfbs object;

setMethod("tfbs.createFromCisBP", signature(cisbp.db="CisBP.db"),
    function(cisbp.db, tf_name=NULL, tf_status=NULL, family_name=NULL, motif_type=NULL, msource_id=NULL , tf.information.type=1, 
    		 expressed.only=TRUE, include.DBID.Missing=TRUE, seq.datatype=NA, 
    		 file.bigwig.plus=NA, file.bigwig.minus=NA, file.bam=NA, file.twoBit=NA, file.gencode.gtf=NA, ncores = 1) 
{
    if( missing(seq.datatype)) seq.datatype <- "GRO-seq";
    if( !(seq.datatype %in% c("GRO-seq", "PRO-seq", "RNA-seq") ) )
    	stop("Only GRO-seq, PRO-seq and RNA-seq data can use this package.");

    cisbp.db@file.tfinfo <- CisBP.active_motif_info( cisbp.db, tf.information.type );

    tbm <- try( read.csv(cisbp.db@file.tfinfo, header=T, sep="\t") );
    if( class(tbm) == "try-error" )
		stop("! Can not open TF information file:", cisbp.db@file.tfinfo, ".\n" );
    
    tbm_f <- c();

    if(!is.null(tf_name))    tbm_f <- c( tbm_f, paste( "tbm$TF_Name=='", tf_name, "'", sep="") );
    if(!is.null(tf_status))  tbm_f <- c( tbm_f, paste( "tbm$TF_Status=='", tf_status, "'", sep="") );
    if(!is.null(family_name))tbm_f <- c( tbm_f, paste( "tbm$Family_Name=='", family_name, "'", sep="") );
    if(!is.null(motif_type)) tbm_f <- c( tbm_f, paste( "tbm$Motif_Type=='", motif_type, "'", sep="") );
    if(!is.null(msource_id)) tbm_f <- c( tbm_f, paste( "tbm$MSource_Identifier=='", msource_id, "'", sep="") );
    
    if(length(tbm_f)>0)
    {
      tbm_f_all <- paste(tbm_f, collapse=" & " );
      nidx <-  eval(parse(text=paste("which(", tbm_f_all, ")")));
      if(length(nidx)<1) return( NULL );
    }
    else
    {
      nidx <- c(1:NROW(tbm));
      tbm_f_all <- "All"; 
    }
    
    cat( " ", length(nidx), "PWM files are defined in the CisBP dataset.\n");

    tmp.dir <- tempdir();
    pwm.files <- c();
	nidx.motif <- c();
	names <- c();
	
	err.noPWM <- 0;
	err.noACGT <- 0;
	
    for (i in nidx)
	{
		motif_id <- as.character(tbm$Motif_ID[i]);
		
		if( as.character(motif_id)==".")
		{
			# cat("! No ID for this motif(.).\n"); 	
			err.noPWM <- err.noPWM + 1;
			next;
        }
        
		pwm.file <- paste( "pwms_all_motifs/", motif_id, ".txt", sep="");
		r.file <- unzip( cisbp.db@zip.file, c(pwm.file), exdir=tmp.dir  );
		if(length(r.file)<1)
		{
			# cat("! Can not find PWM file for motif ID=", motif_id, ".\n" );
			err.noPWM <- err.noPWM + 1;
			next;
		}
        
		tb <- try( read.table(paste(tmp.dir, pwm.file, sep="/"), header=T, sep="\t", row.names=1), TRUE );
		if(class(tb)=="try-error") 
		{
			# cat("! Can not find PWM file for motif ID=", motif_id, ".\n" );
			err.noPWM <- err.noPWM + 1;
			next;
		}
	
		if(NROW(tb)==0)
		{
			# cat("! No A C G T values in the PWM file for motif ID=", motif_id, ".\n" );
			err.noACGT <- err.noACGT + 1
			next;
    	}
    	
    	pwm.files <- c( pwm.files, paste(tmp.dir, pwm.file, sep="/") );
    	nidx.motif<- c( nidx.motif, i );
    	names     <- c( names, motif_id);
    }

	TF_info <- NULL;
	if(length(nidx.motif)>0)
		TF_info <- tbm[ nidx.motif, , drop=F];

	if( err.noPWM + err.noACGT > 0 )
		cat("!", err.noPWM + err.noACGT, "PWM file(s) are failed to be loaded ( File missing:", err.noPWM, ", No A C G T values:", err.noACGT, ").\n");

    if(length(pwm.files)==0)
    {
      cat("! No PWM files to create a tfbs object.\n");
      return(NULL);
    }
    else
      cat("*", length(pwm.files), "PWMs in the tfbs object.\n");

    tfs <- tfbs(
      	 species    = cisbp.db@species,
         filenames  = pwm.files, 
         names      = names, 
         extra_info = TF_info,
         header=T, sep="\t" , row.names=1 );

	if( !missing(file.bigwig.plus) && !missing(file.bigwig.minus) )
	{
		tfs <- tfbs.getExpression(tfs, file.bigwig.plus, file.bigwig.minus, file.bam=file.bam, file.twoBit=file.twoBit, file.gencode.gtf=file.gencode.gtf, seq.datatype=seq.datatype, ncores = ncores  );
		if(!is.null( tfs@expressionlevel ))
		{
			if( expressed.only ) tfs <- tfbs.selectExpressed( tfs, 0.05, include.DBID.Missing );
		}
		else
			cat("! Failed to calculate the gene expression and select expressed TFs.\n");
	}
	
    return(tfs);
})


setMethod("show", "CisBP.db", function(object){

	cat("Species: ", object@species, "\n");
	cat("Data Source: ", object@zip.url, "\n");
	cat("Zip file: ", object@zip.file, "\n");
	cat("TF information: ", object@file.tfinfo, "\n");

});
