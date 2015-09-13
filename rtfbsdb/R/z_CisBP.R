## Create a genral TFBS database class (tfbs.db) 
## 

setClass("tfbs.db", 
	representation(
		species = "character"
	)
)

#' Create a CisBP database class (extending tfbs.db)
#'  
#' Download Link: http://cisbp.ccbr.utoronto.ca/bulk_archive.php
#' paper: Determination and inference of eukaryotic transcription 
#'        factor sequence specificity, Cell, 2014
#'
#' Thress methods for CisBP.db
#'
#' CisBP.group: get the statistical summary( SQL/sum ) according to a specified field
#
#' CisBP.getTFinformation: get the TF inforation with the status indicating the motif is available.
#'
#' tfbs.createFromCisBP : loading PWM files and create a tfbs object
#'
#' show(command)

setClass("CisBP.db",
	representation(
		zip.url     = "character",  ## http url for data download
		zip.file    = "character",  ## zip file including pwm files and TF_information.txt
		zip.date    = "character",  ## download date    
		file.tfinfo = "character"   ## "TF_Information.txt" 
	), contains = "tfbs.db"
)

setMethod("show", "CisBP.db", function(object){
	cat("Species: ", object@species, "\n");
	cat("Data Source: ", object@zip.url, "\n");
	cat("Zip file: ", object@zip.file, "\n");
	cat("Download date: ", object@zip.date, "\n");
	cat("TF information: ", object@file.tfinfo, "\n");

})

#' download zip file and store it to temporary folder. 
#'
#' URL: http://cisbp.ccbr.utoronto.ca/bulk.php 
#' @param species: the values are available on the web page
#' @param URL: the URL of bulk dowbnloads.
#'
#' @return: CisBP.db object;

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
		zip.date= format( Sys.time(), "%m/%d/%Y"),
		zip.url= zip.url);  
}

#' Construct tfbs.DB object from zipped CisBP file.
#'
#' @param zip.file: the zip file downloaded from the CisBP.
#' @param species: the species for note
#'
#' @return: CisBP.db object;

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
		species  = species, 
		file.tfinfo = paste(tmp.dir, "TF_Information.txt", sep="/"),
		zip.file = zip.file,
		zip.date = format( file.info( zip.file )$mtime, "%m/%d/%Y"),
		zip.url  = ""); 
}

#' find matched zip file and return the create date.
#' 
#'
zipfile.match <- function( species )
{
	file.wildcard <- paste( species, "_ver_*.zip", sep="");
	
	path <- system.file("extdata", package="rtfbsdb");
	files <- list.files( path, glob2rx( file.wildcard ) );
	
	if( length(files)!=1 )
		return( NULL )
	else
	{
		zip.file <- system.file("extdata", files[1], package="rtfbsdb");
		r <- file.info( zip.file );

		zip.base <- basename(zip.file);	
		zip.ver <- regexpr("ver_", zip.base, fixed=T)[1];
		zip.date <- substring(zip.base, zip.ver + 4, zip.ver + 13 )

		#return( list( name=zip.file, date=format( r$mtime, "%m/%d/%Y") ) );
		return( list( name=zip.file, date=zip.date ) );
	}
}

#' Construct tfbs.DB object from inner zip file stored in this package.
#'
#' @param species: for human, Homo_sapiens_2015_04_09.zip is used.
#'                 for mouse, Mus_musculus_2015_04_09.zip
#'
#' @return: CisBP.db object;

CisBP.extdata<-function( species )
{
	zip.file <- "";	
	if ( species=="Homo_sapiens" ||  species=="Human" ||  species=="human") 
	{
		zip.file <- zipfile.match( "Homo_sapiens" );
		species <- "Homo_sapiens";
	}
	else if ( species=="Mus_musculus" || species=="Mouse" || species=="mouse") 
	{
		zip.file <- zipfile.match( "Mus_musculus" );
		species <- "Mus_musculus";
	}
	else if ( species=="Drosophila_melanogaster" || species=="dm3" ) 
	{
		zip.file <- zipfile.match( "Drosophila_melanogaster" );
		species <- "Drosophila_melanogaster";
	}
	else
		stop( paste("No pre-installed zip file for ", species, "."));

	if( is.null(zip.file) )
		stop( paste("No pre-installed zip file for ", species, "."))
	else
		cat("  Zip file for", species, "is downloaded on", zip.file$date, "\n");

	tmp.dir <- tempdir();
	r.file <- unzip( zip.file$name, c("TF_Information.txt"), exdir=tmp.dir );
	if(length(r.file)<1)
	{
		cat("! TF_Information.txt can not be found in the zip file(", zip.file, ")\n");
		return(NULL);
	}

	# dont delete this file, the following functions will use it!
	# unlink(r.file);

	new("CisBP.db", 
		species = species, 
		file.tfinfo = paste(tmp.dir, "TF_Information.txt", sep="/"),
		zip.file= zip.file$name,
		zip.date= zip.file$date,
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

setGeneric("CisBP.group", 
	def=function(cisbp.db, 
				group.by = c("tf_name", "tf_species", "tf_status", "family_name", "motif_type", "msource_id"), 
				tf.information.type = NA) {
		standardGeneric("CisBP.group")
	})

setMethod("CisBP.group", signature(cisbp.db="CisBP.db"),
	function(cisbp.db, group.by=c("tf_name", "tf_species", "tf_status", "family_name", "motif_type", "msource_id"), 
			tf.information.type=1 ) {
	
	group_by <- match.arg(group.by);

	file.tfinfo <- cisbp.db@file.tfinfo;
	if(!missing(tf.information.type))
		file.tfinfo <- CisBP.active_motif_info( cisbp.db, tf.information.type );

	tb.motifs <- read.csv( file.tfinfo, header=T, sep="\t");

	col.idx <- which(toupper(group_by) == toupper(colnames(tb.motifs)) )
	gr <- aggregate( tb.motifs["TF_ID"], by=as.data.frame(tb.motifs[,col.idx]), FUN=length)    
	colnames(gr) <- c( group_by, "count");

	gr;	    
})

#' Get the TF_information file with the existing status for each motif
#'
#' @param cisbp.db: CisBP object
#' @param tf.information.type: 1,2 or 3 indicate which motif file will be used.
#'
#' @return: data.frame;

setGeneric("CisBP.getTFinformation", 
	def=function(cisbp.db, tf.information.type = NA ) {
		standardGeneric("CisBP.getTFinformation")
	})

setMethod("CisBP.getTFinformation", signature(cisbp.db="CisBP.db"),
	function(cisbp.db, tf.information.type = NA )
	{
		file.tfinfo <- cisbp.db@file.tfinfo;

		if(!missing(tf.information.type))
			file.tfinfo <- CisBP.active_motif_info( cisbp.db, tf.information.type );

		tbm <- try( read.csv(file.tfinfo , header=T, sep="\t") );
		if( class(tbm) == "try-error" )
			stop("Can not open TF information file:", file.tfinfo );

		tmp.dir <- tempdir();
		motif.existing <- c();

		for (i in 1:NROW(tbm))
		{
			motif_id <- as.character(tbm$Motif_ID[i]);
		
			if( as.character(motif_id)==".")
			{
				motif.existing <- c(motif.existing, FALSE);
				next;
			}

			pwm.file <- paste( "pwms_all_motifs/", motif_id, ".txt", sep="");
			r.file <- unzip( cisbp.db@zip.file, c(pwm.file), exdir=tmp.dir  );
			if(length(r.file)<1)
			{
				motif.existing <- c(motif.existing, FALSE);
				next;
			}

			tb.motif <- try( read.table(paste(tmp.dir, pwm.file, sep="/"), header=T, sep="\t", row.names=1), TRUE );
			if(class(tb.motif)=="try-error") 
			{
				motif.existing <- c(motif.existing, FALSE);
				next;
			}

			if(NROW(tb.motif)==0)
			{
				motif.existing <- c(motif.existing, FALSE);
				next;
			}

			motif.existing <- c(motif.existing, TRUE);
		}
	
	
		cat("  Total motif:", NROW(tbm), ", Missing or empty PWMs:", length(which(!motif.existing)), ".\n");
	
		return( cbind(tbm, Motif_Existing=motif.existing) );
	})

#' read TF_information file and load PWM files
#'
#' @param cisbp.db: CisBP object
#' @param tf_name
#' @param tf_status
#' @param family_name
#' @param motif_type
#' @param msource_id
#' @param tf.information.type
#'
#' @return: tfbs object

setGeneric("tfbs.createFromCisBP", 
		def=function(cisbp.db, tf_name = NULL, tf_status = NULL, family_name = NULL, 
					motif_type = NULL, msource_id =NULL, tf.information.type = 1){
			stopifnot(class(cisbp.db) == "CisBP.db")
			standardGeneric("tfbs.createFromCisBP")
	})

setMethod("tfbs.createFromCisBP", c( cisbp.db = "CisBP.db" ), tfbs_createFromCisBP )
