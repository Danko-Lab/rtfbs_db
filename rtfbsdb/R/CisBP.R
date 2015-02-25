setClass("tfbs.db", 
  representation(
    species = "character"   ## 
    )
  )

setGeneric("tfbs.find", 
    def=function(tfbs.db, ...) {
	  #stopifnot(class(tfbs.db) == "tfbs.db")
	  standardGeneric("tfbs.find")
	})

setGeneric("tfbs.group", 
    def=function(tfbs.db, group_by, ...) {
	  #stopifnot(class(tfbs.db) == "tfbs.db")
	  standardGeneric("tfbs.group")
	})

setClass("CisBP.db", 
  representation(
    zip.url     = "character",  ## http url for data download
    zip.file= "character",      ## zip file including pwm files and TF_information.txt
    file.tfinfo = "character"   ## "TF_Information.txt" 
    ), contains = "tfbs.db"
  )

## download zip file from http://cisbp.ccbr.utoronto.ca/bulk.php
## example: 
#
#
#
CisBP.download <- function( species="Homo_sapiens", url= NULL, ...) 
{
  if(!require(RCurl))
    stop("Package CUrl is necessary to download files.");
  
  #test the HTTP header 
  #if(is.null(url)) url <- "http://127.0.0.1:8888/";
  if( is.null(url) ) url <- "http://cisbp.ccbr.utoronto.ca/bulk_archive.php";
 
  if(!url.exists(url))
     stop( paste("URL is invalid(", url, ")." ));

  cat("  CisBP Download\n");
  cat("* Species =", species, "\n");
  cat("* URL =", url, "\n");
  
  f.encode <- function(x)
  {
    if(x!="Download+Species+Archive%21")
     	return(curlEscape(x))
    else	
		return(x);
  }

  h.curl = getCurlHandle();
  hidden <- curlSetOpt( .opts = list(httpheader = c(Accept="text/html,application/xhtml+xml,application/xml"), verbose = TRUE), curl=h.curl)
  
  t <- postForm( url, curl=h.curl, 
	selSpec   = species,
	"Spec%5B%5D"  = "Logos",
	"Spec%5B%5D"  = "TF_Information",
	"Spec%5B%5D"  = "PWMs",
	submit    = "Download+Species+Archive%21",
	style     = "POST", 
	.contentEncodeFun =f.encode);
	
  library(stringr)

  #get the link of zip file
  zip.path <- str_extract(as.character(t), "tmp/[^ \f\n\r\t\v]+.zip\">Download")
  zip.path <- substring( zip.path, 1, nchar(zip.path)-10 );
	
  zip.url <- paste("http://cisbp.ccbr.utoronto.ca/", zip.path, sep="");
  cat("* Zip file =", zip.url, "\n");
  
  zip.file = tempfile();
  if( download.file( zip.url, zip.file ) != 0 )
    stop( paste("Failed to download file from url=", zip.url, sep="") );

  r.file <- unzip( zip.file, c("TF_Information.txt") );
  if(length(r.file)<1)
  {
    cat("! TF_Information.txt can not be found in the zip file(", zip.file, ")\n");
  	return(NULL);
  }
  
  cat("* TF_Information =", r.file, "\n");

  new("CisBP.db", 
	species = species, 
	family = family,
	file.tfinfo = "TF_Information.txt",
	zip.file= zip.file,
	zip.url= zip.url);  
}

#
#
#
CisBP.zipload <- function( zip.file, species="Homo_sapiens" ) 
{
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
	file.tfinfo = "TF_Information.txt",
	zip.file= zip.file,
	zip.url= ""); 
}

#
#
#
CisBP.extdata<-function( species="Homo_sapiens" )
{
  zip.file <- "";	
  if ( species=="Homo_sapiens") 
    zip.file <- system.file("extdata","Homo_sapiens_2015_02_02_12:09_pm.zip", package="rtfbsdb");

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

#
#
#
CisBP.active_motif_info<-function( tfbs.db, motif_info_type=1 )
{
    motif_infos <- c( "TF_Information.txt", "TF_Information_all_motifs.txt", "TF_Information_all_motifs_plus.txt");
    
    if(motif_info_type<1 || motif_info_type>3)
    {
       cat("! The meta data can be stored in TF_Information.txt(1), TF_Information_all_motifs.txt(2) or TF_Information_all_motifs_plus.txt(3).");
       return(tfbs.db@file.tfinfo);
    }
    
    file.tfinfo <- motif_infos[ motif_info_type ];
    
    tmp.dir <- tempdir();
    r.file <- unzip(tfbs.db@zip.file, c(file.tfinfo), exdir=tmp.dir );
    if(length(r.file)<1)
    {
    
      cat("! The meta file (", file.tfinfo ,") can not be found in the zip file.");
      return( tfbs.db@file.tfinfo );
    }
    
    return( paste(tmp.dir, file.tfinfo, sep="/") );
}

#
#
#
#
setMethod("tfbs.group", signature(tfbs.db="CisBP.db"),
    function(tfbs.db, group_by=c("tf_name", "tf_species", "tf_status", "family_name", "motif_type", "msource_id"), motif_info_type=1 )
    {
      group_by <- match.arg(group_by);
      
      tfbs.db@file.tfinfo <- CisBP.active_motif_info( tfbs.db, motif_info_type );
      tb.motifs <- read.csv(tfbs.db@file.tfinfo, header=T, sep="\t");
      
      col.idx <- which(toupper(group_by) == toupper(colnames(tb.motifs)) )
      gr <- aggregate( tb.motifs["TF_ID"], by=as.data.frame(tb.motifs[,col.idx]), FUN=length)    
      gr;	    
})


#
#
#
#
setMethod("tfbs.find", signature(tfbs.db="CisBP.db"),
    function(tfbs.db, tf_name=NULL, tf_status=NULL, family_name=NULL, motif_type=NULL, msource_id=NULL , motif_info_type=1) 
{

    tfbs.db@file.tfinfo <- CisBP.active_motif_info( tfbs.db, motif_info_type );

    tbm <- read.csv(tfbs.db@file.tfinfo, header=T, sep="\t");
    
    tbm_f <- c();

    if(!is.null(tf_name))    tbm_f <- c( tbm_f, paste( "tbm$TF_Name=='", tf_name, "'", sep="") );
    if(!is.null(tf_status))  tbm_f <- c( tbm_f, paste( "tbm$TF_Status=='", tf_status, "'", sep="") );
    if(!is.null(family_name))tbm_f <- c( tbm_f, paste( "tbm$Family_Name=='", family_name, "'", sep="") );
    if(!is.null(motif_type)) tbm_f <- c( tbm_f, paste( "tbm$Motif_Type=='", motif_type, "'", sep="") );
    if(!is.null(msource_id)) tbm_f <- c( tbm_f, paste( "tbm$MSource_Identifier=='", msource_id, "'", sep="") );

    if(length(tbm_f)>0)
    {
      tbm_f_all <- paste(tbm_f, collapse=" & " );
      #nidx <- which( x$TF_Name == family_name);
      nidx <-  eval(parse(text=paste("which(", tbm_f_all, ")")));
      if(length(nidx)<1) return( NULL );
    }
    else
    {
      nidx <- c(1:NROW(tbm));
      tbm_f_all <- "All"; 
    }
    
    cat( length(nidx), "PWM files are selected!\n");

    tmp.dir <- tempdir();
    pwm.files <- c();

	nidx.motif <- c();

    for (i in nidx)
	{
		if( as.character(tbm$Motif_ID[i])==".") next;
        
		pwm.file <- paste( "pwms_all_motifs/", tbm$Motif_ID[i], ".txt", sep="");
		r.file <- unzip( tfbs.db@zip.file, c(pwm.file), exdir=tmp.dir  );
		if(length(r.file)<1)
		{
			cat("Can not find PWM file:", tbm$Motif_ID[i], ".txt","\n" );
			next;
		}
        
		tb <- try( read.table(paste(tmp.dir, pwm.file, sep="/"), header=T, sep="\t", row.names=1), TRUE );
		if(class(tb)=="try-error") 
		{
			cat("Can not find PWM file:", tbm$Motif_ID[i], ".txt","\n" );
			next;
		}
	
		if(NROW(tb)==0) next;
    
    	pwm.files <- c( pwm.files, paste(tmp.dir, pwm.file, sep="/") );
    	nidx.motif <- c( nidx.motif, i );
    }

	TF_info <- NULL;
	if(length(nidx.motif)>0)
		TF_info <- tbm[ nidx.motif, , drop=F];

    if(length(pwm.files)>0)
      tfbs(
         filenames  = pwm.files, 
         names      = tbm_f_all, 
         extra_info = TF_info,
         header=T, sep="\t" , row.names=1 )
    else
      return(NULL);
})

