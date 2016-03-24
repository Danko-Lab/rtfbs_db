#' Creates a new tfbs object.  Reads files

tfbs <- function(filenames=NULL, names=NULL, species="Homo_sapiens", tf_info=NULL, tf_missing=NULL, ...) 
{
	stopifnot( length(names) == length(filenames) );
	
	pwms <- list();
	idx.pwm <- c();
	mgisymbols <- c();
	if( !is.null(filenames) )
	{
		for(i in 1:length(filenames)) {
			pwms[[i]] <- read.motif(filenames[i], ...)
		}
	
		idx.pwm <- which( unlist( lapply( pwms, function(x) !is.null(x) ) ) );
		if( length(idx.pwm)==0 )
			return( NULL );

		pwm        = pwms[ idx.pwm ];
		filename   = filenames[ idx.pwm ];
		mgisymbols = names[ idx.pwm ];
	}
	else
		filenames <- c();
	
	if( is.null(tf_info) )
		tf_info  <- as.data.frame(NULL)
	else
		tf_info  <- tf_info[idx.pwm,,drop=F];
	
	if( is.null(tf_missing) )
		tf_missing  <- as.data.frame(NULL)

	new("tfbs", 
		species        = species,
		ntfs           = as.integer(length(filenames)),
		pwm            = pwms,
		filename       = as.character(filenames), 
		mgisymbols     = as.character(mgisymbols), 
		tf_info        = tf_info,
		tf_missing     = tf_missing,
		distancematrix = matrix(, nrow=0, ncol=0),
		cluster        = matrix(, nrow=0, ncol=0),
		expressionlevel= as.data.frame(NULL) )
}

#' Create a new tfbs object from all the PWMs found in the supplied folders.
#' Optionally recursively descends into subfolders.

tfbs.dirs <- function(..., species="Homo_sapiens", args.read.motif = NULL, pattern = glob2rx("*.pwm"), recursive = FALSE) 
{
	pwms = list()
	paths = list(...)
	filenames = NULL
	pwm.names = NULL

	k = 1
	for (path in paths) {
		fnames = list.files(path, pattern = pattern, full.names = TRUE, recursive = recursive)

		for (fname in fnames) {
			pwm = do.call("read.motif", c(list(fname), args.read.motif))
			pwms[[k]] = pwm
			k = k + 1
		}

		filenames = c(filenames, fnames)
		pwm.names = c(pwm.names, sapply(fnames, function(str) {
			tmp = basename(str)
			substr(tmp, 1, nchar(tmp) - 4) # drop .pwm from name
		}))
	}

	idx.pwm <- which( unlist( lapply( pwms, function(x) !is.null(x) ) ) );
	if( length(idx.pwm)==0 )
		return( NULL );

	names(pwm.names) <- NULL # clear filenames

	# build object instance
	new("tfbs", 
		#usemotifs      = as.integer(1:length(filenames)),
		species         = species,
		ntfs            = as.integer(length(idx.pwm)),
		filename        = filenames[idx.pwm], 
		mgisymbols      = pwm.names[idx.pwm], 
		pwm             = pwms[idx.pwm],
		distancematrix  = matrix(, nrow=0, ncol=0),
		cluster         = matrix(, nrow=0, ncol=0),
		expressionlevel = as.data.frame(NULL) );
}

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

tfbs_createFromCisBP <- function ( cisbp.db, 
					motif_id    = NULL, 
					tf_name     = NULL, 
					tf_status   = NULL,
					family_name = NULL, 
					motif_type  = NULL, 
					msource_id  = NULL, 
					tf.information.type = 1 )
{
	cisbp.db@file.tfinfo <- CisBP.active_motif_info( cisbp.db, tf.information.type );

	tbm <- try( read.csv(cisbp.db@file.tfinfo, header=T, sep="\t") );
	if( class(tbm) == "try-error" )
		stop("Can not open TF information file:", cisbp.db@file.tfinfo, ".\n" );

	tbm_f <- c();

	if(!is.null(motif_id))   tbm_f <- c( tbm_f, paste( "tbm$Motif_ID=='", motif_id, "'", sep="") );
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

	tbm <- tbm[nidx,, drop=F];
	cat( " ", NROW(tbm), "PWM(s) are defined in the CisBP dataset.\n");

	tmp.dir <- tempdir();
	pwm.files <- c();
	nidx.motif <- c();
	names <- c();
	
	err.missing <- 0;
	err.empty <- 0;
	
	for (i in 1:NROW(tbm) )
	{
		motif_id <- as.character(tbm$Motif_ID[i]);
		
		if( as.character(motif_id)==".")
		{
			# cat("! No ID for this motif(.).\n"); 	
			err.missing <- err.missing + 1;
			next;
		}

		pwm.file <- paste( "pwms_all_motifs/", motif_id, ".txt", sep="");
		r.file <- unzip( cisbp.db@zip.file, c(pwm.file), exdir=tmp.dir  );
		if(length(r.file)<1)
		{
			# cat("! Can not find PWM file for motif ID=", motif_id, ".\n" );
			err.missing <- err.missing + 1;
			next;
		}

		tb <- try( read.table(paste(tmp.dir, pwm.file, sep="/"), header=T, sep="\t", row.names=1), TRUE );
		if(class(tb)=="try-error") 
		{
			# cat("! Can not find PWM file for motif ID=", motif_id, ".\n" );
			err.missing <- err.missing + 1;
			next;
		}
	
		if(NROW(tb)==0)
		{
			# cat("! No A C G T values in the PWM file for motif ID=", motif_id, ".\n" );
			err.empty <- err.empty + 1;
			next;
		}

		pwm.files <- c( pwm.files, paste(tmp.dir, pwm.file, sep="/") );
		nidx.motif<- c( nidx.motif, i );
		names     <- c( names, motif_id);
	}

	TF_info <- NULL;
	TF_missing <- NULL;
	if(length(nidx.motif)>0)
	{
		TF_info <- tbm[ nidx.motif, , drop=F];
		rownames(TF_info)<- c(1:NROW(TF_info));
		TF_missing <- tbm[ -nidx.motif, , drop=F];
		rownames(TF_missing)<- c(1:NROW(TF_missing));
	}
	
	if( err.empty + err.missing > 0 )
		cat("!", err.empty + err.missing, "PWM file(s) are failed to be loaded ( Missing PWMs :", err.missing, ", Empty PWMs :", err.empty, ").\n");

	if(length(pwm.files)==0)
	{
		cat("! No PWM file to create a tfbs object.\n");
		return(NULL);
	}
	else
		cat("*", length(pwm.files), "PWM(s) in the tfbs object.\n");

	tfs <- tfbs(
		species    = cisbp.db@species,
		filenames  = pwm.files, 
		names      = names, 
		tf_info    = TF_info,
		tf_missing = TF_missing,
		header=T, sep="\t" , row.names=1 );


	return(tfs);
}
