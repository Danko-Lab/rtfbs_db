

#' Find the subset by querying the motif table
#'
#' @param motifDB:        MotifList, MotifDb or subset of MotifDb
#' @param organism:       string, the meta filed in the MotifDb
#' @param geneSymbol:     string, the meta filed in the MotifDb
#' @param tfFamily:       string, the meta filed in the MotifDb
#' @param providerName:   string, the meta filed in the MotifDb
#' @param providerId:     string, the meta filed in the MotifDb
#' @param dataSource:     string, the meta filed in the MotifDb
#' @param geneId:         string, the meta filed in the MotifDb
#' @param geneIdType:     string, the meta filed in the MotifDb
#' @param proteinId:      string, the meta filed in the MotifDb
#' @param proteinIdType:  string, the meta filed in the MotifDb
#' @param sequenceCount:  string, the meta filed in the MotifDb
#' @param bindingSequence:string, the meta filed in the MotifDb
#' @param bindingDomain:  string, the meta filed in the MotifDb
#' @param experimentType: string, the meta filed in the MotifDb
#' @param pubmedID:       string, the meta filed in the MotifDb
#' @param pseudocount:    -7, replace -INF for PWM elements.
#'
#' @return: NULL or tfbs object;


tfbs.createFromMotifDb <- function ( 
					motifDB        = NULL,
					organism       = "Hsapiens",
					geneSymbol     = NULL,
					tfFamily       = NULL,
					providerName   = NULL,    
					providerId     = NULL,
					dataSource     = NULL,
					geneId         = NULL,
					geneIdType     = NULL,
					proteinId      = NULL,
					proteinIdType  = NULL,    
					sequenceCount  = NULL, 
					bindingSequence= NULL,
					bindingDomain  = NULL,
					experimentType = NULL,
					pubmedID       = NULL,
					pseudocount    = -7)
{
	if(is.null(motifDB))
	{
		if(!requireNamespace("MotifDb"))
			stop("The MotifDb package is not installed");
		
		motifDB <- MotifDb::MotifDb;
	}
	else
		stopifnot(class(motifDB)=="MotifList");
	
	tbm <- data.frame(motifDB@elementMetadata@listData, stringsAsFactors =FALSE);

	if( length(which(is.na(tbm$geneSymbol))) > 0 )
		tbm[is.na(tbm$geneSymbol),  "geneSymbol"] <- paste("ID_UNK_", 1:length(which(is.na(tbm$geneSymbol))), sep="");
	if( length(which(is.na(tbm$tfFamily))) > 0 )
		tbm[is.na(tbm$tfFamily),   "tfFamily"]   <- paste("TF_UNK_", 1:length(which(is.na(tbm$tfFamily))), sep="");
 	
	tbm_f <- c();
	if( !is.null(organism       )) tbm_f <- c( tbm_f, paste( "tbm$organism=='",        organism, "'", sep="") );
	if( !is.null(geneSymbol     )) tbm_f <- c( tbm_f, paste( "tbm$geneSymbol=='",      geneSymbol, "'", sep="") );
	if( !is.null(tfFamily       )) tbm_f <- c( tbm_f, paste( "tbm$tfFamily=='",        tfFamily, "'", sep="") );
	if( !is.null(providerName   )) tbm_f <- c( tbm_f, paste( "tbm$providerName=='",    providerName, "'", sep="") );
	if( !is.null(providerId     )) tbm_f <- c( tbm_f, paste( "tbm$providerId=='",      providerId, "'", sep="") );
	if( !is.null(dataSource     )) tbm_f <- c( tbm_f, paste( "tbm$dataSource=='",      dataSource, "'", sep="") );
	if( !is.null(geneId         )) tbm_f <- c( tbm_f, paste( "tbm$geneId=='",          geneId, "'", sep="") );
	if( !is.null(geneIdType     )) tbm_f <- c( tbm_f, paste( "tbm$geneIdType=='",      geneIdType, "'", sep="") );
	if( !is.null(proteinId      )) tbm_f <- c( tbm_f, paste( "tbm$proteinId=='",       proteinId, "'", sep="") );
	if( !is.null(proteinIdType  )) tbm_f <- c( tbm_f, paste( "tbm$proteinIdType=='",   proteinIdType, "'", sep="") );
	if( !is.null(sequenceCount  )) tbm_f <- c( tbm_f, paste( "tbm$sequenceCount=='",   sequenceCount, "'", sep="") );
	if( !is.null(bindingSequence)) tbm_f <- c( tbm_f, paste( "tbm$bindingSequence=='", bindingSequence, "'", sep="") );
	if( !is.null(bindingDomain  )) tbm_f <- c( tbm_f, paste( "tbm$bindingDomain=='",   bindingDomain, "'", sep="") );
	if( !is.null(experimentType )) tbm_f <- c( tbm_f, paste( "tbm$experimentType=='",  experimentType, "'", sep="") );
	if( !is.null(pubmedID       )) tbm_f <- c( tbm_f, paste( "tbm$pubmedID=='",        pubmedID, "'", sep="") );

	if(length(tbm_f)>0)
	{
		tbm_f_all <- paste( tbm_f, collapse=" & " );
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

	pwms <- lapply(nidx, function(i) {
		pfm <- t( motifDB@listData[[i]] ) ; 
		pwm <- log(pfm/rowSums(pfm))
		pwm [ pwm==-Inf ] <- pseudocount;
		return(pwm); 
	});

	tbm <- data.frame(Motif_ID=tbm$geneSymbol, TF_Name=tbm$tfFamily, tbm);
	
	tfs <- new("tfbs", 
		species        = as.character(if(is.null(organism)) "ALL" else organism),
		ntfs           = as.integer(NROW(tbm)),
		pwm            = pwms,
		filename       = as.character(""), 
		mgisymbols     = as.character(tbm$Motif_ID), 
		tf_info        = tbm,
		tf_missing     = as.data.frame(NULL),
		distancematrix = matrix(, nrow=0, ncol=0),
		cluster        = matrix(, nrow=0, ncol=0),
		expressionlevel= as.data.frame(NULL) )

	return(tfs);
}
