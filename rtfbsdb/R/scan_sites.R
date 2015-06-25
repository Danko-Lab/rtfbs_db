## 
## find_sites_rtfbs.R
##
## Finds PWMs in DNAse-1 peaks using rtfbs package availiable on CRAN.
##

#' Extend BED file
#'
#' @param bed data.frame with bed regions
#' @param len extension in bp
#' @return extended bed data.frame
extend.bed <- function(bed, len) {
  starts = as.integer(bed[,2] - len)
  ends = as.integer(bed[,3] + len)
  
  N = dim(bed)[2]
  if (N == 3) {
    data.frame(bed[,1], starts, ends)
  } else {
    data.frame(bed[,1], starts, ends, bed[, 4:N])
  }
}

#` Returns the posterior probability of TF binding, conditional on the data passed as part of the function.
#`
#` Optionally takes additional parameters passed through to score.ms.
#`
#` @param tf_name name of the TF.
#` @param tre.bed bed-formatted peak information.
#` @param motif_path path to the motif PWM file.
#` @param divide_num List of parameters for all data types, for the model representing no TF binding.
#` @return List structure representing the match score to the motif.
scan_rtfbs <- function(tf_name, file.twoBit, tre.bed, motif_path, return.posteriors=TRUE, ...) {
  ## Read the pwm and sequence file.
  motif <- read.motif(motif_path, header=TRUE) 
  
  ## Write out new fasta file, adding on half-width of the motif, to correctly align motif center ...
  half_width <- ceiling( (NROW(motif)-1)/2 );
  extBed = extend.bed( tre.bed, half_width - 1 );
  dnase_peaks = read.seqfile.from.bed( extBed, file.twoBit );

  ## Swiched rtfbs to returning posteriors.
  bgModel <- build.mm(dnase_peaks, 3);

  ## Switch to rtfbsPost to return proper data structure.
  binding <- score.ms(dnase_peaks, motif, bgModel, return_posteriors=return.posteriors, ...) ;
  
  if(return.posteriors == FALSE) { ## Parse binding site into genomic coordinates.
    spl <- strsplit(as.character(binding$seqname), ":|-")
    peak_chrom <- as.character(sapply(c(1:NROW(binding)), function(x) {spl[[x]][[1]]}))
    peak_start <- as.integer(sapply(c(1:NROW(binding)), function(x) {spl[[x]][[2]]}))
	binding <- data.frame(chrom= peak_chrom, 
							chromStart= peak_start+ binding$start, 
							chromEnd= peak_start+ binding$end, 
							name= binding$seqname,
							score= binding$score,
							strand= binding$strand);
  }
  
  return(binding);
}

get_binding_site <- function( bgModel1, seq.ms, PWM, return.posteriors, score.threshold=6, fdr.threshold=NA, gc.groups=NA, background.order = 2, background.length = 100000)
{
    if( is.na( gc.groups) )	
    {
	    if(is.na(fdr.threshold))
            binding <- score.ms( seq.ms, PWM, bgModel1, return_posteriors=return.posteriors, threshold=score.threshold )
	    else
	    {
            seq.score <- score.ms( seq.ms, PWM, bgModel1, return_posteriors=return.posteriors, threshold=0 );	    
			simu.ms   <- simulate.ms( bgModel1, background.length );
			simu.score<- score.ms( simu.ms, PWM, bgModel1, threshold=0 );		
            fdrMap    <- calc.fdr( seq.ms, seq.score, simu.ms, simu.score );
            binding   <- output.sites( seq.score, fdrScoreMap  = fdrMap, fdrThreshold = fdr.threshold);
	    }	    	
    }
    else
    {
		if( is.na(score.threshold) ) score.threshold <- 6;

        msGroups <- groupByGC.ms( seq.ms, gc.groups);

		bgModels <- lapply(1:length(msGroups), 
						function(i) { build.mm(msGroups[[i]], background.order) } );

		seq.score <- lapply(1:length(msGroups),
						function(i) { score.ms(msGroups[[i]], PWM, bgModels[[i]], return_posteriors=return.posteriors, threshold=ifelse(!is.na(fdr.threshold), 0, score.threshold));});

		if(!is.na(fdr.threshold))
		{
			simu.ms <- lapply( 1:length(msGroups), 
							function(i){ simulate.ms(bgModels[[i]], background.length)});
			simu.score <- lapply( 1:length(msGroups), 
							function(i){ score.ms(simu.ms[[i]], PWM, bgModels[[i]])});
			fdrMap  <- lapply( 1:length(msGroups), 
							function(i) { calc.fdr(msGroups[[i]], seq.score[[i]], simu.ms[[i]], simu.score[[i]])});
			binding <- lapply( 1:length(msGroups), 
							function(i) { output.sites(seq.score[[i]], fdrScoreMap  = fdrMap[[i]], fdrThreshold = fdr.threshold);} );

			binding <- do.call("rbind", binding);
		}
		else
		{
			binding <- do.call("rbind", seq.score);
		}
	
	}
	
	return( binding );
}

## return.type --> c("matches", "posteriors", "maxposterior", or "writedb")
##
## matches 		-- returns all matching motifs.
## writedb 		-- writes a bed file with matches.  Assuems that sort-bed and starch tools are availiable in $PATH
## posteriors 	-- returns the posteriors at each position.
## maxposterior	-- returns the max(posterior) in each dnase-1 peak.

scanDb_rtfbs <- function(tfbs, file.twoBit, tre.bed, return.type = "matches", file.prefix = NA, usemotifs = NA, ncores = 3, fdr.threshold = NA, score.threshold = 6, gc.groups = NA, background.order = 2, background.length = 100000,...) {
  stopifnot(class(tfbs) == "tfbs")

  if( !is.na(file.prefix))
  	if( !check_folder_writable( file.prefix ) ) 
  	  stop(paste("Can not create files starting with the prefix:", file.prefix));

  ## Read in the DNAse-1 peaks ...
  half_width=15 ## Max size of TF in set of 1800 is 30 (half-width = 15).
  options("scipen"=100, "digits"=4)

  extBed  <- extend.bed( tre.bed, half_width - 1)
  seq.ms  <- read.seqfile.from.bed( extBed, file.twoBit);
  bgModel <- build.mm( seq.ms, 3);

  ## Swiched rtfbs to returning posteriors.
  return.posteriors <- ( return.type=="posteriors" | return.type=="maxposterior" )

  binding_all <- mclapply(usemotifs, function(i, ...) {
	  
	  binding<- NULL;
	  
	  ## Read the pwm and sequence file.
	  PWM <- tfbs@pwm[[i]];
		
	  ## Parse binding site into genomic coordinates.	
      if( return.type == "maxposterior" || return.type == "posteriors") 
      { 
         binding <- score.ms( seq.ms, PWM, bgModel, return_posteriors=TRUE, threshold = score.threshold );

    	 if( return.type == "maxposterior"&& NROW(binding) > 0)
	    	 return(sapply(1:NROW(binding), function(x) { 
	        	max(c(binding[[x]]$MotifModel$Forward - binding[[x]]$Background, binding[[x]]$MotifModel$Reverse - binding[[x]]$Background)) }));
      }

	  ## Parse binding site into genomic coordinates.
	  if(return.posteriors == FALSE ) 
	  { 
		 binding <- get_binding_site( bgModel, seq.ms, PWM, FALSE, score.threshold = score.threshold, fdr.threshold = fdr.threshold, gc.groups=gc.groups, background.order = 2, background.length = 100000 );

         if( NROW(binding) > 0 )
         {
			 spl <- strsplit(as.character(binding$seqname), ":|-")
			 peak_chrom <- as.character(sapply(c(1:NROW(binding)), function(x) {spl[[x]][[1]]}))
			 peak_start <- as.integer(sapply(c(1:NROW(binding)), function(x) {spl[[x]][[2]]}))

			 binding <- data.frame(chrom= peak_chrom, 
									chromStart= peak_start+ binding$start- 1,  ## -1 determined empirically.
									chromEnd= peak_start+ binding$end, 
									name= tfbs@mgisymbols[i], # binding$motif_id
									score= binding$score,
									strand= binding$strand)

			 if(return.type == "writedb") {
				write.table(binding, file= pipe(paste(" sort-bed - | starch - > ",file.prefix,i,".bed.tmp.starch", sep="")),  quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
				return(paste(file.prefix,i,".bed.tmp.starch ", sep=""))
			 }
		 }
	  }
	  
	  return(binding);
  }, mc.cores= ncores)

  if(return.type == "writedb") 
  {
     cat_files <- paste(unlist(binding_all), collapse=" ")
     system(paste("starchcat ", cat_files, " > ",file.prefix,".db.starch", sep=""))
     system(paste("rm ",file.prefix,"*.bed.tmp.starch",sep=""))
     binding_all <- paste(file.prefix,".db.starch", sep="")
  }
  
  if(return.type == "maxposterior") 
  {
     binding_all <- matrix(unlist(binding_all), nrow= NROW(tre.bed), ncol= NROW(usemotifs))
  }

  return(binding_all)
}

# ncores=3 for 4 cores CPU.

tfbs_scanTFsite<-function( tfbs, file.twoBit, tre.bed = NULL, return.type="matches", file.prefix = NA,  usemotifs = NA, ncores = 3, fdr = NA, threshold = 6, gc.groups = NA, background.order = 2, background.length = 100000 )
{
    stopifnot(class(tfbs) == "tfbs")

	if( missing(tre.bed) )
	{
		chromInfo <- get_chromosome_size( file.twoBit );

		offset_dist <- 250;
		chromInfo <- chromInfo[grep("_|chrM|chrY|chrX", chromInfo[,1], invert=TRUE),];
		tre.bed <- data.frame(chrom=chromInfo[,1], chromStart=rep(0)+offset_dist, chromEnd=(chromInfo[,2]-1-offset_dist));
	}
	
	if( missing(file.prefix) ) file.prefix="scan.db";
	if( missing(return.type) ) return.type="matches";
	if( missing(ncores) ) ncores= 3;
	if( missing(threshold) && missing(fdr) ) threshold= 6;
	if( missing(usemotifs)) usemotifs =c(1:tfbs@ntfs);

	r.ret <- scanDb_rtfbs( tfbs, file.twoBit, tre.bed, file.prefix = file.prefix, return.type = return.type, usemotifs = usemotifs, ncores = ncores, 
	                       fdr.threshold = fdr, score.threshold = threshold, gc.groups = gc.groups, background.order = background.order, background.length = background.length ); 
	
	r.parm <- list(file.twoBit = file.twoBit, 
				  file.prefix = file.prefix,  
				  return.type = return.type, 
				  usemotifs   = usemotifs, 
				  ncores      = ncores,  
				  fdr         = fdr, 
				  threshold   = threshold, 
				  gc.groups   = gc.groups, 
				  background.order = background.order, 
				  background.length = background.length);
	
	
	sum.match <- NULL;
	if( return.type == "matches" )
	{
		sum.match <- do.call("rbind", lapply(r.ret, function(x){ 
				tf.name <- "";
				tf.idx  <- which( as.character(tfbs@extra_info$Motif_ID) == as.character(x$name[1]) );
				if( length(tf.idx)>0)
					tf.name <- tfbs@extra_info$TF_Name[tf.idx[1]];
				return( data.frame( x$name[1], tf.name, NROW(x) ) );
		} ) );
		
		colnames(sum.match) <- c("TF_Name", "Motif_ID", "Count");
	}

	
	r.scan <- list( parm = r.parm, bed = tre.bed, result = r.ret, summary=sum.match );
	class( r.scan ) <- c( class(r.scan), "tfbs.finding");
	
	return( r.scan );
}

print.tfbs.finding<-function(x, ...)
{
	r.scan <- x;
	 
	cat("Return type: ", r.scan$parm$return.type, "\n");
	cat("FDR threshold: ", r.scan$parm$fdr, "\n");
	cat("Score threshold: ", r.scan$parm$threshold, "\n");

	if( r.scan$parm$return.type == "matches" )
	{
		df.allfinding <- do.call("rbind", r.scan$result );
		cat("Motifs count: ", length(r.scan$result), "\n");
		cat("Binding sites: ", NROW(df.allfinding), "\n");

		summary <- r.scan$summary[order(r.scan$summary[,3], decreasing = TRUE ),];
		if( NROW( summary )>20)
			summary <- summary[c(1:20),];
		show( summary );		
	}
	
	if( r.scan$parm$return.type == "writedb" )
		cat("Binary Bed file: ", r.scan$result, "\n");
		
	if( r.scan$parm$return.type == "posteriors" )
	{
		df.allfinding <- do.call("rbind", r.scan$result );
		cat("Motifs count: ", length(r.scan$result), "\n");
		cat("Binding sites: ", NROW(df.allfinding), "\n");
	}
	
	if(r.scan$parm$return.type=="maxposterior")
		cat("Matrix posterior: ", NROW(r.scan$result), "*", NCOL(r.scan$result), "\n");
}

summary.tfbs.finding<-function( object, ... )
{
	r.scan <- object;

	if( r.scan$parm$return.type == "matches" )
		return( r.scan$summary )
	else	
		return(NULL);
}

tfbs.reportFinding<-function( tfbs, r.scan, file.pdf = NA, report.size = "letter", report.title = "" )
{
	if( r.scan$parm$return.type != "matches" )
		cat( "! No summary information for the report.\n" )
	else	
	{
		summary <- r.scan$summary[order(r.scan$summary[,3], decreasing = TRUE ),];

		r.scan.sum <- data.frame( No=c(1:NROW(summary)), summary, summary[,1] );

		df.style <- data.frame(position=numeric(0), width=numeric(0), header=character(0), hjust=character(0), style=character(0), extra1=character(0), extra2=character(0), extra2=character(0));
		df.style <- rbind(df.style, data.frame( position=0.00, width=0.04, header="No.",       hjust="left",   style="text", extra1="0",  extra2="0",  extra3="0" ) );
		df.style <- rbind(df.style, data.frame( position=0.04, width=0.10, header="Motif ID",  hjust="left",   style="text", extra1="0",  extra2="0",  extra3="0" ) );
		df.style <- rbind(df.style, data.frame( position=0.14, width=0.10, header="TF Name",   hjust="left",   style="text", extra1="0",  extra2="0",  extra3="0" ) );
		df.style <- rbind(df.style, data.frame( position=0.24, width=0.10, header="Count",     hjust="centre", style="text", extra1="0",  extra2="0",  extra3="0" ) );
		df.style <- rbind(df.style, data.frame( position=0.34, width=0.49, header="Motif Logo",hjust="centre", style="logo", extra1="0",  extra2="0",  extra3="0" ) );

		output_motif_report( tfbs, r.scan.sum, file.pdf, report.size, report.title, df.style );
	}
}

