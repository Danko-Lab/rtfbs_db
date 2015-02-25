get_reads_from_bigwig <- function( bw.plus, bw.minus, twoBit_path )
{
	file.tmp <- tempfile();
		
	system(paste("twoBitInfo", twoBit_path, file.tmp, sep=" "));
	chromInfo <- read.table( file.tmp );
	unlink( file.tmp );

	offset_dist <- 250;
	chromInfo <- chromInfo[grep("_|chrM|chrY|chrX", chromInfo[,1], invert=TRUE),]
	chromInfo <- data.frame(chrom=chromInfo[,1], chromStart=rep(0)+offset_dist, chromEnd=(chromInfo[,2]-1-offset_dist))
	
	r.bed <- data.frame( seqnames=chromInfo[,1],
	  	starts=chromInfo[,2],
	  	ends=chromInfo[,3],
	  	names=".",
	  	scores=".",
	  	strands="+")
	  	
   r.plus <- sum(bed6.region.bpQuery.bigWig( bw.plus, bw.minus, r.bed));

   r.bed <- data.frame( seqnames=chromInfo[,1],
	  	starts=chromInfo[,2],
	  	ends=chromInfo[,3],
	  	names=".",
	  	scores=".",
	  	strands="-")

   r.minus <- sum(bed6.region.bpQuery.bigWig( bw.plus, bw.minus, r.bed));

   return(c(r.plus,r.minus));	  	
}

simple_reduce_bed<-function( r.bed )
{
		if ( !( all( r.bed$strand == "+" ) || all( r.bed$strand == "-" ) ) )
		{
			show(r.bed);
			warning("The BED ranges with two different strands will be forced to merge.");
		}
		
		if ( length(unique(r.bed$chr))>1 )
		{
			warning("The BED ranges in two more chromosomes will be forced to merge.");
		}

		r.reduce <- c(r.bed$chr[1], min(r.bed$start), max(r.bed$end), r.bed$id[1], r.bed$score[1], r.bed$strand[1]);
		return(r.reduce); 		
}

#' Gets expression level of target TF.
#' USE extra_info$DBID to find gene information encoded by GENCODE V21
#'
#' @param tfbs: tfbs object
#' @param bed:  NOT USED
#' @param file_bigwig_plus:  bigwig file for strand plus(+)
#' @param file_bigwig_minus:  bigwig file for strand minus(-)
#' @param twoBit_path:  hg19 Human Genome or other species
#'
#' @return: tfbs.db object with changed expressionlevel;

tfbs_getExpression <- function(tfbs, bed, file_bigwig_plus, file_bigwig_minus, twoBit_path) 
{
	stopifnot(!is.null(tfbs@extra_info));
    
    # load pre-installed database(GENCODE HUMAN V21)
    load( system.file("extdata", "gencode_human21_transcript_ext.rdata", package="rtfbsdb"), environment() );
    
    bw.plus  <- load.bigWig( file_bigwig_plus );
    bw.minus <- load.bigWig( file_bigwig_minus ) ;
    bw.reads <- get_reads_from_bigwig( bw.plus, bw.minus, twoBit_path);
    cat("*", sum(bw.reads), "Reads in", file_bigwig_plus, "and", file_bigwig_minus,"\n");
	
	r.bed.list <- c();
	r.dbid.idx <- c();
	
	for( i in 1:NROW(tfbs@extra_info))
	{
	    dbid <- as.character(tfbs@extra_info$DBID[i]);
    
		r.bed.idx <- which(gencode_human21_transcript_ext$gene_id==dbid);
	    if(length(r.bed.idx)<1)
	    {
	    	cat ("Can't find DBID=", dbid, ", missing data.\n");
	    	next;
	    }
    
    	r.bed <- gencode_human21_transcript_ext[r.bed.idx, c(1,4,5,6,9,7), drop=F ];
	    colnames(r.bed) <- c('chr','start','end','id','score','strand');
	
	    # Here is a easy method to merege a BED data.frame, but it depends on other library.
	    # library(GenomicRanges)
	    # g.bed <- with(r.bed, GRanges(chr, IRanges(start, end), strand, score, id=id))
		# r.bed.rudece <- reduce(g.bed);
		
		r.bed.rudece <- simple_reduce_bed( r.bed );

	    r.bed.list  <- rbind( r.bed.list,  
	    	 c( seqnames = r.bed.rudece[1],
		  		starts   = as.numeric(r.bed.rudece[2])-1,
		  	  	ends     = as.numeric(r.bed.rudece[3]),
		  		names    = ".",
		  		scores   = ".",
		  		strands  = r.bed.rudece[6] ) );
		  
	    r.dbid.idx <- c(r.dbid.idx, i); 
	}

	r.bed.list <- as.data.frame(r.bed.list);
	r.bed.list$starts <- as.numeric(as.character(r.bed.list$starts))
	r.bed.list$ends <- as.numeric(as.character(r.bed.list$ends))
	
	r.reads  <- bed6.region.bpQuery.bigWig( bw.plus, bw.minus, r.bed.list );
    r.lambda <- 0.04 * sum(bw.reads)/10751533/1000 ;

	p.pois <- ppois( abs(r.reads), abs( r.lambda*( r.bed.list[,3] - r.bed.list[,2] ) ), lower.tail=F);
	
	df.exp <- data.frame( 
				Motif_ID = tfbs@extra_info$Motif_ID,
				DBID     = tfbs@extra_info$DBID, 
				chr      = NA, 
				start    = NA, 
				end      = NA, 
				strand   = NA, 
				bed_length= NA,
				reads    = NA,
				lambda   = NA,
				prob     = NA);
				
	df.exp$chr[r.dbid.idx]       <- r.bed.list[,1];
	df.exp$start[r.dbid.idx]     <- r.bed.list[,2];
	df.exp$end[r.dbid.idx]       <- r.bed.list[,3];
	df.exp$strand[r.dbid.idx]    <- as.character(r.bed.list[,6]);
	df.exp$bed_length[r.dbid.idx]<- abs( r.bed.list[,3] - r.bed.list[,2] );
	df.exp$reads[r.dbid.idx]     <- r.reads;
	df.exp$lambda[r.dbid.idx]    <- r.lambda;
	df.exp$prob[r.dbid.idx]      <- p.pois;
	
	# expressionlevel is data.frame including all information.
	tfbs@expressionlevel <- df.exp;
	
	return(tfbs);
}