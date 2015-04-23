# save the chromosme size into one ftable file
#
# save_chromosome_size("mm10.2bit", "chrom_info_mm19.tab")
#
save_chromosome_size <- function( file.twoBit, file.tab)
{
	write.table(get_chromosome_size(file.twoBit), file=file.tab, row.names=F, col.names=F, quote=F);
}


# get chromosome size based on genemo data of 2 bit file.
#
# @file.twoBit, e.g. hh19.2bit, mm10.2bit
#
get_chromosome_size <- function(file.twoBit)
{
	file.tmp <- tempfile();

	system(paste("twoBitInfo", file.twoBit, file.tmp, sep=" "));
	chromInfo <- read.table( file.tmp );
	unlink( file.tmp );

	return(chromInfo);
}

# query read counts of all chromosomes from bigWig file.
#
# @bw.plus,     bigWig object
# @bw.minus,    bigWig object
# @chromInfo    data.frame with 2 columns(chr, size);
#
# @return       vector of reads in plus and minus file.
get_reads_from_bigwig <- function( bw.plus, bw.minus, chromInfo )
{
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
#' @param file.bigwig.plus:  bigwig file for strand plus(+)
#' @param file.bigwig.minus:  bigwig file for strand minus(-)
#' @param file.twoBit:  hg19 Human Genome or other species
#                       hg19.2bit: contains the complete hg19 Human Genome in the 2bit format
#                                  http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit(08-Mar-2009 15:29)
#                       mm10.2bit  http://hgdownload-test.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit(07-Feb-2012 10:52)  
#' @param gencode.ext.rdata:  Genecode RDATA file Encoded by import_gencode
#'
#' @return: tfbs.db object with changed expressionlevel;

tfbs_getExpression <- function(tfbs, file.bigwig.plus, file.bigwig.minus, file.twoBit=NA, gencode.ext.rdata=NA) 
{
	stopifnot(!is.null(tfbs@extra_info));
    
    # load gencode RDATA file and set the table of gencode_transcript_ext
    gencode_transcript_ext <- NULL;
    if( missing(gencode.ext.rdata) )
    {
	    # load pre-installed database(GENCODE HUMAN V21)
	    if (tfbs@species=="Homo_sapiens" || tfbs@species=="human" || tfbs@species=="Human" )
	    	load( system.file("extdata", "gencode_v21_transcript_ext.rdata", package="rtfbsdb"), environment() )
	    # load pre-installed database(GENCODE MOUSE V3)
	    else if (tfbs@species=="Mus_musculus" || tfbs@species=="mouse"  || tfbs@species=="Mouse" )
	    	load( system.file("extdata", "gencode_vM3_transcript_ext.rdata", package="rtfbsdb"), environment() )
	    else
	    	stop("The tfbs object is not calculated from human or mouse data, you need to provide the gencode data for this species.");
    }
    else
    {
	    load( gencode.ext.rdata, environment() );
    	if(is.null(gencode_transcript_ext))
    		stop("Gencode data can not be found in the RDATA file specified by the parameter of gencode.ext.rdata.");
    }
    
    # Load chromosome sizes or calculate the size using file.twoBit
    chromInfo <- NULL;
    if(missing(file.twoBit) )
    {
    	if(tfbs@species=="Homo_sapiens" )
			chromInfo <- read.table ( system.file("extdata", "chrom_info_hg19.tab", package="rtfbsdb"), header=F )
    	else if(tfbs@species=="Mus_musculus")
			chromInfo <- read.table ( system.file("extdata", "chrom_info_mm19.tab", package="rtfbsdb"), header=F )
		else
			warning("The p-value can not be calculated because of no 2bit file. ", call. = FALSE, immediate. = TRUE);
    }
    else
		chromInfo <- get_chromosome_size(file.twoBit);

	# Load bigWig files(minus and plus)
    bw.plus  <- try( load.bigWig( file.bigwig.plus ) );
    bw.minus <- try( load.bigWig( file.bigwig.minus ) );
    if( class(bw.plus)=="try-error" || class(bw.minus)=="try-error" )
    	stop("Failed to load bigwig files.");
   
	# Count the rads for all chromosomes if possible
    bw.reads <- NULL;
	if(!is.null(chromInfo))
	   	bw.reads  <- get_reads_from_bigwig( bw.plus, bw.minus, chromInfo);
   	cat("*", sum(abs(bw.reads)), "Reads in", file.bigwig.plus, "and", file.bigwig.minus,"\n");
	
	# Count the rads for all chromosomes if possible
	r.bed.list <- c();
	r.dbid.idx <- c();
	
	# make a union range for each motif
	for( i in 1:NROW(tfbs@extra_info))
	{
	    dbid <- as.character(tfbs@extra_info$DBID[i]);
    
		r.bed.idx1 <- which(gencode_transcript_ext$gene_id==dbid);
		r.bed.idx2 <- grep( paste(dbid, ".", sep=""), gencode_transcript_ext$gene_id );
		r.bed.idx <- unique( c(r.bed.idx1 ,r.bed.idx2) );

	    if(length(r.bed.idx)<1)
	    {
	    	cat ("Can't find DBID=", dbid, ", missing data.\n");
	    	next;
	    }
    
    	r.bed <- gencode_transcript_ext[r.bed.idx, c(1,4,5,6,9,7), drop=F ];
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
	
	# Query reads for each motif
	r.reads  <- bed6.region.bpQuery.bigWig( bw.plus, bw.minus, r.bed.list );

	r.lambda <- NA;
	p.pois   <- NA;
	if(!is.null(bw.reads))
	{
		r.lambda <- 0.04 * sum(abs(bw.reads))/10751533/1000 ;
		p.pois <- ppois( abs(r.reads), abs( r.lambda*( r.bed.list[,3] - r.bed.list[,2] ) ), lower.tail=F);
	}
	
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

#
# Extract the gencode information and save it into RDATA file
#
# e.g. 
# import_gencode( "human", "gencode.v22.annotation.gtf", "gencode.v22.rdata")

import_gencode <-function( species, file.gtf, file.rdata, ncores = 1)
{
	if(!require(sqldf))
		stop("Package sqldf is necessary to do data operation, please install it.");

	library(parallel)
	
	f <- file( file.gtf)
	temp.db <-tempfile() 
	
	bigdf <- sqldf("select * from f", dbname = temp.db, 
				   file.format = list(header = F, row.names = F, sep="\t", skip=6));
	
	ids.multi <- strsplit(bigdf$V9, ";");

	cat(length(ids.multi), "Items are loaded, starting parse...\n");
	
	g.names <- c();
	dfs <- mclapply( ids.multi, function(x) {
		x2<-strsplit(x, " ");
		lnames <- c()
		values <-c();
		for(i in 1:length(x2) )
		{
			name.f <- 0;
			values.f <- c();
			for(j in 1:length(x2[[i]]))
			{
				if (x2[[i]][j]=="") next;
				if(name.f==0)
				{
					lnames<-c(lnames, x2[[i]][j]);
					name.f <-1;
				}
				else
					values.f <- c(values.f, x2[[i]][j]);
			}
			
			values<- c( values, paste(values.f, collapse=" ") );
		}
	
		df <- t(data.frame(values));
		colnames(df) <- lnames;
		
		g.names <<- c( g.names, lnames);

		df;},  mc.cores = ncores )
	
	g.names <-unique(g.names);	
	
	cat("Parsing is finished.\n");

	select_value<-function(dfs, tag.id)
	{
		tags.value <- lapply(dfs, function(x) {
			idx <- which( colnames(x)==tag.id);
			if (length(idx)<1) return(NA) else return(x[1, idx[1]]);
		})

		tags.value <- gsub("\"", "", tags.value);
		tags.value;
	}

	df.ext <- data.frame(
		gene_id                  = select_value( dfs, "gene_id"), 
		transcript_id            = select_value( dfs, "transcript_id"), 
		gene_type                = select_value( dfs, "gene_type"), 
		gene_status              = select_value( dfs, "gene_status"), 
		gene_name                = select_value( dfs, "gene_name"),
		transcript_type          = select_value( dfs, "transcript_type"),
		transcript_status        = select_value( dfs, "transcript_status"),
		transcript_name          = select_value( dfs, "transcript_name"),
		level                    = select_value( dfs, "level"),      
		tag                      = select_value( dfs, "tag"),      
		transcript_support_level = select_value( dfs, "transcript_support_level"),
		havana_gene              = select_value( dfs, "havana_gene"),
		havana_transcript        = select_value( dfs, "havana_transcript"),
		ont                      = select_value( dfs, "ont"),
		protein_id               = select_value( dfs, "protein_id"),
		ccdsid                   = select_value( dfs, "ccdsid") );


	gencode_transcript_ext <- data.frame(bigdf[,c(1:8)], df.ext);

	save(gencode_transcript_ext, file=file.rdata);

	cat("Outputing RDATA file:", file.rdata, "\n");
}