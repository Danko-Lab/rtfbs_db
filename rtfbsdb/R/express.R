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
	r.bed <- r.bed[ order( r.bed$strand, r.bed$chr, r.bed$start, r.bed$end),,drop=F]
	
	merged <- TRUE;
	while(merged)
	{
		merged <- FALSE;
		if(NROW(r.bed)>1)
		for(i in 1:(NROW(r.bed)-1))
		{
			if( r.bed[i,1] == r.bed[i+1,1] && r.bed[i,6] == r.bed[i+1,6] )
			{
				# inclusion
				if( r.bed[i,2] <= r.bed[i+1,2] && r.bed[i+1,3] <= r.bed[i,3] )
				{
					r.bed <- r.bed[-(i+1),,drop=F]
					merged = TRUE;
					break;
				}	
				#intersect
				else if( r.bed[i,2] <= r.bed[i+1,2] && r.bed[i+1,2] <= r.bed[i,3] && r.bed[i,3] <= r.bed[i+1,3] )
				{
					r.bed[i,3] <- r.bed[i+1,3]
					r.bed <- r.bed[-(i+1),,drop=F]
					merged = TRUE;
					break;
				}
				#Adjacent
				else if( r.bed[i,3] == r.bed[i+1,2] )
				{
					r.bed[i,3] <- r.bed[i+1,3]
					r.bed <- r.bed[-(i+1),,drop=F]
					merged = TRUE;
					break;
				}
			}		
		}
	}
	
	return(r.bed); 		
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
#' @param file.gencode.gtf:  Genecode GTF file downloaded from GENCODE website
#' @param seq.datatype:  three values: GRO-seq, PRO-seq, RNA-seq
#'
#' @return: tfbs.db object with changed expressionlevel;

tfbs_getExpression <- function(tfbs, file.bigwig.plus, file.bigwig.minus, file.twoBit=NA, file.gencode.gtf=NA,  seq.datatype=NA, ncores = 3 ) 
{
	stopifnot(!is.null(tfbs@extra_info));
    
    # load gencode RDATA file and set the table of gencode_transcript_ext
    gencode_transcript_ext <- NULL;
    if( missing(file.gencode.gtf) || is.na(file.gencode.gtf) )
    {
	    # load pre-installed database(GENCODE HUMAN V21)
	    if (tfbs@species=="Homo_sapiens" || tfbs@species=="human" || tfbs@species=="Human" )
	    {
	    	cat("* For Homo_sapiens species, the pre-installed Gencode V21 is used to find expresed TFs.\n");
	    	load( system.file("extdata", "gencode_v21_transcript_ext.rdata", package="rtfbsdb"), environment() )
	    }
	    # load pre-installed database(GENCODE MOUSE V3)
	    else if (tfbs@species=="Mus_musculus" || tfbs@species=="mouse"  || tfbs@species=="Mouse" )
	    {
	    	cat("* For Mus_musculus species, the pre-installed Gencode vM3 is used to find expresed TFs.\n");
	    	load( system.file("extdata", "gencode_vM3_transcript_ext.rdata", package="rtfbsdb"), environment() )
	    }
	    else
	    	stop("The tfbs object is not calculated from human or mouse data, you need to provide the gencode data for this species.");
    }
    else
    {
		gencode_transcript_ext <- try( import_gencode( tfbs@species, file.gencode.gtf, ncores = ncores) );
    	if(is.null(gencode_transcript_ext) || class(gencode_transcript_ext)=="try-error")
    		stop("Gencode data can not be found in the GTF file specified by the parameter of file.gencode.gtf.");

		#save(gencode_transcript_ext, file=file.rdata);
	    #load( gencode.ext.rdata, environment() );
		#cat("Outputing RDATA file:", file.rdata, "\n");
    }

    if(seq.datatype=="GRO-seq" ||seq.datatype=="PRO-seq")
    	gencode_transcript_ext <- gencode_transcript_ext[ which(gencode_transcript_ext$V3=="transcript"),];	
	
    if(seq.datatype=="RNA-seq")
    	gencode_transcript_ext <- gencode_transcript_ext[ which(gencode_transcript_ext$V3=="exon"),];	

   	cat("  For", seq.datatype, ",", NROW(gencode_transcript_ext), "items are selected from GENCODE dataset.\n");
    
    # Load chromosome sizes or calculate the size using file.twoBit
    chromInfo <- NULL;
    if( missing(file.twoBit) || is.na(file.twoBit) )
    {
    	if(tfbs@species=="Homo_sapiens" )
	    {
	    	cat("* The pre-installed chrmosome information of hg19 is loaded.\n");
			chromInfo <- read.table ( system.file("extdata", "chrom_info_hg19.tab", package="rtfbsdb"), header=F )
    	}
    	else if(tfbs@species=="Mus_musculus")
	    {
	    	cat("* The pre-installed chrmosome information of mm10 is loaded.\n");
			chromInfo <- read.table ( system.file("extdata", "chrom_info_mm10.tab", package="rtfbsdb"), header=F )
		}
		else
		{
			cat("! The p-value can not be calculated without 2bit file for ", tfbs@species, ".\n");
			stop("The p-value can not be calculated without 2bit file. ");
		}			
    }
    else
	{
		chromInfo <- try(get_chromosome_size(file.twoBit));
		if( class(chromInfo)=="try-error")
			stop("Failed to use twoBitInfo command to get the chromosome information from 2bit file.");
	}	

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

	r.lambda <- NA;
	if(!is.null(bw.reads))
		r.lambda <- 0.04 * sum(abs(bw.reads))/10751533/1000 ;

	DBIDs <- unique( as.character(tfbs@extra_info$DBID) );
	
if(0)
{
	# Count the reads for all chromosomes if possible
	if(length(DBIDs)>0)
	r.bed.list <- mclapply( 1:length(DBIDs), function(i) {
	    dbid <- as.character(DBIDs[i]);
    
		r.bed.idx1 <- which(gencode_transcript_ext$gene_id==dbid);
		r.bed.idx2 <- grep( paste(dbid, ".", sep=""), gencode_transcript_ext$gene_id );
		r.bed.idx <- unique( c(r.bed.idx1 ,r.bed.idx2) );

	    if(length(r.bed.idx)<1)
	    {
	    	cat ("Can't find DBID=", dbid, ", missing data.\n");
	    	return(c(dbid=dbid, rep(NA, 8)));
	    }
    
    	r.bed <- gencode_transcript_ext[r.bed.idx, c(1,4,5,6,9,7), drop=F ];
	    colnames(r.bed) <- c('chr','start','end','id','score','strand');
	
	    # Here is a easy method to merege a BED data.frame, but it depends on other library.
	    # library(GenomicRanges)
	    # g.bed <- with(r.bed, GRanges(chr, IRanges(start, end), strand, score, id=id))
		# r.bed.rudece <- reduce(g.bed);
		
		r.bed.rudece <- simple_reduce_bed( r.bed );

		# Query reads for each motif
		r.reads  <- bed6.region.bpQuery.bigWig( bw.plus, bw.minus, r.bed.rudece );

		p.pois <- ppois( abs(r.reads), abs( r.lambda * sum( r.bed.rudece[,3] - r.bed.rudece[,2] ) ), lower.tail=F);

	    r.df   <- c( 
	    		dbid		= dbid,
	    		chr 		= r.bed.rudece[1, 1],
		  		start   	= min(r.bed.rudece [,2]),
		  	  	end     	= max(r.bed.rudece [,3]),
		  	  	length     	= sum(r.bed.rudece[,3] - r.bed.rudece[,2] ),
		  		strand  	= r.bed.rudece[1, 6],
		  		reads	    = sum(abs(r.reads)),
		  		lambda      = r.lambda,
		  		p.pois      = p.pois );
		r.df;}, mc.cores = ncores );
}

	if(length(DBIDs)>0)
	r.bed.list <- mclapply( 1:length(DBIDs), function(i) {
	    dbid <- as.character(DBIDs[i]);
    
		r.bed.idx1 <- which(gencode_transcript_ext$gene_id==dbid);
		r.bed.idx2 <- grep( paste(dbid, ".", sep=""), gencode_transcript_ext$gene_id );
		r.bed.idx <- unique( c(r.bed.idx1 ,r.bed.idx2) );

	    if(length(r.bed.idx)<1)
	    {
	    	cat ("! Can't find DBID=", dbid, ", missing data.\n");
	    	return(c(dbid=dbid, rep(NA, 8)));
	    }
    
    	r.bed <- gencode_transcript_ext[r.bed.idx, c(1,4,5,6,9,7), drop=F ];
	    colnames(r.bed) <- c('chr','start','end','id','score','strand');
		
		idx.bed.size1 <- which(r.bed$start >= r.bed$end);
		if( length(idx.bed.size1) > 0 )
		{
			temp.end <- r.bed$end[idx.bed.size1];
			r.bed$end[idx.bed.size1] <- r.bed$start[idx.bed.size1] + 1;
			r.bed$start[idx.bed.size1] <- temp.end;
		}
		
		# Query reads for each motif
		r.reads  <- try(bed6.region.bpQuery.bigWig( bw.plus, bw.minus, r.bed ));
		if( class(r.reads) == "try-error")
	    	return(c(dbid=dbid, rep(NA, 8)));

		bed.max <- which.max(abs(r.reads/((r.bed[,3] - r.bed[,2])+1)));
		
		p.pois <- ppois( abs(r.reads[bed.max]), r.lambda*abs(r.bed[bed.max,3] - r.bed[bed.max,2]) , lower.tail=F);

	    r.df   <- c( 
	    		dbid		= dbid,
	    		chr 		= r.bed[bed.max, 1],
		  		start   	= r.bed [bed.max,2],
		  	  	end     	= r.bed [bed.max,3],
		  	  	length     	= abs(r.bed [bed.max,3] - r.bed [bed.max,2] ),
		  		strand  	= r.bed [bed.max,6],
		  		reads	    = abs(r.reads[bed.max]),
		  		lambda      = r.lambda,
		  		p.pois      = p.pois );
		r.df;}, mc.cores = 1 );  # maybe bigwig cuase the crash if multi-cores are applied


	df.exp <- transform( do.call(rbind, r.bed.list) );
	
	df.idx <- c();
	for(i in 1:NROW(tfbs@extra_info))
	{
		idx <- which(as.character(df.exp$dbid)==as.character(tfbs@extra_info$DBID[i]));
		df.idx <-c(df.idx, idx[1]);
	}
	
	df.exp <- cbind( tfbs@extra_info$Motif_ID, df.exp[df.idx,,drop=F] );
	colnames(df.exp) <- c("Motif_ID", "DBID", "chr", "start", "end", "length", "strand", "reads", "lambda", "p.pois" );
	
	df.exp <- transform( df.exp, 
				start  = as.numeric(as.character(start) ),
				end    = as.numeric(as.character(end) ),
				length = as.numeric(as.character(length) ),
				reads  = as.numeric(as.character(reads) ),
				lambda = as.numeric(as.character(lambda) ),
				p.pois = as.numeric(as.character(p.pois) ) );
	
	# expressionlevel is data.frame including all information.
	tfbs@expressionlevel <- df.exp;

    try( unload.bigWig( bw.plus ) );
    try( unload.bigWig( bw.minus ) );
	
	return(tfbs);
}

#
# Extract the gencode information and save it into RDATA file
#
# e.g. 
# import_gencode( "human", "gencode.v22.annotation.gtf", "gencode.v22.rdata")

import_gencode <-function( species, file.gencode.gtf, ncores = 1)
{
	f.gtf <- try (file( file.gencode.gtf ) );
	if(class(f.gtf)=="try-error")
	{
		cat("! Failed to open GTF file(", file.gencode.gtf, ").\n");
		return(NULL);
	}
	
	
	bigdf <- read.table(f.gtf, header = F, sep="\t", skip=6);
	# sqldf package
	# bigdf <- sqldf("select * from f", dbname = tempfile(), 
	#			   file.format = list(header = F, row.names = F, sep="\t", skip=6));

	ids.multi <- strsplit(as.character(bigdf$V9), ";");

	cat(" ", length(ids.multi), "Items are loaded, it will take long time to parse...\n" );
	
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
	
	cat("  Parsing is finished.\n");

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
	
	return(gencode_transcript_ext);
}


tfbs.selectExpressed<-function( tfs, prob.sig=0.05, include.DBID.Missing=FALSE )
{
	if( include.DBID.Missing )
		tf.expresed <- which( tfs@expressionlevel$p.pois<=prob.sig | is.na(tfs@expressionlevel$p.pois) )
	else
		tf.expresed <- which( tfs@expressionlevel$p.pois<=prob.sig );
		
	if(length(tf.expresed)>0)
	{
		tfs@expressionlevel <- tfs@expressionlevel[ tf.expresed,,drop=F ];
		tfs@extra_info      <- tfs@extra_info[ tf.expresed,,drop=F ];
		tfs@mgisymbols      <- tfs@mgisymbols[ tf.expresed ];
		tfs@filename        <- tfs@filename[ tf.expresed ];
		tfs@pwm             <- tfs@pwm[ tf.expresed ];
		tfs@TFID            <- tfs@TFID[ tf.expresed ];

		cat("* After filtering by the gene expression,", length(tf.expresed), "expressed TFs are selected from", tfs@ntfs, "PWM files.\n"); 
		tfs@ntfs            <- length(tf.expresed);
	}
	
	return(tfs);
}
