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

	r.plus <- sum(abs(bigWig::bed6.region.bpQuery.bigWig( bw.plus, bw.minus, r.bed)));

	r.bed <- data.frame( seqnames=chromInfo[,1],
		starts=chromInfo[,2],
		ends=chromInfo[,3],
		names=".",
		scores=".",
		strands="-")

	r.minus <- sum(abs(bigWig::bed6.region.bpQuery.bigWig( bw.plus, bw.minus, r.bed)));

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

lambda_estimate_in_bam <- function( file.bam, file.twoBit, file.gencode.gtf, sample.size=100*1000, win.size=1000, pseudo_lambda=0.01, use.strand)
{
	import_gene_loci <-function(file.gencode.gtf)
	{
		pipe.cmd <- paste( "awk '{print $1,$2,$3,$4,$5 }' ", file.gencode.gtf, sep="")

		# for gzipped GTF file
		if(file_ext(file.gencode.gtf)=="gz" )
			pipe.cmd <- paste( "zcat ", ifelse( get_os()=="osx", " < ", " " ), file.gencode.gtf," | awk '{print $1,$2,$3,$4,$5 }' ", sep="");

		bigdf <- try( read.table( pipe(pipe.cmd), header = F ), silent = T);
		if (class(bigdf) =="try-error")
			check_command_error( bigdf, c("zcat", "awk") );

		if( exists("bigdf") && !is.null(bigdf) && !is.na(bigdf) )
		{
			colnames(bigdf) <- c("V1", "V2", "V3", "V4", "V5");

			#to make sure the ending postion > start position
			bigdf <- bigdf[ bigdf[,5] > bigdf[,4], ,drop=F];

			return(bigdf[,c(1,4,5 )]);
		}
		else
		{
			cat("! Failed to call awk command to get the reduced gencode file.\n");
			return(NULL);
		}
	}

	get_complement<-function( bed )
	{
		tmp.bed <- tempfile();
		write.table(bed, file=tmp.bed, quote=F, sep="\t", row.names=FALSE, col.names=FALSE,);

		tmp.genome <- tempfile();
		chromInfo <- get_chromosome_size(file.twoBit);
		chromInfo <- chromInfo[grep("_", chromInfo[,1], invert=TRUE),]
		chromInfo <- chromInfo[order( chromInfo[,1]),]
		write.table( chromInfo, file=tmp.genome, quote=F, sep="\t", row.names=FALSE, col.names=FALSE,);

		pipe.cmd <- paste("sort-bed ", tmp.bed, " | bedtools complement -i - -g ", tmp.genome );
		df.comp  <- try( read.table( pipe(pipe.cmd), header = F ), silent = T);
		if (class(df.comp) =="try-error")
			check_command_error( df.comp, c("sort-bed", "bedtools") );

		return(df.comp);
	}

	split_bed<-function( bed, winsize )
	{
		tmp.bed <- tempfile();
		write.table(bed, file=tmp.bed, quote=F, sep="\t", row.names=FALSE, col.names=FALSE,);

		pipe.cmd <- paste("bedtools makewindows -b ", tmp.bed , " -w ", winsize);
		df.split  <- try( read.table( pipe(pipe.cmd), header = F ), silent=T);
		if (class(df.split ) =="try-error")
			check_command_error( df.split, c("bedtools") );

		return(df.split);
	}

	Stats_mode <- function(x) {
		ux <- unique(x)
		ux[which.max(tabulate(match(x, ux)))]
	}

	df.bed <- import_gene_loci(file.gencode.gtf);
	if( is.null( df.bed ) ) return( NULL );

	df.comp <- get_complement( df.bed );
	df.comp[,2] <- df.comp[,2] + 1000;
	df.comp[,3] <- df.comp[,3] - 1000;
	df.comp.width <- df.comp[,3] - df.comp[,2];
	df.comp <- df.comp [ which( df.comp.width >= win.size), ];

	## split bed loci by 'bedtools' command
	df.comp.win <- split_bed( df.comp, win.size)
	df.comp.win <- df.comp.win [ which( df.comp.win[,3] - df.comp.win[,2] >= win.size), ];

	if ( sample.size >= NROW( df.comp.win ) )
		sample.size <- NROW( df.comp.win );

	sub.sample <- sample( 1:NROW(df.comp.win) )[ 1:sample.size ];
	df.comp.win <- df.comp.win[sub.sample, ];

	df.comp.strand <- rbind( cbind( df.comp.win, V4=".", V5=".", V6="+"),
	 				 		 cbind( df.comp.win, V4=".", V5=".", V6="-") );
	if(use.strand)
		reads <- get_bam_reads( file.bam, df.bed = df.comp.strand, use.strand=use.strand )
	else
		reads <- get_bam_reads( file.bam, df.bed = df.comp.win, use.strand=use.strand )

	## NO BAM file or BAM file is not indexed
	if( all(is.na(reads)))  return( NULL );

	## remove the ouliner reads
	if( quantile(reads, 0.99) > 0)
		reads <- reads[ which( reads <= quantile(reads, 0.99) ) ];

	# option 1: Mode
	lambda1 <- Stats_mode(reads);
	# option 2: Average
	lambda2 <- mean(reads);

	lambda <- max( lambda1, lambda2, pseudo_lambda);

	return( lambda );
}

get_bam_reads<-function(file.bam, df.bed = NULL, file.twoBit = NULL, use.strand = FALSE )
{
	options("scipen"=100, "digits"=4)

	read_bam <- function(df.bed, use.strand)
	{
		file.df3 <- tempfile();
		write.table( df.bed, file=file.df3, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t");
		pipe.cmd <- paste( "sort-bed ", file.df3, " | bedtools multicov -bams", file.bam, " -bed - ");
		if(use.strand)
			pipe.cmd <- paste( pipe.cmd, "-s");

		ts <- try(read.table( pipe(pipe.cmd) ), silent=T);
		unlink(file.df3);

		if(class(ts) == "try-error")
			check_command_error( ts, c("bedtools", "sort-bed") )
		else
			return(ts[, NCOL(ts) ] )
	}

	filename <- tempfile();
	## suppress samtools command, keep codes here
	if(is.null(df.bed) && is.null(file.twoBit))
	{
		## Mapped reads only(-F 4), Unmapped reads only(-f 4)
		pipe.cmd <- paste( "samtools view -c ", file.bam, " -F 4 ", sep=" ");
		ts <- try( read.table( pipe(pipe.cmd) ), silent=T );
		if(class(ts)== "try-error")
			check_command_error( ts, c("samtools") )
		else
			return(ts[1,1])

	}

	if(!is.null(df.bed))
	{
		if( !use.strand )
			return( read_bam(df.bed[, c(1:3)], FALSE) )
		else
			return( read_bam( df.bed[, c(1:6)], TRUE) );
	}
	else
	if( !is.null(file.twoBit) )
	{
		tb.gen <- get_chromosome_size( file.twoBit );
		df.bed <- data.frame(tb.gen[,1], 1, tb.gen[,2]);

		if( !use.strand )
			return( sum(read_bam(df.bed, FALSE) ) )
		else
		{
			reads.plus <- read_bam( data.frame(df.bed, ".", ".", "+"), TRUE);
			reads.minus <- read_bam( data.frame(df.bed, ".", ".", "-"), TRUE);
			return( sum(reads.plus + reads.minus) );
		}
	}
	else
		return(NA);
}


cpu.PRO.seq <- function(DBIDs, i.from, i.to, gencode_transcript_ext, file.bigwig.plus, file.bigwig.minus, reads.total, r.lambda, use.strand)
{
	r.df.len <- 11;
	bw.plus  <- try( bigWig::load.bigWig( file.bigwig.plus ) );
	bw.minus <- try( bigWig::load.bigWig( file.bigwig.minus ) );

	r.bed.list <- lapply( i.from:i.to, function(i)
	{
		dbid <- as.character(DBIDs[i]);

		r.bed.idx1 <- which(gencode_transcript_ext$gene_id==dbid);
		r.bed.idx2 <- grep( paste(dbid, ".", sep=""), gencode_transcript_ext$gene_id );
		r.bed.idx <- unique( c(r.bed.idx1 ,r.bed.idx2) );

		if(length(r.bed.idx)<1) return(c(dbid=dbid, rep(NA, r.df.len)));

		r.bed <- gencode_transcript_ext[r.bed.idx, c(1,4,5,12,6,7), drop=F ];
		colnames(r.bed) <- c('chr','start','end','id','score','strand');

		idx.bed.size1 <- which(r.bed$start >= r.bed$end);
		if( length(idx.bed.size1) > 0 )
		{
			temp.end <- r.bed$end[idx.bed.size1];
			r.bed$end[idx.bed.size1] <- r.bed$start[idx.bed.size1] + 1;
			r.bed$start[idx.bed.size1] <- temp.end;
		}

		# Query reads for each motif
		r.reads  <- try( bigWig::bed6.region.bpQuery.bigWig( bw.plus, bw.minus, r.bed ) )

		if( class(r.reads) == "try-error")
			r.df   <- c(dbid, rep(NA, r.df.len))
		else
		{
			bed.max <- which.max(abs( r.reads /((r.bed[,3] - r.bed[,2])+1)));
			p.pois  <- ppois( abs(r.reads[bed.max]), r.lambda * abs( r.bed[bed.max,3] - r.bed[bed.max,2] ) , lower.tail=F);
			length  <- abs(r.bed [bed.max,3] - r.bed [bed.max,2] );
			reads   <- abs(r.reads[bed.max]);
			RPKM    <- 10^9 * reads / length / reads.total;
			lambda0      <- r.lambda * abs( r.bed[bed.max,3] - r.bed[bed.max,2] );
			lambda.RPKM  <- 10^9 * r.lambda / reads.total;

			r.df   <- c(
					"DBID"      = dbid,
					"txID"      = as.character(r.bed [bed.max, 4]),
					"chr"       = as.character(r.bed [bed.max, 1]),
					"start"     = r.bed [bed.max, 2],
					"end"       = r.bed [bed.max, 3],
					"length"    = length,
					"strand"    = as.character(r.bed [bed.max,6]),
					"reads"	    = reads,
					"lambda"    = lambda0,
					"RPKM"      = RPKM,
					"lambda.RPKM"= lambda.RPKM,
					"p.pois"    = p.pois );
		}

		return(r.df);
	});

	if( !is.null(bw.plus) )	try( bigWig::unload.bigWig( bw.plus ) );
	if( !is.null(bw.minus) ) try( bigWig::unload.bigWig( bw.minus ) );

	return(r.bed.list);
}

cpu.RNA.seq <- function(DBIDs, i.from, i.to, gencode_transcript_ext, gencode_exon_ext, file.bam, reads.total, r.lambda, use.strand)
{
	r.df.len = 11;

	r.bed.list <- lapply( i.from:i.to, function(i) 	{
		dbid <- as.character(DBIDs[i]);

		r.bed.idx1 <- which(gencode_exon_ext$gene_id==dbid);
		r.bed.idx2 <- grep( paste(dbid, ".", sep=""), gencode_exon_ext$gene_id );
		r.bed.idx <- unique( c(r.bed.idx1 ,r.bed.idx2) );

		if(length(r.bed.idx)<1) return(c(dbid=dbid, rep(NA, r.df.len)));

		# chr|start|stop|id|score|strand
		r.bed.exon <- gencode_exon_ext[r.bed.idx, c(1,4,5,12,6,7), drop=F ];
		colnames(r.bed.exon) <- c('chr','start','end','id','score','strand');

		idx.bed.size1 <- which(r.bed.exon$start >= r.bed.exon$end);
		if( length(idx.bed.size1) > 0 )
		{
			temp.end <- r.bed.exon$end[idx.bed.size1];
			r.bed.exon$end[idx.bed.size1] <- r.bed.exon$start[idx.bed.size1] + 1;
			r.bed.exon$start[idx.bed.size1] <- temp.end;
		}

		r.reads  <- try( get_bam_reads( file.bam, df.bed=r.bed.exon, use.strand=use.strand ) );
		if( class(r.reads) == "try-error")
			r.df   <- c( dbid, rep(NA, r.df.len) )
		else
		{
			r.rna.idx1 <- which(gencode_transcript_ext$gene_id==dbid);
			r.rna.idx2 <- grep( paste(dbid, ".", sep=""), gencode_transcript_ext$gene_id );
			r.rna.idx  <- unique( c(r.rna.idx1 ,r.rna.idx2) );
			if(length(r.rna.idx)<1)
				r.df   <- c( dbid, rep(NA, r.df.len) )
			else
			{
				r.bed.rna <- gencode_transcript_ext[r.rna.idx, c(1,4,5,12,6,7), drop=F ];

				rna.reads  <- c();
				rna.length <- c();

				for(i in 1:length(r.rna.idx))
				{
					exon.idx   <- which(  r.bed.exon[,2]<r.bed.rna[i,3] &
										  r.bed.exon[,3]>r.bed.rna[i,2] &
										  as.character(r.bed.rna[i,6])==as.character(r.bed.exon[,6]) );

					if(length(exon.idx)>0)
					{
						rna.reads  <- c( rna.reads, sum(r.reads[exon.idx]));
						rna.length <- c( rna.length, sum( r.bed.exon[exon.idx,3]- r.bed.exon[exon.idx,2]) );
					}
					else
					{
						rna.reads  <- c( rna.reads, NA );
						rna.length <- c( rna.length, NA );
					}
				}

				max.idx <- which.max( rna.reads/rna.length );
				if(length(max.idx )<1)
					r.df   <- c(dbid, rep(NA, r.df.len))
				else
				{
					p.pois  <- ppois( abs(rna.reads[max.idx]), r.lambda * rna.length[max.idx], lower.tail=F);
					length  <- rna.length[max.idx];
					reads   <- rna.reads[max.idx];
					RPKM    <- 10^9 * reads / length / reads.total;
					lambda0      <- r.lambda * length ;
					lambda.RPKM  <- 10^9 * lambda0/ length / reads.total;

					r.df   <- c(
						"DBID"      = dbid,
						"txID"      = as.character(r.bed.rna [max.idx, 4]),
						"chr"       = as.character(r.bed.rna [max.idx, 1]),
						"start"     = r.bed.rna [max.idx, 2],
						"end"       = r.bed.rna [max.idx, 3],
						"length"    = length,
						"strand"    = if(use.strand) as.character(r.bed.rna [max.idx,6]) else ".",
						"reads"	    = reads,
						"lambda"    = lambda0,
						"RPKM"      = RPKM,
						"lambda.RPKM"= lambda.RPKM,
						"p.pois"    = p.pois );

				}
			}
		}

		return(r.df);
	});

	return(r.bed.list);
}

#' Gets expression level of target TF.
#' USE tf_info$DBID to find gene information encoded by GENCODE V21
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

tfbs_getExpression <- function(tfbs,
								file.twoBit,
								file.gencode.gtf,
								file.bigwig.plus=NA,
								file.bigwig.minus=NA,
								file.bam=NA,
								seq.datatype=c("GRO-seq", "PRO-seq", "RNA-seq" ),
								use.strand = FALSE,
								ncores = 1 )
{
	stopifnot( NROW(tfbs@tf_info)>0 );

	if( !missing(seq.datatype))
		seq.datatype <- match.arg( seq.datatype );

	if( ( missing(seq.datatype) || is.na(seq.datatype) ) && !missing(file.bam) )
		seq.datatype <- "RNA-seq";

	if( missing(seq.datatype) || is.na(seq.datatype))
		seq.datatype <- "GRO-seq";

	if( seq.datatype=="RNA-seq" )
	{
		if ( (is.na(file.bam) || missing(file.bam) ))
			stop("Not specify the indexed BAM file for RNA-seq data.");
		if ( !file.exists(file.bam) )
			stop(paste("BAM file is not accessible, File=", file.bam));
	}

	if( seq.datatype!="RNA-seq" )
	{
		if( is.na(file.bigwig.plus) || is.na(file.bigwig.minus) )
			stop("Not specify the plus and minus Bigwig files for GRO-seq or PRO-seq data.");

		if( !file.exists (file.bigwig.plus)  )
			stop(paste("Bigwig file is not accessible, File=", file.bigwig.plus));

		if( !file.exists (file.bigwig.minus)  )
			stop(paste("Bigwig file is not accessible, File=", file.bigwig.minus));
	}

	if( !file.exists (file.twoBit)  )
		stop(paste("Twobit file is not accessible, File=", file.twoBit));

	if( !file.exists (file.gencode.gtf)  )
		stop(paste("Gencode file is not accessible, File=", file.gencode.gtf));

	reads.total  <- 0;
	reads.lambda.kb <- 0;
	win.size <- 2000;
	pseudo_lambda <- 0.02;

	# load gencode RDATA file and set the table of gencode_transcript_ext
	gencode_transcript_ext <- try( import_gencode( tfbs@species, file.gencode.gtf, "transcript") );
	if(is.null(gencode_transcript_ext) || class(gencode_transcript_ext)=="try-error")
		stop("Gencode data can not be found in the GTF file specified by the parameter of file.gencode.gtf.")
	else
		cat(" ", NROW(gencode_transcript_ext), "transcripts are selected from GENCODE dataset for", seq.datatype, ".\n");

	if(seq.datatype=="GRO-seq" ||seq.datatype=="PRO-seq" )
	{
		if(!requireNamespace("bigWig", quietly = TRUE))
			stop("Package bigWig is required to calculate reads for GRO-seq or PRO-seq data");

		# Load bigWig files(minus and plus)
		bw.plus  <- try( bigWig::load.bigWig( file.bigwig.plus ) );
		bw.minus <- try( bigWig::load.bigWig( file.bigwig.minus ) );

		if( class(bw.plus)=="try-error" || class(bw.minus)=="try-error" )
			stop("Failed to load bigwig files.");

		reads.total <- sum(abs(c(bw.plus$primaryDataSize, bw.minus$primaryDataSize)) );

		if( !is.null(bw.plus) )	try( bigWig::unload.bigWig( bw.plus ) );
		if( !is.null(bw.minus) ) try( bigWig::unload.bigWig( bw.minus ) );

		cat(" ", reads.total, "Reads in", file.bigwig.plus, "and", file.bigwig.minus,"\n");
	}
	else if(seq.datatype=="RNA-seq")
	{
		# load gencode GTF file and set the table of gencode_exon_ext
		gencode_exon_ext <- try( import_gencode( tfbs@species, file.gencode.gtf, "exon") );
		if(is.null(gencode_exon_ext) || class(gencode_exon_ext)=="try-error")
			stop("Gencode data can not be found in the GTF file specified by the parameter of file.gencode.gtf.")
		else
			cat(" ", NROW(gencode_exon_ext), "exons are selected from GENCODE dataset.\n");

		reads.total <- get_bam_reads( file.bam, file.twoBit = file.twoBit, use.strand=use.strand);
		cat(" ", reads.total, "mapped reads in", file.bam, "\n");

		reads.lambda.kb <- lambda_estimate_in_bam( file.bam, file.twoBit, file.gencode.gtf, win.size=win.size, pseudo_lambda=pseudo_lambda, use.strand=use.strand );
		if( is.null(reads.lambda.kb) )
			stop("Failed to load the BAM file.");

		if( reads.lambda.kb < pseudo_lambda )
			stop("Lambda of Poisson distribution is too samll( = 0 reads/kb).");

		cat("*", "Lambda of Poisson distribution is estimated (=", round(reads.lambda.kb*1000/win.size, 2), "reads/kb).\n");
	}
	else
		stop("The function tfbs.getExpression only supports GRO-seq, PRO-seq an RNA-seq data.");

	if (seq.datatype=="RNA-seq")
		r.lambda <- reads.lambda.kb/win.size
	else
		r.lambda <- 0.04 * reads.total/10751533/1000 ;


	DBIDs <- unique( as.character(tfbs@tf_info$DBID) );
	if( length(DBIDs)==0 )
	{
		cat("!  No DBID found in the Gencode file.\n");
		return(tfbs);
	}

	if( ncores > length(DBIDs) ) ncores <- length(DBIDs);

	df.exp <- mclapply( 1:ncores, function(i){
		sect <- round( seq( 1, length(DBIDs)+1, length.out = ncores+1) );

		r.bed.list <- c();
		if(seq.datatype=="RNA-seq")
			r.bed.list <- cpu.RNA.seq( DBIDs, sect[i], sect[i+1]-1, gencode_transcript_ext, gencode_exon_ext, file.bam, reads.total, r.lambda, use.strand)
		else
			r.bed.list <- cpu.PRO.seq( DBIDs, sect[i], sect[i+1]-1, gencode_transcript_ext, file.bigwig.plus, file.bigwig.minus, reads.total, r.lambda, use.strand );

		df.exp <- transform( do.call(rbind, r.bed.list) );
		colnames(df.exp) <- c("DBID", "txID", "chr", "txStart", "txEnd", "txLength", "strand", "reads","lambda", "RPKM", "lambda.RPKM", "p.pois");

		return( df.exp );
	}, mc.cores = ncores );


	df.exp0 <- do.call( rbind, df.exp );
	df.idx  <- match( as.character(tfbs@tf_info$DBID), as.character(df.exp0$DBID) )
	df.exp  <- cbind( tfbs@tf_info$Motif_ID, df.exp0[df.idx,,drop=F] );
	colnames(df.exp) <- c("Motif_ID", "DBID", "txID", "chr", "txStart", "txEnd", "txLength", "strand", "reads", "lambda", "RPKM", "lambda.RPKM", "p.pois" );

	# these dummy statements are setup here to pass R CMD check rtfbsdb --as-cran
	txStart  <-NA;
	txEnd    <-NA;
	txLength <-NA;
	reads    <-NA;
	lambda   <- NA;
	p.pois   <-NA;
	RPKM     <- NA;
	lambda.RPKM <- NA;

	df.exp <- transform( df.exp,
				"txStart"  = as.numeric(as.character(txStart) ),
				"txEnd"    = as.numeric(as.character(txEnd) ),
				"txLength" = as.numeric(as.character(txLength) ),
				"reads"    = as.numeric(as.character(reads) ),
				"lambda"   = as.numeric(as.character(lambda) ),
				"RPKM"     = as.numeric(as.character(RPKM) ),
				"lambda.RPKM" = as.numeric(as.character(lambda.RPKM) ),
				"p.pois"   = as.numeric(as.character(p.pois) ) );

	## for RNA-seq, txLength is sum of exon length, not the length of transcipt
	if(seq.datatype=="RNA-seq")
		colnames(df.exp)[7] <- "exonLen";

	if(NROW(df.exp)>0)  rownames(df.exp) <- c(1:NROW(df.exp));

	# expressionlevel is data.frame including all information.
	tfbs@expressionlevel <- df.exp;

	DBD.missing <- length(which(is.na(df.exp$reads)));
	if( DBD.missing>0 )
		cat( "*", DBD.missing, "motifs did not find DBID in the Gencode file.\n");

	return(tfbs);
}

#
# Extract the gencode information and save it into RDATA file
#
# e.g.
# import_gencode( "human", "gencode.v22.annotation.gtf", V3.type="GRO-seq")

import_gencode <-function( species, file.gencode.gtf, V3.type=NA )
{
	if( missing(V3.type)) V3.type <- "transcript";

	awk.cmd <- paste( "awk '($3==\"", V3.type, "\"){gsub( /\\\";?/, \"\", $10);gsub( /\\\";?/, \"\", $12);print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ", file.gencode.gtf, sep="");

	# for gzipped GTF file
	if(file_ext(file.gencode.gtf)=="gz" )
		awk.cmd <- paste( "zcat ", ifelse( get_os()=="osx", " < ", " " ),  file.gencode.gtf," | awk '($3==\"", V3.type, "\"){gsub( /\\\";?/, \"\", $10);gsub( /\\\";?/, \"\", $12);print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ", sep="");

	bigdf <- try( read.table( pipe(awk.cmd), header = F ), silent=T);
	if (class(bigdf) =="try-error")
		check_command_error( bigdf, c("zcat", "awk") );

	if( exists("bigdf") && !is.null(bigdf) && !is.na(bigdf) )
	{
		colnames(bigdf) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "gene_id_label", "gene_id", "tx_id_label", "tx_id" );
		return(bigdf);
	}
	else
	{
		cat("! Failed to call awk command to get the reduced gencode file.\n");
		return(NULL);
	}
}

# make_gencode_rdata( "human", "gencode.v22.annotation.gtf", "gencode.v22.rdata")

make_gencode_rdata <-function( species, file.gencode.gtf, ncores = 1)
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

tfbs_selectExpressedMotifs <- function( tfbs,
							file.twoBit,
							file.gencode.gtf,
							file.bigwig.plus = NA,
							file.bigwig.minus = NA,
							file.bam = NA,
							seq.datatype=c("GRO-seq", "PRO-seq", "RNA-seq" ),
							pvalue.threshold = 0.05,
							lowest.reads.RPKM = NA,
							include.DBID.missing = TRUE,
							use.strand = FALSE,
							ncores = 1)
{
	selectExpressed<-function( tfs, prob.sig=0.05, include.DBID.Missing=FALSE )
	{
		if( include.DBID.Missing )
			tf.expresed <- which( tfs@expressionlevel$p.pois<=prob.sig | is.na(tfs@expressionlevel$p.pois) )
		else
			tf.expresed <- which( tfs@expressionlevel$p.pois<=prob.sig );

		if (!is.na(lowest.reads.RPKM) && length(tf.expresed)>0 )
		{
			rpkm <- tfs@expressionlevel[tf.expresed, 'RPKM'];
			if( include.DBID.Missing )
				tf.expresed <- tf.expresed[ which(rpkm >= lowest.reads.RPKM | is.na(rpkm) ) ]
			else
				tf.expresed <- tf.expresed[ which(rpkm >= lowest.reads.RPKM) ]
		}

		if(length(tf.expresed)>0)
		{
			tfs@expressionlevel <- tfs@expressionlevel[ tf.expresed,,drop=F ];
			tfs@tf_info         <- tfs@tf_info[ tf.expresed,,drop=F ];
			tfs@mgisymbols      <- tfs@mgisymbols[ tf.expresed ];
			tfs@filename        <- tfs@filename[ tf.expresed ];
			tfs@pwm             <- tfs@pwm[ tf.expresed ];
			tfs@distancematrix  <- matrix(, nrow=0, ncol=0);
			tfs@cluster         <- matrix(, nrow=0, ncol=0);

			cat("*", length(tf.expresed),
				"expressed TFs are selected from", tfs@ntfs, "motifs after filtering by the gene expression.\n");

			tfs@ntfs            <- length(tf.expresed);
		}

		return(tfs);
	}

	tfs <- tfbs_getExpression( tfbs,
								file.twoBit,
								file.gencode.gtf,
								file.bigwig.plus,
								file.bigwig.minus,
								file.bam = file.bam,
								seq.datatype = seq.datatype,
								use.strand = use.strand,
								ncores = ncores  );

	if( NROW( tfs@expressionlevel ) > 0 )
		tfs <- selectExpressed( tfs, pvalue.threshold, include.DBID.missing )
	else
		cat("! Failed to calculate the gene expression and select expressed TFs.\n");

	return(tfs);
}