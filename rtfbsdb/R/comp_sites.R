##
## Comparative TFBS enrichment
##
##

#' Create GC bins that span the joint set GC content distribution.
#'
#' @param set1.ms first set of sequences
#' @param set2.ms second set of sequences
#' @param n integer number of bins (default 4)
#' @return list with two integer vectors, bins1 and bins2, indicating bin assignments for each set.

joint.gc.quantile.bins <- function(set1.ms, set2.ms, n = 4) 
{
	gc1 = gcContent.ms(set1.ms)
	gc2 = gcContent.ms(set2.ms)

	cutoffs = quantile(c(gc1, gc2), seq(from = 0, to = 1, length.out = n + 1), names = FALSE, na.rm=T)
	bin1 <- findInterval(gc1, cutoffs, all.inside = TRUE)
	bin2 <- findInterval(gc2, cutoffs, all.inside = TRUE)

	# check each bin has at least one sequence
	n.seqs.per.bin = sapply(1:n, function(i) min(sum(bin1 == i, na.rm=T), sum(bin2 == i, na.rm=T)))

	stopifnot(all(n.seqs.per.bin > 0))

	return(list(bins1 = bin1, bins2 = bin2))
}

# NOTE: empirical p-values are computed under the assumption that the negative set
#       is not bound by the specified factor (but is accessible)
comparative_scan_rtfbs.ms <- function(pwm, 
						positive.ms, 
						negative.ms, 
						background.ms, 
						background.mm, 
						threshold = NA, 
						threshold.type = NA, 
						calc.empirical.pvalue = FALSE ) 
{
	if (threshold.type == "fdr") 
	{
		pos.sites = score.ms(positive.ms, pwm, background.mm, threshold = 0)
		neg.sites = score.ms(negative.ms, pwm, background.mm, threshold = 0)
		bg.sites  = score.ms(background.ms, pwm, background.mm, threshold = 0)
		fdrtbl    = calc.fdr(concat.ms(positive.ms, negative.ms), rbind(pos.sites, neg.sites), background.ms, bg.sites)

		## if something wrong
		if (is.null(fdrtbl))
			return(NULL)

		# find first score threshold below FDR threshold
		idxs = which(fdrtbl$fdr <= threshold)
		if (length(idxs) == 0)
			return(NULL)

		thresh = fdrtbl[idxs[length(idxs)], 1];
		pos.mask = pos.sites$score >= thresh;
		neg.mask = neg.sites$score >= thresh;

		Npos = length(unique(pos.sites$seqname[pos.mask]));
		Nneg = length(unique(neg.sites$seqname[neg.mask]));

		# NOTE: empirical p-values are computed under the assumption that the negative set
		#       is not bound by the specified factor (but is accessible)
		if (calc.empirical.pvalue) {
			epvals = sapply(pos.sites$score[pos.mask], function(thresh) {
				K = sum(neg.sites$score >= thresh)
				Z = length(neg.sites$score)

				(K + 1)/(Z + 1) # empirical p-value
			})

			return(list(Nneg   = Nneg, 
						Npos   = Npos, 
						thresh = thresh, 
						sites  = pos.sites[pos.mask, ], 
						empirical.pvalues = epvals))
		}

		return(list(Nneg = Nneg, Npos = Npos, thresh = thresh, sites = pos.sites[pos.mask, ]))
	} 
	# threshold.type == score
	else 
	{
		pos.sites = score.ms(positive.ms, pwm, background.mm, threshold = threshold)
		neg.sites = score.ms(negative.ms, pwm, background.mm, threshold = threshold)

		Npos = length(unique(pos.sites$seqname))
		Nneg = length(unique(neg.sites$seqname))

		if (calc.empirical.pvalue) {
			epvals = sapply(pos.sites$score, function(thresh) {
				K = sum(neg.sites$score >= thresh)
				Z = length(neg.sites$score)

				(K + 1)/(Z + 1) # empirical p-value
			})

			return(list(Nneg   = Nneg, 
						Npos   = Npos, 
						thresh = threshold, 
						sites  = pos.sites, 
						empirical.pvalues = epvals))
		}

		return(list(Nneg = Nneg, Npos = Npos, thresh = threshold, sites = pos.sites))
	}
}

# NOTE: empirical p-values are computed under the assumption that the negative set
#       is not bound by the specified factor (but is accessible)
comparative_scan_rtfbs <- function(pwm, file.twoBit, 
							positive.bed, 
							negative.bed, 
							threshold = NA, 
							threshold.type = NA, 
							gc.groups = 1, 
							background.order = 2, 
							background.length = 100000, 
							calc.empirical.pvalue = FALSE) 
{
	# read sequences
	positive.ms = read.seqfile.from.bed(positive.bed, file.twoBit)
	negative.ms = read.seqfile.from.bed(negative.bed, file.twoBit)

	# compute GC quantile masks
	gcBins = joint.gc.quantile.bins(positive.ms, negative.ms, n=gc.groups)

	# get per GC quantile results
	result.sites = NULL
	Npos   = 0
	Nneg   = 0
	thresh = NULL
	epvals = NULL

	for (idx in 1:gc.groups) 
	{
		# compute background model
		pos.ms.i = positive.ms[gcBins$bins1 == idx]
		neg.ms.i = negative.ms[gcBins$bins2 == idx]
		both.ms  = concat.ms(pos.ms.i, neg.ms.i)
		bg.mm    = build.mm(both.ms, background.order)

		# generate background sequences
		bg.ms <- simulate.ms(bg.mm, background.length)

		# collect sequences
		res = comparative_scan_rtfbs.ms(pwm, 
									pos.ms.i, 
									neg.ms.i, 
									bg.ms, 
									bg.mm, 
									threshold = threshold, 
									threshold.type = threshold.type, 
									calc.empirical.pvalue = calc.empirical.pvalue)

		if (!is.null(res)) {
			result.sites = rbind(result.sites, res$sites)
			Npos   = Npos + res$Npos
			Nneg   = Nneg + res$Nneg
			thresh = c(thresh, res$thresh)
	
			if (calc.empirical.pvalue)
				epvals = c(epvals, res$empirical.pvalues)
		}
	}
	
	# fisher test
	pos.left <- dim(positive.bed)[1] - Npos;
	if( pos.left<0 ) pos.left <- 0;
	neg.left <- dim(negative.bed)[1] - Nneg;
	if( neg.left<0 ) neg.left <- 0;
	
	tbl = rbind(c(Npos, Nneg), c( pos.left, neg.left ));
	pval = fisher.test(tbl)$p.value;
	
	# merge results
	bed = tfbs_to_bed(result.sites, "pwm")
	if (!is.null(bed) && dim(bed)[1] > 0) {
		ord = order(bed[,1], bed[,2])
		bed = bed[ord,]

		if (calc.empirical.pvalue)
			epvals = epvals[ord]
	}

	if (calc.empirical.pvalue)
		result = list(Npos     = Npos, 
					  Nneg     = Nneg, 
					  pvalue   = pval, 
					  thresh   = thresh, 
					  sites    = bed, 
					  empirical.pvalues = epvals)
	else
		result = list(Npos   = Npos, 
					  Nneg   = Nneg, 
					  pvalue = pval, 
					  thresh = thresh, 
					  sites  = bed)
}

tfbs_to_bed <- function(sites, tf.name) 
{
	## Alternative formulation fails if sites is empty.
	if (NROW(sites) == 0)
		return(NULL)
	
	spl <- strsplit(as.character(sites$seqname), ":|-")

	chroms <- sapply(spl, function(pair) pair[1]);
	starts <- as.integer(sapply(spl, function(pair) pair[2]));
	
	bed = data.frame(
			chrom      = chroms,
			chromStart = as.integer(starts + sites$start - 1), # add site offset
			chromEnd   = as.integer(starts + sites$end),
			name       = tf.name,
			score      = sites$score,
			strand     = sites$strand )
	
	return(bed)
}

background.check<-function( gc.pos, gc.neg, gc.correction, file.pdf.vioplot=NA, verbose=TRUE )
{
	pdf.output <- FALSE;
	gc.test <- wilcox.test(gc.pos, gc.neg, conf.int=TRUE, conf.level=0.9 );
	
	# Actually 0.01 is not good for wilcox.test, 
	if( gc.test$p.value<0.01 && verbose)
	{
		cat("! Difference between GC content in negative and positive TREs (use 'gc.correction.pdf' to see pdf figure):\n"); 
		cat("  p-value (Wilcoxon-Mann-Whitney test) =", gc.test$p.value, "\n"); 
		cat("  median/sample size of GC positive =", median(gc.pos), "/", length(gc.pos), "\n"); 
		cat("  median/sample size of GC negative =", median(gc.neg), "/", length(gc.neg), "\n");
		
		if( !is.null(file.pdf.vioplot) && !is.na( file.pdf.vioplot)  )
		{
			r.try <- try ( pdf(file.pdf.vioplot) );
			if( class(r.try) != "try-error" )
			{
				vioplot( gc.pos, gc.neg, names=c("Positive", "Negative") );
				abline(h=median(gc.pos), lty="dotted");
				pdf.output <- TRUE;
			}
			else
				cat("  Failed to output the vioplot figure for the results.\n");
		}		
	}
	
	# No need to do sampling the background data.
	if( gc.test$p.value > 0.01 || gc.correction==FALSE )
	{
		# PDF has not been finished until dev.off(); 
		if( pdf.output )
		{
			dev.off();
			cat("* Please check the vioplot figure to make sure, the vioplot figure: ", file.pdf.vioplot, "\n" );
		}
		return(NULL);
	}
	
	if( length(gc.neg)< 5000 )
	{
		if(verbose) cat("! Failed to make GC correction due to small bed data(size<5000).\n");
		return(NULL);
	}
	
	resample <- function(ref, orig, n=10000, nbins=10) {
		## (1) Break gcContent down into 10 equally sized bins.
		breaks <- seq(0, 1, length.out=nbins) #seq(min(ref), max(ref), length.out=nbins)

		## (2) Get the empirical frequencies of each bin.
		empir <- sapply(1:(NROW(breaks)-1), function(x) {sum(breaks[x]<= ref & ref < breaks[x+1])/NROW(ref)})
		
		## (3) Re-sample TREs in the BG set w/ probability proportional to the bin.
		resamp_prob <- rep(1, NROW(orig))
		for(i in 1:NROW(empir)) {
			incl <- breaks[i] <= orig & orig < breaks[i+1]
			resamp_prob[incl] <- empir[i]/ sum(incl)
		}

		sample(1:NROW(resamp_prob), n, prob= resamp_prob, replace=FALSE)
	}
	
	try.sample<-function(n.sample)
	{
		indx.bgnew <- resample( gc.pos, gc.neg, n=n.sample)
		gc.testx <- wilcox.test(gc.pos, gc.neg[indx.bgnew], conf.int=TRUE, conf.level=0.9 );
		return(gc.testx$p.value);
	}

	## Make n.sample as large as possible, up to 10k sequences.  More BG gains statistical power.
	#n.sample <- min(c(10000, round(length(gc.neg)/10)))

	# In order to get the high p-value, multiple resample sizes are used to do test, find a best one.
	# Need to consider ? 
	
	#ns.sample <- c ( round(length(gc.neg)/c(5,10,15,20,25)), 1000);
	#ns.sample <- ns.sample[ ns.sample>=1000 ];

	ns.sample <- c ( round(length(gc.neg)/c( 2:10,15,20,25)), 1000, 800, 500 );
	ns.sample <- ns.sample[ ns.sample >= 500 ];
	
	ns.pvalue <- unlist(lapply(ns.sample, try.sample));
	n.sample  <- ns.sample[ which.max(ns.pvalue) ];

	indx.bgnew <- resample( gc.pos, gc.neg, n=n.sample)
	gc.test2 <- wilcox.test(gc.pos, gc.neg[indx.bgnew], conf.int=TRUE, conf.level=0.9 );
	
	if(verbose)
	{
		cat("* After the resampling from negative TREs:\n"); 
		cat("  p-value (Wilcoxon-Mann-Whitney test) =", gc.test2$p.value, "\n"); 
		cat("  median/sample size of GC negative =", median(gc.neg[indx.bgnew]), "/", length(indx.bgnew), "\n");
	}
	
	if( pdf.output )
	{
		vioplot(gc.pos, gc.neg, gc.neg[indx.bgnew], names=c("Positive", "Negative", "Negative.resample")); 
		abline(h=median(gc.pos), lty="dotted")
		dev.off();
		cat("* Please check the vioplot figure to make sure, the vioplot figure: ", file.pdf.vioplot, "\n" );
	}

	## return sampling background.
	return( indx.bgnew );
}

background.generate <- function( positive.bed )
{
	pos.range  <- range(positive.bed[,3] - positive.bed[,2]);
	
	chr.size    <- aggregate(positive.bed[,3], list(positive.bed[,1]), max);

	file.chr.size <- tempfile();
	write.table(chr.size, file=file.chr.size, quote=F, row.names=F, col.names=F, sep="\t");	

	file.pos.bed <- tempfile();
	write.table( positive.bed, file=file.pos.bed, quote=F, row.names=F, col.names=F, sep="\t");	
	
	cmd.pipe <- paste("sort-bed ", file.pos.bed, " | bedtools complement -i - -g ", file.chr.size, sep=" ");
	bed.complement <- read.table( pipe( cmd.pipe ), header=F );

	bed.complement <- bed.complement[ -which( bed.complement[,3] - bed.complement[,2] <= pos.range[2] ),]
	bed.width      <- bed.complement[,3] - bed.complement[,2];

	n.neg <- sample( 1:length(bed.width), max(NROW(positive.bed)*2, 50000), replace=TRUE, 
					prob = (as.numeric(bed.width-pos.range[2]))/sum(as.numeric(bed.width-pos.range[2]) ));
	
	n.start <- unlist(lapply( n.neg, function(i){ bed.complement[i,2] + round( runif(1, 1, bed.width[i])) } ) )
	n.stop <-  unlist(lapply( 1:length(n.start), function(i){ n.start[i] + round( runif(1, pos.range[1], pos.range[2])) } ))
	
	unlink( file.chr.size );
	unlink( file.pos.bed );
	
	return(data.frame(chr=bed.complement[n.neg,1], start=n.start, stop=n.stop));
}

comparative_scanDb_rtfbs <- function( tfbs, file.twoBit, 
					positive.bed, 
					negative.bed, 
					file.prefix = NA, 
					cluster.mat = NA, 
					ncores = 1, 
					gc.correction = TRUE, 
					gc.correction.pdf = NA, 
					gc.correction.verbose = TRUE, 
					threshold = NA, 
					threshold.type = NA, 
					gc.groups=1, 
					background.order = 2, 
					background.length = 100000, 
					pv.adj = NA ) 
{
	stopifnot(class(tfbs) == "tfbs")
	
	if( missing(pv.adj) ) pv.adj <- "bonferroni";

	if( !is.na(file.prefix))
		if( !check_folder_writable( file.prefix ) ) 
			stop(paste("Can not create files starting with the prefix:", file.prefix));

	if(NROW(negative.bed) < 50)
		stop(paste("Negative BED file contains only ", NROW(negative.bed), " entries.  Strongly suggest size is increased."))

	# read sequences
	positive.ms = read.seqfile.from.bed(positive.bed, file.twoBit)
	negative.ms = read.seqfile.from.bed(negative.bed, file.twoBit)

	gc.pos <- gcContent.ms(positive.ms);
	gc.neg <- gcContent.ms(negative.ms);
	
	# for unknown reason, some ms objects get NaN values for gc content.
	# here we need to remove these data.
	if( length(which(is.na(gc.pos)))>0 ) 
	{
		positive.ms <- positive.ms [ -which(is.na(gc.pos)) ];
		gc.pos <- gc.pos [ -which(is.na(gc.pos)) ];
	}	

	if( length(which(is.na(gc.neg)))>0 )
	{
		negative.ms <- negative.ms [ -which(is.na(gc.neg)) ];
		gc.neg <- gc.neg [ -which(is.na(gc.neg)) ];
	}
	
	# detect the difference of gcContent between positive and negative TREs, 
	# if the difference is significant, make a correction for the negative TREs 
	# based on the resampling method.	
	
	bg.sample <- background.check( gc.pos, 
								 gc.neg, 
								 gc.correction, 
								 gc.correction.pdf, 
								 verbose=gc.correction.verbose )
	if(!is.null(bg.sample))
	{
		gc.neg <- gc.neg[ bg.sample ];
		negative.ms  <- negative.ms [ bg.sample ];
	}	

	# compute CG quantile masks
	gcBins = joint.gc.quantile.bins(positive.ms, negative.ms, n=gc.groups)
	
	# partition sequences
	seqs.subsets = lapply(1:gc.groups, function(idx) {
		pos.ms.i = positive.ms[gcBins$bins1 == idx]
		neg.ms.i = negative.ms[gcBins$bins2 == idx]
	
		list(pos = pos.ms.i, neg = neg.ms.i)
	})
	
	# compute background models
	gcBacks = lapply(1:gc.groups, function(idx) {
		both.ms = concat.ms(seqs.subsets[[idx]]$pos, seqs.subsets[[idx]]$neg)
		build.mm(both.ms, background.order)
	})

	# compute background sequences
	gcBacks.ms = lapply(1:gc.groups, function(idx) {
		simulate.ms(gcBacks[[idx]], background.length)
	})

	if(missing(cluster.mat)) cluster.mat <- cbind(1:tfbs@ntfs, 1); 

	# iterate over TF set
	binding_all <- mclapply(cluster.mat[,1], function(i, ...) {
		# get PWM information
		pwm = tfbs@pwm[[i]];

		pwm.name = paste(tfbs@mgisymbols[i], "@", i, sep='');

		# scan sequences, per GC bin
		result.sites = NULL;
		Npos = 0;
		Nneg = 0;

		for (idx in 1:gc.groups) {
			pos.ms.i = seqs.subsets[[idx]]$pos;
			neg.ms.i = seqs.subsets[[idx]]$neg;
			bg.i.mm  = gcBacks[[idx]];
			bg.i.ms  = gcBacks.ms[[idx]];

			res = comparative_scan_rtfbs.ms( pwm, 
											 pos.ms.i, 
											 neg.ms.i, 
											 bg.i.ms, 
											 bg.i.mm, 
											 threshold = threshold, 
											 threshold.type = threshold.type );

			if (!is.null(res)) {
				result.sites = rbind(result.sites, res$sites);
				Npos = Npos + res$Npos;
				Nneg = Nneg + res$Nneg;
			}
		}
	
		# process result
		result.bed = tfbs_to_bed(result.sites, pwm.name)
	
		pos.left <- length(gc.pos) - Npos;
		if( pos.left<0 ) pos.left <- 0;
		neg.left <- length(gc.neg) - Nneg;
		if( neg.left<0 ) neg.left <- 0;

		tbl  = rbind(c(Npos, Nneg), c(pos.left, neg.left))
		pval = fisher.test(tbl)$p.value;
	
		# save sites
		starch.file <- NA;
		if( !is.na(file.prefix) )
		{
			starch.file = paste(file.prefix, ".", i, ".bed.starch", sep='')
			write.starchbed(result.bed, starch.file)
		}
	
		Npos.tmp <- Npos;
		if(Npos.tmp==0) Npos.tmp <- 1;
		Nneg.tmp <- Nneg;
		if(Nneg.tmp==0) Nneg.tmp <- 1;
	
		Nneg.expected <- Nneg * length(gc.pos)/length(gc.neg); 
		fe.ratio <- (Npos.tmp/length(gc.pos) )/(Nneg.tmp/length(gc.neg));
	
		# return info
		return(data.frame('motif.id' = "", 
						  'tf.name'  = pwm.name, 
						  'Npos'     = Npos, 
						  'Nneg'     = Nneg, 
						  'expected' = round( Nneg.expected, 1), 
						  'fe.ratio' = fe.ratio, 
						  'pvalue'   = pval, 
						  'pv.adj'   = pval, 
						  'starch'   = starch.file,
						  'gc.pos'   = length(gc.pos),
						  'gc.neg'   = length(gc.neg) ) );
	}, mc.cores = ncores)

	# recombine results
	r.df <- do.call("rbind", binding_all);

	#put this dummy variable here to pass R CMD check rtfbsdb --as-cran	
	starch <- NULL;

	if (NROW(r.df)>0)	
	{
		r.df$pv.adj <- adjust.pvale( r.df$pvalue, cluster.mat, pv.adj );

		# the following 'starch' is a column in r.df, not equal NULL;
		if( missing(file.prefix) || is.na(file.prefix) )  
			r.df <- subset( r.df, select = -starch);
	}

	return( list(ret=r.df, bg.sample = bg.sample ) );
}

tfbs_enrichmentTest<-function( tfbs, file.twoBit, 
					positive.bed, 
					negative.bed = NA, 
					file.prefix = NA, 
					use.cluster = FALSE, 
					ncores = 1,
					gc.correction = TRUE, 
					gc.correction.pdf = NA, 
					gc.robust.rep = NA,
					threshold = 6, 
					threshold.type = c("score", "fdr"), 
					gc.groups = 1, 
					background.order = 2, 
					background.length = 100000, 
					pv.adj = p.adjust.methods) 
{
	stopifnot(class(tfbs) == "tfbs")

	if( !check_bed( positive.bed ) ) 
		stop("Wrong format in the parameter of 'positive.bed', at least three columns including chromosome, strat, stop.");
	if( !missing(negative.bed) && !check_bed( negative.bed ) ) 
		stop("Wrong format in the parameter of 'negative.bed', at least three columns including chromosome, strat, stop.");

	if( missing( threshold.type ) ) threshold.type <- "score";
	if( threshold.type == "score" && missing( threshold ) ) threshold <- 6;
	if( threshold.type == "fdr" && missing( threshold ) ) threshold <- 0.1;

	if( !missing(pv.adj)) pv.adj <- match.arg(pv.adj)
	if( missing(pv.adj) ) pv.adj <- "bonferroni";
	if( pv.adj == "fdr" ) pv.adj <- "BH";

	if( missing(gc.robust.rep) || is.na(gc.robust.rep)) gc.robust.rep <- 1;
	if( gc.robust.rep >1 && gc.robust.rep<3 ) gc.robust.rep <- 3;
	
	if( missing( ncores) ) ncores <- 1;
	if( missing( gc.groups) ) gc.groups <- 1;
	if( missing( gc.correction) ) gc.correction <- FALSE;
	if( missing( background.order ) ) background.order <- 2;
	if( missing( background.length ) ) background.length <- 100000;
	if( missing( use.cluster) ) use.cluster <- FALSE;
	if( use.cluster && NROW(tfbs@cluster)==0 ) 
		stop("No cluster information in the tfbs object");
	if( use.cluster && NCOL(tfbs@cluster)!=3 ) 
		stop("No selected motif for each cluster in the tfbs object");
	
	if( use.cluster ) 
	{
		r.mat <- range( tfbs@cluster[,1] );
		if( r.mat[1]<1 || r.mat[2]>tfbs@ntfs )
			stop("The first column of 'tfbs@cluster' exceeds the range of motif data set.");

		cluster.mat <- tfbs@cluster[ tfbs@cluster[,3]==1, c(1,2), drop=F];
	}
	else
		cluster.mat <- cbind( 1:tfbs@ntfs, 1);
	
	if( missing(negative.bed) ) 
	{
		negative.bed <- background.generate( positive.bed );
		cat("*", NROW(negative.bed),  "GC negative loci are randomly generated.\n");  
	}
	
	r.comp <- comparative_scanDb_rtfbs( tfbs, 
							file.twoBit, 
							positive.bed, 
							negative.bed, 
							file.prefix, 
							cluster.mat = cluster.mat,
							ncores = ncores,
							gc.correction = gc.correction,
							gc.correction.pdf = gc.correction.pdf, 
							gc.correction.verbose = TRUE, 
							threshold = threshold , 
							threshold.type = threshold.type, 
							gc.groups = gc.groups,
							background.order = background.order, 
							background.length = background.length, 
							pv.adj=pv.adj ); 

	if( is.null(r.comp) ) return(NULL);
	ret <- r.comp$ret;
	
	if( !is.null(r.comp$bg.sample)  && gc.robust.rep > 1 ) 
	{
		ret.list <- list();
		ret.list[[1]] <- ret;
		
		for(i in 2:gc.robust.rep)
		{
			cat("* GC robust replication for background resampling, loop=", i, "\n");  
			r.comp <- comparative_scanDb_rtfbs( tfbs, 
								file.twoBit, 
								positive.bed, 
								negative.bed, 
								file.prefix, 
								cluster.mat = cluster.mat,
								ncores = ncores,
								gc.correction = gc.correction,
								gc.correction.pdf = gc.correction.pdf, 
								gc.correction.verbose = FALSE, 
								threshold = threshold , 
								threshold.type = threshold.type, 
								gc.groups = gc.groups,
								background.order = background.order, 
								background.length = background.length, 
								pv.adj=pv.adj ); 
			ret.list[[i]] <- r.comp$ret;								
		}						
		
		for(i in 1:NROW(ret))
		{
			Npos.list <- c();
			Nneg.list <- c();
			
			for(k in 1:gc.robust.rep)
			{
				Npos.list <- c( Npos.list, ret.list[[k]][i,'Npos'] );
				Nneg.list <- c( Nneg.list, ret.list[[k]][i,'Nneg'] );
			}
			
			Npos <- median( Npos.list );
			Nneg <- median( Nneg.list );
			n.gc.pos <- ret.list[[1]][i,'gc.pos'];
			n.gc.neg <- ret.list[[1]][i,'gc.neg'];
			Nneg.expected <- Nneg * n.gc.pos/ n.gc.neg ; 
			
			Npos <- ifelse( Npos==0, 1, Npos );
			Nneg <- ifelse( Nneg==0, 1, Nneg );
			
			fe.ratio <- ( Npos / n.gc.pos )/( Nneg / n.gc.neg );

			pos.left <-  n.gc.pos - Npos;
			if( pos.left<0 ) pos.left  <- 0;
			neg.left <- n.gc.neg - Nneg;
			if( neg.left<0 ) neg.left  <- 0;

			tbl  = rbind( c( Npos, Nneg ), c( pos.left, neg.left ) )
			pval = fisher.test(tbl)$p.value;
			
			ret[i, 'Npos']     <- Npos;
			ret[i, 'Nneg']     <- Nneg;
			ret[i, 'expected'] <- round( Nneg.expected, 1);
			ret[i, 'fe.ratio'] <- fe.ratio;
			ret[i, 'pvalue']   <- pval;
		}
		
		ret$pv.adj <- adjust.pvale( ret$pvalue, cluster.mat, pv.adj );
	}
	
	ret$Nneg   <- NULL;
	ret$gc.pos <- NULL;
	ret$gc.neg <- NULL;
	
	if(!is.null(ret))	
	{
		if( NROW(tfbs@tf_info) > 0 )
		{
			tf.motifid   <- unlist(strsplit(as.character( ret$tf.name ), "@"))[seq(1, length(ret$tf.name)*2-1, 2)];
			tf.idx       <- match( tf.motifid, tfbs@tf_info$Motif_ID );
			ret$tf.name  <- tfbs@tf_info$TF_Name[ tf.idx ];
			ret$motif.id <- tfbs@tf_info$Motif_ID[ tf.idx ];
		}
	}
	
	r.parm <- list( file.twoBit   = file.twoBit, 
				file.prefix       = file.prefix, 
				use.cluster       = use.cluster, 
				cluster.mat       = cluster.mat, 
				ncores            = ncores,
				gc.robust.rep     = gc.robust.rep,
				threshold.type    = threshold.type, 
				threshold         = threshold, 
				pv.adj            = pv.adj, 
				background.order  = background.order, 
				background.length = background.length, 
				gc.correction     = gc.correction );
	
	r <- list( result = ret, parm = r.parm);
	
	class(r) <- c("tfbs.enrichment");

	return(r);
}

print.tfbs.enrichment<-function( x, ..., pv.threshold=0.05, pv.adj=NA )
{
	r.comp <- x;

	if (is.na(pv.adj)) pv.adj <- r.comp$parm$pv.adj;

	cat("GC correction:", r.comp$parm$gc.correction, "\n");
	cat("p-value correction:",  pv.adj, "\n");
	cat("Significant p-value:", pv.threshold, "\n");
	cat("Threshold type :", r.comp$parm$threshold.type, "\n");
	cat("Threshold:", r.comp$parm$threshold, "\n");
	cat("Background.order:", r.comp$parm$background.order, "\n");
	cat("Background.length:", r.comp$parm$background.length, "\n");
	cat("GC robust replication:", r.comp$parm$gc.robust.rep, "\n");

	cat("Total Motif:", NROW(r.comp$result), "\n");
	cat("\nSignificant Motifs(or top 20):\n");

	r.comp.sig <- summary.tfbs.enrichment(r.comp, pv.threshold=pv.threshold, pv.adj=pv.adj);
	if(NROW(r.comp.sig)>20)
		r.comp.sig <- r.comp.sig[c(1:20),];

	show(r.comp.sig);		
}

adjust.pvale<-function( r.pvalue, cluster.mat, pv.adj )
{
	# If the cluster is used,...No cluster info ==>  cluster.mat[,2]=1
	cluster.id <- unique( cluster.mat[,2] );
	for(i in 1:length(cluster.id))
	{
		cluster.set <- which( cluster.mat[,2] == cluster.id[i] )

		# if the cluster index is called from the compare function, the p.adjust will be used to the cluster range, not all results.
		r.pvalue[cluster.set] <- p.adjust( r.pvalue[cluster.set], method=pv.adj );
	}
		
	return(r.pvalue);
}

summary.tfbs.enrichment<-function( object, pv.threshold=0.05, pv.adj=NA, ...)
{
	r.comp <- object;

	r <- r.comp$result;
	if(!is.na(pv.adj))
		r$pv.adj <- adjust.pvale( r$pvalue, r.comp$parm$cluster.mat, pv.adj )
	
	r.comp.sum <- r[ order( r$pv.adj), c("motif.id","tf.name","Npos","expected","pv.adj","fe.ratio") ];
	
	colnames(r.comp.sum) <- c("Motif id","TF name","Npos","Expected","pv.adj","Fold enrichment");

	return( r.comp.sum );
}

tfbs.reportEnrichment<-function( tfbs, 
								 r.comp, 
								 file.pdf = NA, 
								 report.size = "letter", 
								 report.title = "", 
								 enrichment.type = c ("both", "enriched", "depleted"),
								 sig.only = TRUE, 
								 pv.threshold = 0.05, 
								 pv.adj = NA,
								 sorted = c ("pvalue", "enrich.ratio"))
{
	stopifnot(class(tfbs) == "tfbs" && class(r.comp) == "tfbs.enrichment" )

	if(!missing(pv.adj))  
		r.comp$result$pv.adj <- adjust.pvale( r.comp$result$pvalue, r.comp$parm$cluster.mat, pv.adj )

	if( !missing(sorted)) 
		sorted <- match.arg( sorted )
	else			
		sorted <- "pvalue";

	r.comp.sel <- r.comp$result;
	if( sorted == "enrich.ratio") 
		r.comp.sel <- r.comp.sel[ order(r.comp.sel$fe.ratio, decreasing = TRUE), c("motif.id","tf.name","Npos","expected","pv.adj","fe.ratio"), drop=F ]
	else
		r.comp.sel <- r.comp.sel[ order(r.comp.sel$pvalue), c("motif.id","tf.name","Npos","expected","pv.adj","fe.ratio"), drop=F ];

	if( !missing(enrichment.type)) 
		enrichment.type <- match.arg(enrichment.type)
	else			
		enrichment.type <- "both";

	idx.sel <- 1:NROW(r.comp.sel);
	if( sig.only )
	{
		idx.sel <- which( r.comp.sel$pv.adj <= pv.threshold );

		if (enrichment.type=="enriched")
			idx.sel <- which( r.comp.sel$pv.adj <= pv.threshold & r.comp.sel$fe.ratio >= 1);

		if (enrichment.type=="depleted")
			idx.sel <- which( r.comp.sel$pv.adj <= pv.threshold & r.comp.sel$fe.ratio < 1);
	}
	else
	{
		if (enrichment.type=="enriched")
			idx.sel <- which( r.comp.sel$fe.ratio >= 1);

		if (enrichment.type=="depleted")
			idx.sel <- which( r.comp.sel$fe.ratio < 1);
	}
	
	if( length(idx.sel) > 0 )
		r.comp.sel <- r.comp.sel[ idx.sel, ,drop=F ]
	else	
	{
		cat( "! No motif information for the report.\n" );
		return(invisible(NULL));
	}	

	r.comp.new <- data.frame( No=c(1:NROW(r.comp.sel)),
							  r.comp.sel[, c("motif.id","tf.name","Npos","expected","pv.adj","fe.ratio","motif.id"), drop=F ] );
	
	df.style <- data.frame(
			position = numeric(0), 
			width    = numeric(0), 
			header   = character(0),     
			hjust    = character(0), 
			style    = character(0), 
			extra1   = character(0), 
			extra2   = character(0), 
			extra3   = character(0),
			extra4   = character(0));
			
	df.style <- rbind(df.style, 
			data.frame( position = 0.00, 
						width    = 0.04, 
						header   = "No.",        
						hjust    = "left",   
						style    = "text", 
						extra1   = "0",  
						extra2   = "0",  
						extra3   = "0",
						extra4   = "0"));
	df.style <- rbind(df.style, 
			data.frame( position = 0.04, 
						width    = 0.10, 
						header   = "Motif ID",   
						hjust    = "left",   
						style    = "text", 
						extra1   = "0",  
						extra2   = "0",  
						extra3   = "0",
						extra4   = "0"));
	df.style <- rbind(df.style, 
			data.frame( position = 0.14, 
						width    = 0.09, 
						header   = "TF Name",    
						hjust    = "left",   
						style    = "text", 
						extra1   = "0",  
						extra2   = "0",  
						extra3   = "0",
						extra4   = "0"));
	df.style <- rbind(df.style, 
			data.frame( position = 0.23, 
						width    = 0.06, 
						header   = "N. Pos.",    
						hjust    = "centre", 
						style    = "text", 
						extra1   = "0",  
						extra2   = "0",  
						extra3   = "0",
						extra4   = "0"));
	df.style <- rbind(df.style, 
			data.frame( position = 0.29, 
						width    = 0.06, 
						header   = "Expected.",  
						hjust    = "centre", 
						style    = "text", 
						extra1   = "0",  
						extra2   = "0",  
						extra3   = "0",
						extra4   = "0"));
	df.style <- rbind(df.style, 
			data.frame( position = 0.35, 
						width    = 0.08, 
						header   = "p-value",    
						hjust    = "centre", 
						style    = "bar",  
						extra1   = "1e-6", 
						extra2   = "1",
						extra3   = "1",
						extra4   = "1"));
	df.style <- rbind(df.style, 
			data.frame( position = 0.43, 
						width    = 0.08, 
						header   = "*Enrichment",
						hjust    = "centre", 
						style    = "bar",  
						extra1   = "-10",
						extra2   = "10", 
						extra3   = "0",
						extra4   = "2"));
	df.style <- rbind(df.style, 
			data.frame( position = 0.51, 
						width    = 0.49,
						header   = "Motif Logo", 
						hjust    = "centre", 
						style    = "logo", 
						extra1   = "0",  
						extra2   = "0",  
						extra3   = "0",
						extra4   = "0"));
	
	output_motif_report( tfbs, 
						r.comp.new, 
						file.pdf, 
						report.size, 
						report.title, 
						df.style, 
						"*Enrichment: Enrichment ratio for positive reads against negative reads." );
}
