##
## Comparative TFBS search
##
##

#' Create GC bins that span the joint set GC content distribution.
#'
#' @param set1.ms first set of sequences
#' @param set2.ms second set of sequences
#' @param n integer number of bins (default 4)
#' @return list with two integer vectors, bins1 and bins2, indicating bin assignments for each set.
joint.gc.quantile.bins <- function(set1.ms, set2.ms, n = 4) {
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
comparative_scan_rtfbs.ms <- function(pwm, positive.ms, negative.ms, background.ms, background.mm, fdr.threshold = 0.1, score.threshold = NA, calc.empirical.pvalue = FALSE) 
{
  if (!is.na(fdr.threshold)) 
  {
    pos.sites = score.ms(positive.ms, pwm, background.mm, threshold = 0)
    neg.sites = score.ms(negative.ms, pwm, background.mm, threshold = 0)
    bg.sites = score.ms(background.ms, pwm, background.mm, threshold = 0)
    
    fdrtbl = calc.fdr(concat.ms(positive.ms, negative.ms), rbind(pos.sites, neg.sites), background.ms, bg.sites)
    
    ##
    if (is.null(fdrtbl))
      return(NULL)
    
    # find first score threshold below FDR threshold
    idxs = which(fdrtbl$fdr <= fdr.threshold)
    if (length(idxs) == 0)
      return(NULL)
    
    thresh = fdrtbl[idxs[length(idxs)], 1]
    
    pos.mask = pos.sites$score >= thresh
    neg.mask = neg.sites$score >= thresh
    
    Npos = length(unique(pos.sites$seqname[pos.mask]))
    Nneg = length(unique(neg.sites$seqname[neg.mask]))
    
    # NOTE: empirical p-values are computed under the assumption that the negative set
    #       is not bound by the specified factor (but is accessible)
    if (calc.empirical.pvalue) {
      epvals = sapply(pos.sites$score[pos.mask], function(thresh) {
        K = sum(neg.sites$score >= thresh)
        Z = length(neg.sites$score)

        (K + 1)/(Z + 1) # empirical p-value
      })
      
      return(list(Nneg = Nneg, Npos = Npos, thresh = thresh, sites = pos.sites[pos.mask, ], empirical.pvalues = epvals))
    }
    
    return(list(Nneg = Nneg, Npos = Npos, thresh = thresh, sites = pos.sites[pos.mask, ]))
  } 
  else 
  {
    pos.sites = score.ms(positive.ms, pwm, background.mm, threshold = score.threshold)
    neg.sites = score.ms(negative.ms, pwm, background.mm, threshold = score.threshold)
    
    Npos = length(unique(pos.sites$seqname))
    Nneg = length(unique(neg.sites$seqname))
  
    if (calc.empirical.pvalue) {
      epvals = sapply(pos.sites$score, function(thresh) {
        K = sum(neg.sites$score >= thresh)
        Z = length(neg.sites$score)

        (K + 1)/(Z + 1) # empirical p-value
      })
      
      return(list(Nneg = Nneg, Npos = Npos, thresh = score.threshold, sites = pos.sites, empirical.pvalues = epvals))
    }
    
    return(list(Nneg = Nneg, Npos = Npos, thresh = score.threshold, sites = pos.sites))
  }
}

# NOTE: empirical p-values are computed under the assumption that the negative set
#       is not bound by the specified factor (but is accessible)
comparative_scan_rtfbs <- function(pwm, file.twoBit, positive.bed, negative.bed, fdr.threshold = 0.1, score.threshold = NA, gc.groups=4, background.order = 2, background.length = 100000, calc.empirical.pvalue = FALSE) {
  # read sequences
  positive.ms = read.seqfile.from.bed(positive.bed, file.twoBit)
  negative.ms = read.seqfile.from.bed(negative.bed, file.twoBit)

  # compute GC quantile masks
  gcBins = joint.gc.quantile.bins(positive.ms, negative.ms, n=gc.groups)

  # get per GC quantile results
  result.sites = NULL
  Npos = 0
  Nneg = 0
  thresh = NULL
  epvals = NULL
  
  for (idx in 1:gc.groups) 
  {
    # compute background model
    pos.ms.i = positive.ms[gcBins$bins1 == idx]
    neg.ms.i = negative.ms[gcBins$bins2 == idx]
    both.ms = concat.ms(pos.ms.i, neg.ms.i)
        
    bg.mm = build.mm(both.ms, background.order)
    
    # generate background sequences
    bg.ms <- simulate.ms(bg.mm, background.length)
  
    # collect sequences
    res = comparative_scan_rtfbs.ms(pwm, pos.ms.i, neg.ms.i, bg.ms, bg.mm, fdr.threshold = fdr.threshold, score.threshold = score.threshold, calc.empirical.pvalue = calc.empirical.pvalue)
    
    if (!is.null(res)) {
      result.sites = rbind(result.sites, res$sites)
      Npos = Npos + res$Npos
      Nneg = Nneg + res$Nneg
      thresh = c(thresh, res$thresh)
      
      if (calc.empirical.pvalue)
        epvals = c(epvals, res$empirical.pvalues)
    }
  }
  
  # fisher test
  tbl = rbind(c(Npos, Nneg), c(dim(positive.bed)[1] - Npos, dim(negative.bed)[1] - Nneg))
  pval = fisher.test(tbl)$p.value
  
  # merge results
  bed = tfbs_to_bed(result.sites, "pwm")
  if (!is.null(bed) && dim(bed)[1] > 0) {
    ord = order(bed[,1], bed[,2])
    bed = bed[ord,]

    if (calc.empirical.pvalue)
        epvals = epvals[ord]
  }

  if (calc.empirical.pvalue)
    result = list(Npos = Npos, Nneg = Nneg, pvalue = pval, thresh = thresh, sites = bed, empirical.pvalues = epvals)
  else
    result = list(Npos = Npos, Nneg = Nneg, pvalue = pval, thresh = thresh, sites = bed)
}

tfbs_to_bed <- function(sites, tf.name) {
  if (NROW(sites) == 0) ##dim(sites)[1] == 0) ## Alternative formulation fails if sites is empty.
      return(NULL)
  
  spl <- strsplit(as.character(sites$seqname), ":|-")

  chroms <- sapply(spl, function(pair) pair[1]);
  starts <- as.integer(sapply(spl, function(pair) pair[2]));
  
  bed = data.frame(
    chrom = chroms,
    chromStart = as.integer(starts + sites$start - 1), # add site offset
    chromEnd = as.integer(starts + sites$end),
    name = tf.name,
    score = sites$score,
    strand = sites$strand
  )
  
  return(bed)
}

write.starchbed <- function(bed, filename) {
  # pipe bed into starch file
  write.table(bed, file = pipe(paste("sort-bed - | starch - >", filename)), 
    quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

comparative_scanDb_rtfbs <- function( tfbs, file.twoBit, positive.bed, negative.bed, file.prefix = NA, use.cluster = NA, ncores = 3, 
	gc.correction = FALSE, fdr.threshold = 0.1, score.threshold = NA, gc.groups=4, background.order = 2, background.length = 100000, pv.adj = NA ) {
  
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
  
  r.bgchk <- background.check( gc.pos, gc.neg, gc.correction, file.prefix )
  if(!is.null(r.bgchk))
  {
     gc.neg <- gc.neg[ r.bgchk ];
     negative.ms  <- negative.ms [ r.bgchk ];
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
   
  if(missing(use.cluster)) use.cluster <- cbind(1:tfbs@ntfs, 1); 
  
  # iterate over TF set
  binding_all <- mclapply(use.cluster[,1], function(i, ...) {
	# get PWM information
	pwm = tfbs@pwm[[i]]
    
    pwm.name = paste(tfbs@mgisymbols[i], "@", i, sep='')

    # scan sequences, per GC bin
    result.sites = NULL
    Npos = 0
    Nneg = 0
  
    for (idx in 1:gc.groups) {
      pos.ms.i = seqs.subsets[[idx]]$pos
      neg.ms.i = seqs.subsets[[idx]]$neg
      bg.i.mm = gcBacks[[idx]]
      bg.i.ms = gcBacks.ms[[idx]]
      
      res = comparative_scan_rtfbs.ms(pwm, pos.ms.i, neg.ms.i, bg.i.ms, bg.i.mm, fdr.threshold = fdr.threshold, score.threshold = score.threshold)
      
      if (!is.null(res)) {
        result.sites = rbind(result.sites, res$sites)
        Npos = Npos + res$Npos
        Nneg = Nneg + res$Nneg
      }
    }
    
    # process result
    result.bed = tfbs_to_bed(result.sites, pwm.name)
    
    tbl = rbind(c(Npos, Nneg), c(length(gc.pos) - Npos, length(gc.neg) - Nneg))
    pval = fisher.test(tbl)$p.value
    
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
  	es.ratio <- (Npos.tmp/length(gc.pos) )/(Nneg.tmp/length(gc.neg));
    
    # return info
    return(data.frame('motif.id'="", 'tf.name' = pwm.name, 'Npos' = Npos, 'Nneg' = Nneg, 'es.ratio'= es.ratio, 'pvalue' = pval, 'pv.adj'= pval, 'starch' = starch.file))
  }, mc.cores = ncores)
  
  # recombine results
  r.df <- do.call("rbind", binding_all);

  #put this dummy variable here to pass R CMD check rtfbsdb --as-cran	
  starch <- NULL;
  
  if (NROW(r.df)>0)	
  {
	 r.df$pv.adj <- adjust.pvale( r.df$pvalue, use.cluster, pv.adj );
  	 
  	 # the following 'starch' is a column in r.df, not equal NULL;
     if( missing(file.prefix) || is.na(file.prefix) )  r.df <- subset( r.df, select = -starch);
  }
  
  return(r.df);
}

tfbs_compareTFsite<-function( tfbs, file.twoBit, positive.bed, negative.bed, file.prefix = NA, use.cluster = NA, ncores = 3,
	gc.correction = FALSE, fdr = 0.1, threshold = NA, gc.groups=4, background.order = 2, background.length = 100000, pv.adj=p.adjust.methods) 
{
    stopifnot(class(tfbs) == "tfbs")
  
 	if( !is.valid.bed( positive.bed ) ) stop("Wrong format in the parameter of 'positive.bed', at least three columns including chromosome, strat, stop.");
    if( !is.valid.bed( negative.bed ) ) stop("Wrong format in the parameter of 'negative.bed', at least three columns including chromosome, strat, stop.");
 
    if( !missing(pv.adj)) pv.adj <- match.arg(pv.adj)
    if( missing(pv.adj) ) pv.adj <- "bonferroni";
    if( pv.adj == "fdr" ) pv.adj <- "BH";

	if( missing( ncores) ) ncores <- 3;
	if( missing( gc.correction) ) gc.correction <- FALSE;
	if( missing( fdr) ) fdr <- 0.1;
	if( missing( threshold ) ) threshold <- NA;
	if( missing( background.order ) ) background.order <- 2;
	if( missing( background.length ) ) background.length <- 100000;

	if( !missing( use.cluster) ) 
	{
		mat.cluster <- as.matrix( use.cluster[, c(1,2),drop=F ] );
		r.mat <- range( mat.cluster[,1] )
		if( r.mat[1]<1 || r.mat[2]> tfbs@ntfs )
			stop("The first column of 'use.cluster' exceed the range of motif data set.");
	}
	else
		mat.cluster <- cbind( 1:tfbs@ntfs, NA);
	
	ret <- comparative_scanDb_rtfbs( tfbs, 
		file.twoBit, 
		positive.bed, 
		negative.bed, 
		file.prefix, 
		use.cluster = mat.cluster,
		ncores = ncores,
		gc.correction = gc.correction,
		fdr.threshold = fdr , 
		score.threshold = threshold , 
		gc.groups = gc.groups,
		background.order = background.order, 
		background.length = background.length, 
		pv.adj=pv.adj ); 
		
	if(!is.null(ret))	
	{
		if(!is.null(tfbs@extra_info))
		{
			tf.motifid <- unlist(strsplit(as.character( ret$tf.name ), "@"))[seq(1, length(ret$tf.name)*2-1, 2)];
			tf.idx <- match( tf.motifid, tfbs@extra_info$Motif_ID );
			ret$tf.name  <- tfbs@extra_info$TF_Name[ tf.idx ];
			ret$motif.id <- tfbs@extra_info$Motif_ID[ tf.idx ];
		}
	}
	
	r.parm <- list( file.twoBit = file.twoBit, 
				file.prefix = file.prefix, 
				use.cluster  = mat.cluster, 
				ncores      = ncores,
				fdr         = fdr, 
				threshold   = threshold, 
				pv.adj      = pv.adj, 
				background.order    = background.order, 
				background.length   = background.length, 
				gc.correction = gc.correction );
	
	r <- list( result = ret, parm = r.parm);
	
	class(r) <- c("tfbs.comparson");

	return(r);
}

background.check<-function( gc.pos, gc.neg, gc.correction, file.prefix=NA )
{
	gc.test <- wilcox.test(gc.pos, gc.neg, conf.int=TRUE, conf.level=0.9 );
	
	# Actually 0.01 is not good for wilcox.test, 
	if(gc.test$p.value<0.01 )
	{
		cat("! The difference between negative and positive TREs is significant, p-value of Wilcox test:", gc.test$p.value, "\n" );
		
		if( require(vioplot) )
		{
			pdf.file <- paste("vioplot.before.correct", file.prefix, "pdf", sep=".");
			cat("* Please check the vioplot figure to make sure, the vioplot figure: ", pdf.file, "\n" );

			r.try <- try ( pdf(pdf.file) );
			if( class(r.try) != "try-error" )
			{
				vioplot(gc.pos, gc.neg, names=c("Positive", "Negative"));
				abline(h=median(gc.pos), lty="dotted")
				dev.off();
			}
			else
				cat("  Failed to output the vioplot figure for the results.\n");
		}
	}
	
	#No need to do sampling the background data.
	if( gc.correction==FALSE || gc.test$p.value > 0.01)
		return(NULL);

	if( length(gc.neg)< 5000 )
	{
		cat("! Failed to make background correction due to small bed data(size<5000).\n");
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
	ns.sample <- c ( round(length(gc.neg)/c(5,10,15,20,25)), 1000);
	ns.sample <- ns.sample[ ns.sample>=1000 ];
	
	ns.pvalue <- unlist(lapply(ns.sample, try.sample));
	n.sample  <- ns.sample[ which.max(ns.pvalue) ];

	indx.bgnew <- resample( gc.pos, gc.neg, n=n.sample)
	gc.test2 <- wilcox.test(gc.pos, gc.neg[indx.bgnew], conf.int=TRUE, conf.level=0.9 );
	cat("* After the resampling from negative TREs, sampe size:", n.sample, "p-value of Wilcox test:", gc.test2$p.value, "\n" );

	if( require(vioplot) )
	{
		pdf.file <- paste("vioplot.after.correct", file.prefix, "pdf", sep=".");
		cat("* The vioplot figure after correction:", pdf.file, "\n" );

		r.try <- try ( pdf(pdf.file) );
		if( class(r.try) != "try-error")
		{
			vioplot(gc.pos, gc.neg, gc.neg[indx.bgnew], names=c("Positive", "Negative", "Negative.resample")); 
			abline(h=median(gc.pos), lty="dotted")
			dev.off()
		}
		else
			cat("  Failed to output the vioplot figure for the results of correction.\n");
		
	}
	
	## return sampling background.
	return( indx.bgnew );
}

print.tfbs.comparson<-function( x, ..., pv.cutoff=0.05, pv.adj=NA )
{
	r.comp <- x;

	if (is.na(pv.adj)) pv.adj <- r.comp$parm$pv.adj;

	cat("GC correction:", r.comp$parm$gc.correction, "\n");
	cat("p-value correction:",  pv.adj, "\n");
	cat("Significant p-value:", pv.cutoff, "\n");
	cat("TF binding FDR threshold:", r.comp$parm$fdr, "\n");
	cat("TF binding score threshold:", r.comp$parm$threshold, "\n");
	cat("TF binding background.order:", r.comp$parm$background.order, "\n");
	cat("TF binding background.length:", r.comp$parm$background.length, "\n");

	cat("Total Motif:", NROW(r.comp$result), "\n");
	cat("\nSignificant Motifs(or top 20):\n");

	r.comp.sig <- summary.tfbs.comparson(r.comp, pv.cutoff=pv.cutoff, pv.adj=pv.adj);
	if(NROW(r.comp.sig)>20)
		r.comp.sig <- r.comp.sig[c(1:20),];

	show(r.comp.sig);		
}

adjust.pvale<-function( r.pvalue, use.cluster, pv.adj )
{
	# If the cluster is used,...No cluster info ==>  use.cluster[,2]=NA
	cluster.id <- unique( use.cluster[,2] );
	for(i in 1:length(cluster.id))
	{
		cluster.set <- which( use.cluster[,2] == cluster.id[i] )

		# if the cluster index is called from the compare function, the p.adjust will be used to the cluster range, not all results.
		r.pvalue[cluster.set] <- p.adjust( r.pvalue[cluster.set], method=pv.adj );
	}
		
	return(r.pvalue);
}

summary.tfbs.comparson<-function( object, pv.cutoff=0.05, pv.adj=NA, ...)
{
    r.comp <- object;
    
	r <- r.comp$result;
	if(!is.na(pv.adj))
		r$pv.adj <- adjust.pvale( r$pvalue, r.comp$parm$use.cluster, pv.adj )
	
	r.comp.sum <- r[ order( r$pv.adj), c("motif.id","tf.name","Npos","Nneg","pv.adj","es.ratio") ];

	return( r.comp.sum );
}

tfbs.reportComparson<-function( tfbs, r.comp, file.pdf=NA, report.size="letter", report.title="", sig.only=TRUE, pv.cutoff=0.05, pv.adj=NA )
{
	if(!is.na(pv.adj))  
		r.comp$result$pv.adj <- adjust.pvale( r.comp$result$pvalue, r.comp$parm$use.cluster, pv.adj )
	
	r.comp.sel <- r.comp$result;
	r.comp.sel <- r.comp.sel[ order(r.comp.sel$pvalue), c("motif.id","tf.name","Npos","Nneg","pv.adj","es.ratio"), drop=F ];
	
	if( sig.only )
	{
		idx.sel <- which( r.comp.sel$pv.adj <= pv.cutoff );
		if( length(idx.sel) > 0 )
			r.comp.sel <- r.comp.sel[ idx.sel, ]
		else
			r.comp.sel <- r.comp.sel[ -c(1:NROW(r.comp.sel)), , drop=F];
	}
	
	if(NROW(r.comp.sel)==0)
	{
		cat( "! No motif information for the report.\n" );
		return(invisible(NULL));
	}	

	r.comp.new <- data.frame( No=c(1:NROW(r.comp.sel)), r.comp.sel[, c("motif.id","tf.name","Npos","Nneg","pv.adj","es.ratio","motif.id"), drop=F ] );
	
	df.style <- data.frame(position=numeric(0), width=numeric(0), header=character(0), hjust=character(0), style=character(0), extra1=character(0), extra2=character(0), extra2=character(0));
	df.style <- rbind(df.style, data.frame(position=0.00, width=0.04, header="No.",        hjust="left",   style="text", extra1="0",  extra2="0",  extra3="0"));
	df.style <- rbind(df.style, data.frame(position=0.04, width=0.10, header="Motif ID",   hjust="left",   style="text", extra1="0",  extra2="0",  extra3="0"));
	df.style <- rbind(df.style, data.frame(position=0.14, width=0.09, header="TF Name",    hjust="left",   style="text", extra1="0",  extra2="0",  extra3="0"));
	df.style <- rbind(df.style, data.frame(position=0.23, width=0.06, header="N. Pos.",    hjust="centre", style="text", extra1="0",  extra2="0",  extra3="0"));
	df.style <- rbind(df.style, data.frame(position=0.29, width=0.06, header="N. Neg.",    hjust="centre", style="text", extra1="0",  extra2="0",  extra3="0"));
	df.style <- rbind(df.style, data.frame(position=0.35, width=0.08, header="p-value",    hjust="centre", style="bar",  extra1="1e-6", extra2="1",extra3="1"));
	df.style <- rbind(df.style, data.frame(position=0.43, width=0.08, header="*Ratio",     hjust="centre", style="bar",  extra1="-10",extra2="10", extra3="0"));
	df.style <- rbind(df.style, data.frame(position=0.51, width=0.49, header="Motif Logo", hjust="centre", style="logo", extra1="0",  extra2="0",  extra3="0"));
	
	output_motif_report( tfbs, r.comp.new, file.pdf, report.size, report.title, df.style, "*Ratio: Enrichment ratio for positive reads against negative reads." );
}
