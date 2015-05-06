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
comparative_scan_rtfbs.ms <- function(pwm, positive.ms, negative.ms, background.ms, background.mm, fdr = 0.1, threshold = NA, calc.empirical.pvalue = FALSE) {
  if (is.na(threshold)) {
    pos.sites = score.ms(positive.ms, pwm, background.mm, threshold = 0)
    neg.sites = score.ms(negative.ms, pwm, background.mm, threshold = 0)
    bg.sites = score.ms(background.ms, pwm, background.mm, threshold = 0)
    
    fdrtbl = calc.fdr(concat.ms(positive.ms, negative.ms), rbind(pos.sites, neg.sites), background.ms, bg.sites)
    
    ##
    if (is.null(fdrtbl))
      return(NULL)
    
    # find first score threshold below FDR threshold
    idxs = which(fdrtbl$fdr <= fdr)
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
  } else {
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
      
      return(list(Nneg = Nneg, Npos = Npos, thresh = threshold, sites = pos.sites, empirical.pvalues = epvals))
    }
    return(list(Nneg = Nneg, Npos = Npos, thresh = threshold, sites = pos.sites))
  }
}

# NOTE: empirical p-values are computed under the assumption that the negative set
#       is not bound by the specified factor (but is accessible)
comparative_scan_rtfbs <- function(pwm, file.twoBit, positive.bed, negative.bed, fdr = 0.1, threshold = NA, background.order = 2, background.length = 100000, calc.empirical.pvalue = FALSE) {
  # read sequences
  positive.ms = read.seqfile.from.bed(positive.bed, file.twoBit)
  negative.ms = read.seqfile.from.bed(negative.bed, file.twoBit)

  # compute GC quantile masks
  gcBins = joint.gc.quantile.bins(positive.ms, negative.ms)

  # get per GC quantile results
  result.sites = NULL
  Npos = 0
  Nneg = 0
  thresh = NULL
  epvals = NULL
  
  for (idx in 1:4) {
    # compute background model
    pos.ms.i = positive.ms[gcBins$bins1 == idx]
    neg.ms.i = negative.ms[gcBins$bins2 == idx]
    both.ms = concat.ms(pos.ms.i, neg.ms.i)
        
    bg.mm = build.mm(both.ms, background.order)
    
    # generate background sequences
    bg.ms <- simulate.ms(bg.mm, background.length)
  
    # collect sequences
    res = comparative_scan_rtfbs.ms(pwm, pos.ms.i, neg.ms.i, bg.ms, bg.mm, fdr = fdr, threshold = threshold, calc.empirical.pvalue = calc.empirical.pvalue)
    
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
    result = list(Npos = Npos, Nneg = Nneg, assoc.pvalue = pval, thresh = thresh, sites = bed, empirical.pvalues = epvals)
  else
    result = list(Npos = Npos, Nneg = Nneg, assoc.pvalue = pval, thresh = thresh, sites = bed)
}

tfbs_to_bed <- function(sites, tf.name) {
  if (NROW(sites) == 0) ##dim(sites)[1] == 0) ## Alternative formulation fails if sites is empty.
      return(NULL)
  
  spl <- strsplit(as.character(sites$seqname), ":|-")

  chroms <- sapply(spl, function(pair) pair[1])
	starts <- as.integer(sapply(spl, function(pair) pair[2]))
  
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

comparative_scanDb_rtfbs <- function( tfbs, file.twoBit, positive.bed, negative.bed, file.prefix=NA, 
	usemotifs=NA, background.correction=FALSE, fdr = 0.1, threshold = NA, background.order = 2, background.length = 100000, ncores = 3) {
  
  stopifnot(class(tfbs) == "tfbs")
  
  if( !is.na(file.prefix))
  	if( !check_folder_writable( file.prefix ) ) 
  	  stop(paste("Can not create files starting with the prefix:", file.prefix));

  if(NROW(negative.bed) < 50)
	stop(paste("Negative BED file contains only ", NROW(negative.bed), " entries.  Strongly suggest size is increased."))
  
  # read sequences
  positive.ms = read.seqfile.from.bed(positive.bed, file.twoBit)
  negative.ms = read.seqfile.from.bed(negative.bed, file.twoBit)
	
  # detect the difference of gcContent between positive and negative TREs, 
  # if the difference is significant, make a correction for the negative TREs 
  # based on the resampling method.	
  r.bgchk <- background.check( positive.ms, negative.ms, background.correction, file.prefix )
  if(!is.null(r.bgchk))
  {
    negative.bed <- negative.bed[ r.bgchk, ]; 
  	negative.ms  <- read.seqfile.from.bed( negative.bed, file.twoBit);
  }	

  # compute CG quantile masks
  gcBins = joint.gc.quantile.bins(positive.ms, negative.ms)
  
  # partition sequences
  seqs.subsets = lapply(1:4, function(idx) {
    pos.ms.i = positive.ms[gcBins$bins1 == idx]
    neg.ms.i = negative.ms[gcBins$bins2 == idx]
    
    list(pos = pos.ms.i, neg = neg.ms.i)
  })
  
  # compute background models
  gcBacks = lapply(1:4, function(idx) {
    both.ms = concat.ms(seqs.subsets[[idx]]$pos, seqs.subsets[[idx]]$neg)
    
    build.mm(both.ms, background.order)
  })
  
  # compute background sequences
  gcBacks.ms = lapply(1:4, function(idx) {
    simulate.ms(gcBacks[[idx]], background.length)
  })
   
  if(missing(usemotifs)) usemotifs<-c(1:tfbs@ntfs); 
  
  # iterate over TF set
  binding_all <- mclapply(usemotifs, function(i, ...) {
	  # get PWM information
	  pwm = tfbs@pwm[[i]]
    pwm.name = paste(tfbs@mgisymbols[i], "@", i, sep='')

    # scan sequences, per GC bin
    result.sites = NULL
    Npos = 0
    Nneg = 0
  
    for (idx in 1:4) {
      pos.ms.i = seqs.subsets[[idx]]$pos
      neg.ms.i = seqs.subsets[[idx]]$neg
      bg.i.mm = gcBacks[[idx]]
      bg.i.ms = gcBacks.ms[[idx]]
      
      res = comparative_scan_rtfbs.ms(pwm, pos.ms.i, neg.ms.i, bg.i.ms, bg.i.mm, fdr = fdr, threshold = threshold)
      
      if (!is.null(res)) {
        result.sites = rbind(result.sites, res$sites)
        Npos = Npos + res$Npos
        Nneg = Nneg + res$Nneg
      }
    }
    
    # process result
    result.bed = tfbs_to_bed(result.sites, pwm.name)
    
    tbl = rbind(c(Npos, Nneg), c(NROW(positive.bed) - Npos, NROW(negative.bed) - Nneg))
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
  	es.ratio <- (Npos.tmp/NROW(positive.bed) )/(Nneg.tmp/NROW(negative.bed));
    
    pv.bonferroni <- pval*length(usemotifs);
    if( pv.bonferroni > 1)  pv.bonferroni<-1;
    
    # return info
    return(data.frame(tf.name = pwm.name, Npos = Npos, Nneg = Nneg, es.ratio= es.ratio, assoc.pvalue = pval, pv.bonferroni=pv.bonferroni, starch = starch.file))
  }, mc.cores = ncores)
  
  # recombine results
  r.df <- do.call("rbind", binding_all);

  if( missing(file.prefix) || is.na(file.prefix) )  r.df <- subset(r.df, select=-starch);
  
  return(r.df);
}

tfbs_compareTFsite<-function( tfbs, file.twoBit, positive.bed, negative.bed, file.prefix=NA, 
	usemotifs=NA, background.correction = FALSE, fdr = 0.1, threshold = NA, background.order = 2, background.length = 100000, ncores = 3) 
{
    stopifnot(class(tfbs) == "tfbs")

	if( missing( fdr) ) fdr <- 0.1;
	if( missing( threshold ) ) threshold <- NA;
	if( missing( background.order ) ) background.order <- 2;
	if( missing( background.length ) ) background.length <- 100000;
	if( missing( ncores) ) ncores <- 3;
	if( missing( background.correction) ) background.correction <- FALSE;
	if( missing( usemotifs) ) usemotifs <- c(1:tfbs@ntfs);

	r <- comparative_scanDb_rtfbs( tfbs, 
		file.twoBit, 
		positive.bed, 
		negative.bed, 
		file.prefix, 
		usemotifs = usemotifs,
		background.correction = background.correction,
		fdr = fdr , 
		threshold = threshold , 
		background.order = background.order, 
		background.length = background.length, 
		ncores = ncores ); 
		
	if(!is.null(r))	
	{
		if(!is.null(tfbs@extra_info))
		{
			tf.motifid <- unlist(strsplit(as.character(r$tf.name), "@"))[seq(1, length(r$tf.name)*2-1, 2)];
			tf.idx <- match( tf.motifid, tfbs@extra_info$Motif_ID );
			r$tf.name <- tfbs@extra_info$TF_Name[tf.idx];
			r$motif.id <-tfbs@extra_info$Motif_ID[tf.idx];
		}
	}
	
	r;
}

background.check<-function( positive.ms, negative.ms, background.correction, file.prefix=NA )
{
	gc.pos <- gcContent.ms(positive.ms);
	gc.neg <- gcContent.ms(negative.ms);
	
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
	if( background.correction==FALSE || gc.test$p.value > 0.01)
		return(NULL);
		
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

	## Make n.sample as large as possible, up to 10k sequences.  More BG gains statistical power.
	n.sample <- min(c(10000, round(length(gc.neg)/10)))

	indx.bgnew <- resample( gc.pos, gc.neg, n=n.sample)

	gc.test2 <- wilcox.test(gc.pos, gc.neg[indx.bgnew], conf.int=TRUE, conf.level=0.9 );

	cat("* After the resampling from negative TREs, p-value of Wilcox test:", gc.test2$p.value, "\n" );

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

