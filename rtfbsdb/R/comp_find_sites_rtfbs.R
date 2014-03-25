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
  
  cutoffs = quantile(c(gc1, gc2), seq(from = 0, to = 1, length.out = n + 1), names = FALSE)
  bin1 <- findInterval(gc1, cutoffs, all.inside = TRUE)
  bin2 <- findInterval(gc2, cutoffs, all.inside = TRUE)
  
  # check each bin has at least one sequence
  n.seqs.per.bin = sapply(1:n, function(i) min(sum(bin1 == i), sum(bin2 == i)))
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
comparative_scan_rtfbs <- function(pwm, positive.bed, negative.bed, fdr = 0.1, threshold = NA, background.order = 2, background.length = 100000, twoBit_path= "/gbdb/hg19/hg19.2bit", calc.empirical.pvalue = FALSE) {
  # read sequences
  positive.ms = read.seqfile.from.bed(positive.bed, twoBit_path)
  negative.ms = read.seqfile.from.bed(negative.bed, twoBit_path)

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
  ord = order(bed[,1], bed[,2])

  if (calc.empirical.pvalue)
    result = list(Npos = Npos, Nneg = Nneg, assoc.pvalue = pval, thresh = thresh, sites = bed[ord,], empirical.pvalues = epvals[ord])
  else
    result = list(Npos = Npos, Nneg = Nneg, assoc.pvalue = pval, thresh = thresh, sites = bed[ord,])
}

tfbs_to_bed <- function(sites, tf.name) {
  if (dim(sites)[1] == 0)
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
  # sort BED
  ord = order(bed[,1], bed[,2])
  
  # pipe bed into starch file
  write.table(bed[ord, ], file = pipe(paste("starch - >", filename)), 
    quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

comparative_scanDb_rtfbs <- function(tfbs, positive.bed, negative.bed, file_prefix, fdr = 0.1, threshold = NA, background.order = 2, background.length = 100000, twoBit_path= "/gbdb/hg19/hg19.2bit", ncores = 3) {
  stopifnot(class(tfbs) == "tfbs")
  
  # read sequences
  positive.ms = read.seqfile.from.bed(positive.bed, twoBit_path)
  negative.ms = read.seqfile.from.bed(negative.bed, twoBit_path)
  
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
  
  # iterate over TF set
  binding_all <- mclapply(tfbs@usemotifs, function(i, ...) {
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
    
    tbl = rbind(c(Npos, Nneg), c(dim(positive.bed)[1] - Npos, dim(negative.bed)[1] - Nneg))
    pval = fisher.test(tbl)$p.value
    
    # save sites
    starch.file = paste(file_prefix, ".", i, ".bed.starch", sep='')
    write.starchbed(result.bed, starch.file)
    
    # return info
    return(data.frame(tf.name = pwm.name, Npos = Npos, Nneg = Nneg, assoc.pvalue = pval, starch = starch.file))
  }, mc.cores = ncores)
  
  # recombine results
  do.call("rbind", binding_all)
}
