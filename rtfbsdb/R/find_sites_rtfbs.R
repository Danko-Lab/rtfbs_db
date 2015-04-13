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
#` @param dnase_peaks_bed bed-formatted peak information.
#` @param motif_path path to the motif PWM file.
#` @param divide_num List of parameters for all data types, for the model representing no TF binding.
#` @return List structure representing the match score to the motif.
scan_rtfbs <- function(tf_name, file.twoBit, dnase_peaks_bed, motif_path, return.posteriors=TRUE, ...) {
  ## Read the pwm and sequence file.
  motif <- read.motif(motif_path, header=TRUE) 
  
  ## Write out new fasta file, adding on half-width of the motif, to correctly align motif center ...
  half_width <- ceiling((NROW(motif)-1)/2)
  extBed = extend.bed(dnase_peaks_bed, half_width - 1)
  dnase_peaks = read.seqfile.from.bed(extBed, file.twoBit)

  ## Swiched rtfbs to returning posteriors.
  bgModel <- build.mm(dnase_peaks, 3)

  ## Switch to rtfbsPost to return proper data structure.
  binding <- score.ms(dnase_peaks, motif, bgModel, return_posteriors=return.posteriors, ...) 
  
  if(return.posteriors == FALSE) { ## Parse binding site into genomic coordinates.
    spl <- strsplit(as.character(binding$seqname), ":|-")
    peak_chrom <- as.character(sapply(c(1:NROW(binding)), function(x) {spl[[x]][[1]]}))
    peak_start <- as.integer(sapply(c(1:NROW(binding)), function(x) {spl[[x]][[2]]}))
	binding <- data.frame(chrom= peak_chrom, 
							chromStart= peak_start+ binding$start, 
							chromEnd= peak_start+ binding$end, 
							name= binding$seqname,
							score= binding$score,
							strand= binding$strand)
  }
  
  return(binding)
}

## return.type --> c("matches", "posteriors", "maxposterior", or "writedb")
##
## matches 		-- returns all matching motifs.
## writedb 		-- writes a bed file with matches.  Assuems that sort-bed and starch tools are availiable in $PATH
## posteriors 	-- returns the posteriors at each position.
## maxposterior	-- returns the max(posterior) in each dnase-1 peak.

scanDb_rtfbs <- function(tfbs, file.twoBit, dnase.peaks.bed, file.prefix= "data.db", usemotifs=NA, ncores= 3, return.type="matches", threshold= 6, ...) {
  stopifnot(class(tfbs) == "tfbs")

  ## Read in the DNAse-1 peaks ...
  half_width=15 ## Max size of TF in set of 1800 is 30 (half-width = 15).
  options("scipen"=100, "digits"=4)
  extBed = extend.bed(dnase.peaks.bed, half_width - 1)
  dnase_peaks = read.seqfile.from.bed(extBed, file.twoBit)
  
  ## Swiched rtfbs to returning posteriors.
  bgModel <- build.mm(dnase_peaks, 3)
  return.posteriors <- (return.type=="posteriors"|return.type=="maxposterior")

  binding_all <- mclapply(usemotifs, function(i, ...) {
	  ## Read the pwm and sequence file.
	  motif <- tfbs@pwm[[i]]

	  ## Switch to rtfbsPost to return proper data structure.
	  binding <- score.ms(dnase_peaks, motif, bgModel, return_posteriors=return.posteriors, threshold=threshold, ...)

      if(return.type == "maxposterior" & NROW(binding) > 0) { ## Parse binding site into genomic coordinates.
        return(sapply(1:NROW(binding), function(x) { max(c(binding[[x]]$MotifModel$Forward - binding[[x]]$Background, binding[[x]]$MotifModel$Reverse - binding[[x]]$Background)) }))
      }

	  if(return.posteriors == FALSE & NROW(binding) > 0) { ## Parse binding site into genomic coordinates.
		spl <- strsplit(as.character(binding$seqname), ":|-")
		peak_chrom <- as.character(sapply(c(1:NROW(binding)), function(x) {spl[[x]][[1]]}))
		peak_start <- as.integer(sapply(c(1:NROW(binding)), function(x) {spl[[x]][[2]]}))
		binding <- data.frame(chrom= peak_chrom, 
								chromStart= peak_start+ binding$start- 1,  ## -1 determined empirically.
								chromEnd= peak_start+ binding$end, 
								name= paste(tfbs@mgisymbols[i],"@",i,sep=""), # binding$seqname,
								score= binding$score,
								strand= binding$strand)
		if(return.type == "writedb") {
		  write.table(binding, file= pipe(paste(" sort-bed - | starch - > ",file.prefix,i,".bed.tmp.starch", sep="")),  quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
		  return(paste(file.prefix,i,".bed.tmp.starch ", sep=""))
		}
	  }
	  return(binding)
  }, mc.cores= ncores)

  if(return.type == "writedb") {
   cat_files <- paste(unlist(binding_all), collapse=" ")
   system(paste("starchcat ", cat_files, " > ",file.prefix,".db.starch", sep=""))
   system(paste("rm ",file.prefix,"*.bed.tmp.starch",sep=""))
   binding_all <- paste(file.prefix,".db.starch", sep="")
  }
  if(return.type == "maxposterior") {
   binding_all <- matrix(unlist(binding_all), nrow= NROW(dnase_peaks), ncol= NROW(usemotifs))
  }

  return(binding_all)
}

#ncores=3 for 4 cores CPU.
# Under 1 node with 1 task mcapply can not be going will.
 

tfbs_scanTFsite<-function(tfbs, file.twoBit, dnase.peaks.bed=NULL, file.prefix="scan.db",  usemotifs=NA, ncores= 3, return.type="matches", threshold= 6, ... ) 
{

  stopifnot(class(tfbs) == "tfbs")

	if( missing(dnase.peaks.bed) )
	{
		chromInfo <- get_chromosome_size(file.twoBit);

		offset_dist <- 250;
		chromInfo <- chromInfo[grep("_|chrM|chrY|chrX", chromInfo[,1], invert=TRUE),];
		dnase.peaks.bed <- data.frame(chrom=chromInfo[,1], chromStart=rep(0)+offset_dist, chromEnd=(chromInfo[,2]-1-offset_dist));
	}
	
	if( missing(file.prefix) ) file.prefix="data.db";
	if( missing(return.type) ) return.type="matches";
	if( missing(ncores) ) ncores= 3;
	if( missing(threshold) ) threshold= 6;
	if( missing(usemotifs)) usemotifs =c(1:tfbs@ntfs);

	r <- scanDb_rtfbs( tfbs, file.twoBit, dnase.peaks.bed, file.prefix=file.prefix, usemotifs=usemotifs, ncores=ncores, return.type=return.type, threshold= threshold, ...); 
	r;
}
