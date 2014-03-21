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
scan_rtfbs <- function(tf_name, dnase_peaks_bed, motif_path, twoBit_path= "/gbdb/hg19/hg19.2bit", return_posteriors=TRUE, ...) {
  ## Read the pwm and sequence file.
  motif <- read.motif(motif_path, header=TRUE) 
  
  ## Write out new fasta file, adding on half-width of the motif, to correctly align motif center ...
  half_width <- ceiling((NROW(motif)-1)/2)
  extBed = extend.bed(dnase_peaks_bed, half_width - 1)
  dnase_peaks = read.seqfile.from.bed(extBed, twoBit_path)

  ## Swiched rtfbs to returning posteriors.
  bgModel <- build.mm(dnase_peaks, 3)

  ## Switch to rtfbsPost to return proper data structure.
  binding <- score.ms(dnase_peaks, motif, bgModel, return_posteriors=return_posteriors, ...) 
  
  if(return_posteriors == FALSE) { ## Parse binding site into genomic coordinates.
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

## return_type --> c("matches", "posteriors", "maxposterior", or "writedb")
##
## matches 		-- returns all matching motifs.
## writedb 		-- writes a bed file with matches.  Assuems that sort-bed and starch tools are availiable in $PATH
## posteriors 	-- returns the posteriors at each position.
## maxposterior	-- returns the max(posterior) in each dnase-1 peak.
scanDb_rtfbs <- function(tfbs, dnase_peaks_bed, file_prefix= "data.db", twoBit_path= "/gbdb/hg19/hg19.2bit", ncores= 3, return_type="matches", threshold= 6, ...) {
  stopifnot(class(tfbs) == "tfbs")

  ## Read in the DNAse-1 peaks ...
  half_width=15 ## Max size of TF in set of 1800 is 30 (half-width = 15).
  options("scipen"=100, "digits"=4)
  extBed = extend.bed(dnase_peaks_bed, half_width - 1)
  dnase_peaks = read.seqfile.from.bed(extBed, twoBit_path)
  
  ## Swiched rtfbs to returning posteriors.
  bgModel <- build.mm(dnase_peaks, 3)
  return_posteriors <- (return_type=="posteriors"|return_type=="maxposterior")

  binding_all <- mclapply(tfbs@usemotifs, function(i, ...) {
	  ## Read the pwm and sequence file.
	  motif <- tfbs@pwm[[i]]

	  ## Switch to rtfbsPost to return proper data structure.
	  binding <- score.ms(dnase_peaks, motif, bgModel, return_posteriors=return_posteriors, threshold=threshold, ...)

      if(return_type == "maxposterior" & NROW(binding) > 0) { ## Parse binding site into genomic coordinates.
        return(sapply(1:NROW(binding), function(x) { max(c(binding[[x]]$MotifModel$Forward - binding[[x]]$Background, binding[[x]]$MotifModel$Reverse - binding[[x]]$Background)) }))
      }
	  
	  if(return_posteriors == FALSE & NROW(binding) > 0) { ## Parse binding site into genomic coordinates.
		spl <- strsplit(as.character(binding$seqname), ":|-")
		peak_chrom <- as.character(sapply(c(1:NROW(binding)), function(x) {spl[[x]][[1]]}))
		peak_start <- as.integer(sapply(c(1:NROW(binding)), function(x) {spl[[x]][[2]]}))
		binding <- data.frame(chrom= peak_chrom, 
								chromStart= peak_start+ binding$start- 1,  ## -1 determined empirically.
								chromEnd= peak_start+ binding$end, 
								name= paste(tfbs@mgisymbols[i],"@",i,sep=""), # binding$seqname,
								score= binding$score,
								strand= binding$strand)
		if(return_type == "writedb") {
		  write.table(binding, file= pipe(paste(" sort-bed - | starch - > ",file_prefix,i,".bed.tmp.starch", sep="")),  quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
		  return(paste(file_prefix,i,".bed.tmp.starch ", sep=""))
		}
	  }
	  return(binding)
  }, mc.cores= ncores)

  if(return_type == "writedb") {
   cat_files <- paste(unlist(binding_all), collapse=" ")
   system(paste("starchcat ", cat_files, " > ",file_prefix,".db.starch", sep=""))
   system(paste("rm ",file_prefix,"*.bed.tmp.starch",sep=""))
   binding_all <- paste(file_prefix,".db.starch", sep="")
  }
  if(return_type == "maxposterior") {
   binding_all <- matrix(unlist(binding_all), nrow= NROW(dnase_peaks), ncol= NROW(tfbs@usemotifs))
  }

  return(binding_all)
}

