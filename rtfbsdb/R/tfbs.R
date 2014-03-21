 
## Reads in a motif from a PWM using the TRANSFAC format.
## example: 
## setwd("/usr/data/GROseq.parser/pwm_data/teal")
## read.motif("FOXL1.RTAAAYA.pwm", header=TRUE)
##
## For loading JASPAR, Use Andre's: /usr/projects/GROseq.parser/jaspar.load.R
##
read.motif <- function(motif_path, pseudocount= -7, force_even= FALSE, ...) {
  ## Read the pwm and sequence file.
  motif <- tryCatch({ 
	as.matrix(read.table(motif_path, ...))
	#motif <- motif/rowSums(motif) 
	}, error= function(e) {read.pwm(motif_path)})
  if(sum(motif[1,])>0) {
	motif <- log(motif/rowSums(motif)) # Divide by counts.
  }
  motif[motif==-Inf] <- pseudocount ## Set a pseudocount of ~0.1%.
  if((NROW(motif) %% 2) == 1 & force_even) motif <- rbind(motif, log(c(0.25, 0.25, 0.25, 0.25)))  ## Off-by-one bug for odd motifs.  Not sure why?!  For now hack it.
 
  return(motif)
}

## Reverse complement.
## Assume motif columns in order: A, C, G, T
reverse.complement <- function(motif) {
  colnames(motif) <- colnames(motif)[c(4:1)]
  return( motif[c(NROW(motif):1),c(4:1)] ) ## Reverse complement.
}

## Correlation between motifs w/ equal size.
cor.motif.eq.size <- function(motif1, motif2) {
  stopifnot(NROW(motif1) == NROW(motif2))
  return(cor(as.vector(motif1), as.vector(motif2)))
}

## Extends a shorter motif w/ a pre-defined background (BG)= {A, C, G, T}.
extend.motif <- function(motif2, left, right, BG) {
  leftBG  <- t(matrix(rep(BG, left), nrow=4))
  rightBG <- t(matrix(rep(BG, right), nrow=4))
  return(rbind(leftBG, motif2, rightBG))
}

## Compares two motifs.  Returns Pearson's R for the log-values.
compare.motifs <- function(motif1, motif2, BG=log(c(0.25, 0.25, 0.25, 0.25))) {
  ## REQUIRE: NROW(motif1) >= NROW(motif2)
  if(NROW(motif1) < NROW(motif2)) {
    tmp <- motif1
	motif1 <- motif2
	motif2 <- tmp
	remove(tmp)
  }

  ## Create two vectors
  ld <- (NROW(motif1)-NROW(motif2)+1)
  size_small <- NROW(motif2)
  
  ## Compare natural.
  max_plus <- sapply(c(1:ld), function(x){ cor.motif.eq.size(motif1[c(x:(x+size_small-1)),], motif2) })
  
  ## Also reverse complement.
  motif1_rc <- reverse.complement(motif1)
  max_minus <- sapply(c(1:ld), function(x){ cor.motif.eq.size(motif1_rc[c(x:(x+size_small-1)),], motif2) })
  
  ## Get the maximum value.
  if(max(max_plus) >= max(max_minus)) {
    max_pos <- which.max(max_plus)
	maxStr_mot1 <- motif1
  }  else {
    max_pos <- which.max(max_minus)
	maxStr_mot1 <- motif1_rc
  }
  
  ## One we have the optimal placement, extend the motifs with background.  
  motif_score <- cor.motif.eq.size(maxStr_mot1, extend.motif(motif2, (max_pos-1), (ld-max_pos), BG))
  return(motif_score)
}

## Write.seqfile.
## Writes a fasta file containing all sequences in a bed region.
read.seqfile.from.bed <- function(bed, twoBitPath, tmpdir = getwd()) {
  # create temporary filenames
  tmp.seq = tempfile(tmpdir=tmpdir)
  tmp.fa = tempfile(tmpdir=tmpdir)
  
  # write sequence list
  seqList = paste(bed[,1],":", as.integer(bed[,2]), "-", as.integer(bed[,3]), sep="")
  writeLines(seqList, tmp.seq)
  
  # generate fasta file
  cmd = paste("twoBitToFa -seqList=", tmp.seq, " ", twoBitPath, " ", tmp.fa, sep="")
  system(cmd, wait = TRUE)
  
  # read data
  ms_data <- read.ms(tmp.fa)
  
  # clean up
  unlink(c(tmp.seq, tmp.fa))
  
  return(ms_data)
}

## For drawing a levelplot.
yb.sig.pal <- function(n, scale=10) {
 ints<- c(0:(n-1))/(n-1)   ## Linear scale from 0:1 x N values.
 ints<- 1/(1+exp(scale*(0.5-ints)))## Transfer to sigmoidal scale.
 b<- min(ints)
 m<- 2*b/(n-1)
 ints<- ints+(m*c(0:(n-1)) -b)## Transform by linear function to fill colors out to maxes.
 
 ## Transfer to colorspace.
 # Yellow: 255, 255, 0
 # White:  255, 255, 255
 # Blue:   0, 0, 255
 YW <- ints[ints < 0.5] *2
 WB <- (ints[ints >= 0.5]-0.5) *2
 YW[YW<0] <- 0; WB[WB>1] <- 1
 c(rgb(1, 1, YW), rgb(1-WB, 1-WB, 1))
}


##
## S4 class for storing TFBS.  
## sequence preferences of a set of TFs.
setClass("tfbs", 
  representation(
    TFID=     "character",      ## A non-unique ID for TF i.
    ntfs=     "integer",        ## Number of motifs in matrix.
    pwm=        "list",         ## PWM for TF i.
    usemotifs=  "integer",      ## The indices of TFs to be used for analyses, such as scanning DNA sequences.

    filename= "character",      ## The filename of the PWM.
    distancematrix="matrix",    ## Distance matrix between motifs
    cluster=  "integer",        ## The number of the cluster that this TF is included in.
    mgisymbols= "character",    ## Unique gene symbols for TF i.
    expressionlevel= "numeric"  ## Expression level.
  ),
  )


## Creates a new tfbs object.  Reads files
tfbs <- function(filenames, names, ...) {
  pwms <- list()
  for(i in 1:length(filenames)) {
    curr <- read.motif(filenames[i], ...)
    pwms[[i]] <- curr
  }
  new("tfbs", 
	TFID= "", 
	ntfs= as.integer(length(filenames)),
	filename= filenames, 
	distancematrix= matrix(0, length(filenames), length(filenames)), 
	cluster= as.integer(0), 
	usemotifs= as.integer(1:length(filenames)),
	mgisymbols= as.character(names), 
	expressionlevel=numeric(0), 
	pwm= pwms)
}

setGeneric("tfbs.getDistanceMatrix", 
    def=function(tfbs, ncores=3, BG=log(c(0.25, 0.25, 0.25, 0.25))) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.getDistanceMatrix")
	})
setMethod("tfbs.getDistanceMatrix", c(tfbs="tfbs"),
    function(tfbs, ncores=3, BG=log(c(0.25, 0.25, 0.25, 0.25))) {
      for(i in 1:tfbs@ntfs) {
		tfbs@distancematrix[i, ] <- unlist(mclapply(1:tfbs@ntfs, function(j) {compare.motifs(tfbs@pwm[[i]], tfbs@pwm[[j]], BG=BG)}, mc.cores= ncores))
#        for(j in 1:tfbs@ntfs) {
#          motif_cor <- compare.motifs(tfbs@pwm[[i]], tfbs@pwm[[j]])
#          tfbs@distancematrix[i, j] <- motif_cor
#        }
      }
      return(tfbs)
    })

## Clusters TFs based on DNA sequence preferences.
setGeneric("tfbs.clusterMotifs", 
    def=function(tfbs, n_motifs, draw_heatmap= FALSE) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.clusterMotifs")
	})
setMethod("tfbs.clusterMotifs", c(tfbs="tfbs"),
    function(tfbs, n_motifs, draw_heatmap= FALSE) {
      mat <- tfbs@distancematrix
      hc1 <- agnes(as.dist((1-mat)^5), diss=TRUE)
      cuth <- cutree(hc1, k=n_motifs)

      if(draw_heatmap) {
        hc1 <- as.dendrogram(hc1)
        ord.hc1 <- order.dendrogram(hc1)
        hc2 <- reorder(hc1, mat[ord.hc1])
        ord.hc2 <- order.dendrogram(hc2)

        pal100 <- c("#7E291B","#66E52C","#8F66F0","#58DBE8","#396526","#EAABC1","#E1C33C","#3E3668","#EB3F90","#C3E6A8","#E74618","#66A2E9","#3E7774","#DF9056","#3C2C21","#DF40D7","#6CEF92","#8C5A6B","#BC8AE2","#A03B99","#56AC2D","#389C6C","#E26B7E","#706B4B","#D2E374","#A0A560","#7B1C3E","#49F7DB","#C8C6E9","#414FA2","#95A590","#8A669A","#98A62F","#9E792C","#D69489","#547FE5","#DF6340","#849BAB","#E63A61","#8A386D","#DAE338","#263715","#BBDEE3","#3F1324","#A1E03B","#383544","#76C2E8","#794D38","#DD74E2","#D7BF8E","#E366B9","#894D19","#D57221","#D9B764","#B0303D","#D6BFBC","#757A28","#D991CC","#356344","#E3A22D","#223C36","#83B664","#D8E4C5","#AE9ED9","#37576D","#CF7398","#6BEAB2","#6C9266","#B33662","#4B340E","#57E65B","#E23635","#9B464B","#757089","#578BBA","#A6311B","#B2714E","#457F20","#4EA3B1","#B93286","#73D463","#531914","#6B78C2","#E07467","#8D786E","#515018","#361E40","#AA4FCC","#90D2C4","#71469B","#419C4C","#37558A","#B393B0","#9BE090","#856BDA","#66B593","#47CA84","#5D205B","#672F3F","#59D8C2")
        pal500 <- rep(pal100, 5)

        pl <- levelplot((mat)[ord.hc2, ord.hc2], col.regions= yb.sig.pal(100, scale=3), xlab="", ylab="",
        colorkey = list(space="left", labels=list(cex=1.5)), 
        legend = list(
          top = list(fun = dendrogramGrob,
          args = list(x = hc2, ord = ord.hc2, side = "top", #lwd=2,
          size = 7, size.add = 0.5, 
          add = list(rect = list(col = "transparent", fill = pal500[cuth])),
          type = "rectangle"))))
      }
	  tfbs@cluster <- cuth
      return(tfbs)

})

setGeneric("tfbs.setUseMotifs.random", 
    def=function(tfbs) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.setUseMotifs.random")
	})
setMethod("tfbs.setUseMotifs.random", c(tfbs="tfbs"),
    function(tfbs) {
      tfbs@usemotifs <- sapply(1:max(tfbs@cluster), function(x) {
	    a <- which(tfbs@cluster == x)
		if(length(a) > 1) {
		  return(sample(a, 1)) ## DANGEROUS!! If length(a) == 1, samples from 1:a[1].
		} else {
	      return(a)
		}})
		
      return(tfbs)
})


## Gets expression level of target TF.
## TODO: Add the MGI symbol to each TF.  Not 100% sure where to do this?!
setGeneric("tfbs.getExpression", 
    def=function(tfbs, bed, file_plus, file_minus) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.getExpression")
	})
setMethod("tfbs.getExpression", c(tfbs="tfbs"),
    function(tfbs, bed, file_plus, file_minus) {
      bw_plus <- load.bigWig("file_plus")
      bw_minus <- load.bigWig("file_minus")
	  
      indx <- match(tfbs@mgisymbols, bed[,4])
      bw_plus <- load.bigWig(file_plus)
      bw_minus <- load.bigWig(file_minus)
      tfbs@expressionlevel <- bedQuery.bigWig(bed[indx,], bw_plus, bw_minus)
      return(tfbs)
})


################################################
## Basic functions for drawing TFs.

## Draws the logo for a single tf.
setGeneric("tfbs.drawLogo", 
    def=function(tfbs, i) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.drawLogo")
	})
setMethod("tfbs.drawLogo", c(tfbs="tfbs"),
    function(tfbs, i) {
      grid.newpage()
      vp1 <- viewport(x=0, y=0, width=1, height=1, just=c("left","bottom"))
      pushViewport(vp1) 
      seqLogo(makePWM(exp(t(tfbs@pwm[[i]]))), xaxis = TRUE, yaxis = TRUE)
      popViewport()
	})

## Draws all logos for each cluster.
setGeneric("tfbs.drawLogosForClusters", 
    def=function(tfbs) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.drawLogosForClusters")
	})
setMethod("tfbs.drawLogosForClusters", c(tfbs="tfbs"),
    function(tfbs) {
		clusters <- tfbs@cluster
		
		for(i in 1:max(clusters)) {
			motifs_in_cluster <- which(clusters==i)
			nmotifs <- length(motifs_in_cluster)
			#print(nmotifs)
			
			grid.newpage()
			for(j in seq(1, nmotifs, 2)) {
			  vp1 <- viewport(x=0, y=((j-1)/2)/ceiling(nmotifs/2), width=0.5, height=1/ceiling(nmotifs/2), just=c("left","bottom"))
			  pushViewport(vp1) 
			  seqLogo(makePWM(exp(t(tfbs@pwm[[motifs_in_cluster[j]]]))), xaxis = FALSE, yaxis = FALSE)
			  popViewport()
			  
			  if((j+1) <= nmotifs) {
				vp1 <- viewport(x=0.5, y=((j-1)/2)/ceiling(nmotifs/2), width=0.5, height=1/ceiling(nmotifs/2), just=c("left","bottom"))
				pushViewport(vp1) 
				seqLogo(makePWM(exp(t(tfbs@pwm[[motifs_in_cluster[j]]]))), xaxis = FALSE, yaxis = FALSE)
				popViewport()
			  }
			}
		}
	})
