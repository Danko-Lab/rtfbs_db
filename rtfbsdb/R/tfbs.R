 
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
	}, error= function(e) { show(e); read.pwm(motif_path)})
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
    species    = "character",    ## such as, Homo_sapiens or Mus_musculus.
    TFID       = "character",    ## A non-unique ID for TF i.
    ntfs       = "integer",      ## Number of motifs in matrix.
    pwm        = "list",         ## PWM for TF i.
    extra_info = "data.frame",   ## extra finformation for PWMs, it maybe different with motif databse
    #usemotifs = "integer",      ## The indices of TFs to be used for analyses, such as scanning DNA sequences.
    filename   = "character",    ## The filename of the PWM.
    distancematrix="matrix",     ## Distance matrix between motifs
    #cluster    = "matrix",      ## The number of the cluster that this TF is included in.
    mgisymbols = "character",    ## Unique gene symbols for TF i.
    expressionlevel= "data.frame"   ## Expression level.
  ),
  )

## Creates a new tfbs object.  Reads files
tfbs <- function(filenames, names, species="Homo_sapiens", extra_info=NULL, ...) {
  pwms <- list()
  for(i in 1:length(filenames)) {
    curr <- read.motif(filenames[i], ...)
    pwms[[i]] <- curr
  }
  
  if( is.null(extra_info) )
     extra_info  <- as.data.frame(NULL);

  new("tfbs", 
    species    = species,
	TFID       = names, 
	ntfs       = as.integer(length(filenames)),
	extra_info = extra_info,
	filename   = filenames, 
	distancematrix= matrix(0, length(filenames), length(filenames)), 
	#cluster   = as.matrix(NULL), 
	#usemotifs = as.integer(1:length(filenames)),
	mgisymbols = as.character(names), 
	expressionlevel=as.data.frame(NULL), 
	pwm= pwms)
}

# Create a new tfbs object from all the PWM files found in the supplied folders.
# Optionally recursively descends into subfolders.
tfbs.dirs <- function(..., species="Homo_sapiens", args.read.motif = NULL, pattern = glob2rx("*.pwm"), recursive = FALSE) {
  pwms = list()
  paths = list(...)
  filenames = NULL
  pwm.names = NULL

  k = 1
  for (path in paths) {
    fnames = list.files(path, pattern = pattern, full.names = TRUE, recursive = recursive)

    for (fname in fnames) {
      pwm = do.call("read.motif", c(list(fname), args.read.motif))
      pwms[[k]] = pwm
      k = k + 1
    }

    filenames = c(filenames, fnames)
    pwm.names = c(pwm.names, sapply(fnames, function(str) {
      tmp = basename(str)
      substr(tmp, 1, nchar(tmp) - 4) # drop .pwm from name
    }))
  }

  names(pwm.names) <- NULL # clear filenames

  # build object instance
  new("tfbs", 
      TFID     = "", 
      species  = species,
      ntfs     = as.integer(length(filenames)),
      filename = filenames, 
      distancematrix = matrix(0, length(filenames), length(filenames)), 
      #cluster = as.matrix(NULL), 
      #usemotifs = as.integer(1:length(filenames)),
      mgisymbols = pwm.names,
      expressionlevel = as.data.frame(NULL), 
      pwm = pwms)
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
    def=function(tfbs, subset=NA, pdf.heatmap=NA, method=NA, group.k=NA) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.clusterMotifs")
	})

setMethod("tfbs.clusterMotifs", c(tfbs="tfbs"), tfbs_clusterMotifs);

#setMethod("tfbs.clusterMotifs", c(tfbs="tfbs"), 
#    function(tfbs, n_motifs, draw_heatmap= FALSE, subset=NULL) {
#      mat <- tfbs@distancematrix
#      if(!is.null(subset)) 
#         mat <- tfbs@distancematrix[subset,subset,drop=F];
#         
#      hc1 <- agnes(as.dist((1-mat)^5), diss=TRUE)
#      cuth <- cutree(hc1, k=n_motifs)
#
#      if(draw_heatmap) {
#        hc1 <- as.dendrogram(hc1)
#        ord.hc1 <- order.dendrogram(hc1)
#        hc2 <- reorder(hc1, mat[ord.hc1])
#        ord.hc2 <- order.dendrogram(hc2)
#
#        pal100 <- c("#7E291B","#66E52C","#8F66F0","#58DBE8","#396526","#EAABC1","#E1C33C","#3E3668","#EB3F90","#C3E6A8","#E74618","#66A2E9","#3E7774","#DF9056","#3C2C21","#DF40D7","#6CEF92","#8C5A6B","#BC8AE2","#A03B99","#56AC2D","#389C6C","#E26B7E","#706B4B","#D2E374","#A0A560","#7B1C3E","#49F7DB","#C8C6E9","#414FA2","#95A590","#8A669A","#98A62F","#9E792C","#D69489","#547FE5","#DF6340","#849BAB","#E63A61","#8A386D","#DAE338","#263715","#BBDEE3","#3F1324","#A1E03B","#383544","#76C2E8","#794D38","#DD74E2","#D7BF8E","#E366B9","#894D19","#D57221","#D9B764","#B0303D","#D6BFBC","#757A28","#D991CC","#356344","#E3A22D","#223C36","#83B664","#D8E4C5","#AE9ED9","#37576D","#CF7398","#6BEAB2","#6C9266","#B33662","#4B340E","#57E65B","#E23635","#9B464B","#757089","#578BBA","#A6311B","#B2714E","#457F20","#4EA3B1","#B93286","#73D463","#531914","#6B78C2","#E07467","#8D786E","#515018","#361E40","#AA4FCC","#90D2C4","#71469B","#419C4C","#37558A","#B393B0","#9BE090","#856BDA","#66B593","#47CA84","#5D205B","#672F3F","#59D8C2")
#        pal500 <- rep(pal100, 5)
#
#        #pl <- levelplot((mat)[ord.hc2, ord.hc2], col.regions= yb.sig.pal(100, scale=3), xlab="", ylab="",
#        print( levelplot((mat)[ord.hc2, ord.hc2], col.regions= yb.sig.pal(100, scale=3), xlab="", ylab="",
#        colorkey = list(space="left", labels=list(cex=1.5)), 
#        legend = list(
#          top = list(fun = dendrogramGrob,
#          args = list(x = hc2, ord = ord.hc2, side = "top", #lwd=2,
#          size = 7, size.add = 0.5, 
#          add = list(rect = list(col = "transparent", fill = pal500[cuth])),
#          type = "rectangle")))) );
#      }
#	  tfbs@cluster <- cbind(subset, cuth);
#	  return(tfbs);
#})


################################################
## Basic functions for drawing TFs.

## Draws the logo for a single tf.
setGeneric("tfbs.drawLogo", 
    def=function(tfbs, file.pdf=NULL, index=NULL, tf_id=NULL, motif_id=NULL, tf_name=NULL, family_name=NULL, tf_status=NULL, groupby=NULL) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.drawLogo")
	})
setMethod("tfbs.drawLogo", c(tfbs="tfbs"),
    function(tfbs, file.pdf=NULL, index=NULL, tf_id=NULL, motif_id=NULL, tf_name=NULL, family_name=NULL, tf_status=NULL, groupby=NULL) 
    {
      idx.select <- c();
      
	  if( !is.null(index) ) idx.select <- c(idx.select, index);
	  if( !is.null(tfbs@extra_info)  && !is.null(tf_id) ) idx.select <- c( idx.select, which(!is.na(match(tfbs@extra_info$TF_ID, tf_id ))) );
	  if( !is.null(tfbs@extra_info)  && !is.null(motif_id) ) idx.select <- c( idx.select, which(!is.na(match(tfbs@extra_info$Motif_ID, motif_id ))) );
	  if( !is.null(tfbs@extra_info)  && !is.null(tf_name) ) idx.select <- c( idx.select, which(!is.na(match(tfbs@extra_info$TF_Name, tf_name ))) );
	  if( !is.null(tfbs@extra_info)  && !is.null(family_name) ) idx.select <- c( idx.select, which(!is.na(match(tfbs@extra_info$Family_Name, family_name ))) );
	  if( !is.null(tfbs@extra_info)  && !is.null(tf_status) ) idx.select <- c( idx.select, which(!is.na(match(tfbs@extra_info$TF_Status, tf_status ))) );
	  if( is.null(index) && is.null(tf_id) && is.null(motif_id) && is.null(tf_name) && is.null(family_name) && is.null(tf_status)) idx.select <- c(1:tfbs@ntfs);
	  
	  if( length(which(is.na(idx.select)))>0 )
	  {
	  	  cat("!", length(which(is.na(idx.select))), "motifs can not be found in the tfbs object.");
	  	  idx.select <- idx.select[ !is.na(idx.select) ];
	  }
	  
	  idx.select <- sort(unique(idx.select));

	  draw_viewport <- function(i, xaxis = TRUE, yaxis = TRUE, cex = 1 )
	  {
		  tf_name <- tfbs@mgisymbols[i];	
		  if (!is.null(tfbs@extra_info))
			  tf_name <- paste( tfbs@extra_info[ i, "TF_Name"], " (Motif_ID:", tfbs@extra_info[ i, "Motif_ID"], "/DBID:", tfbs@extra_info[ i, "DBID"], ")", sep="");

		  pushViewport( viewport(x=0, y=0.94, width=1, height=0.05, just=c("left","bottom")) );
		  grid.text( tf_name, rot=0, gp=gpar(cex=cex), check.overlap=T);  
		  popViewport();

		  pushViewport( viewport(x=0, y=0, width=1, height=0.94, just=c("left","bottom")) );
		  seqLogo( exp(t(tfbs@pwm[[i]])), xaxis = xaxis, yaxis = yaxis)
		  popViewport();
	  
	  }
	
	  tfbs.grpupdby<-c("Family_Name", "TF_Name", "TF_Status", "Motif_Type");
	  if( !is.null(groupby) && length( which(groupby == tfbs.grpupdby ) ) == 0)
	  {
	  	 cat("The availabe value for groupby are Family_Name, TF_Name, TF_Status and Motif_Type.\n");
	  	 groupby <- NULL;
	  }	 

	  if( !is.null(file.pdf) && !is.na(file.pdf) ) pdf(file.pdf); 
	  
	  if(is.null(groupby))	
	  {
	  	  for(i.motif in idx.select)
	      {
	      	  if( i.motif<=0 || i.motif>tfbs@ntfs)
	      	      next;
      	  
			  grid.newpage()
			  vp1 <- viewport(x=0, y=0, width=1, height=1, just=c("left","bottom"))
			  pushViewport(vp1);
              draw_viewport(i.motif);
              popViewport();
      	  }
      }
      else
      {
      	  groups <- unique( as.character( tfbs@extra_info[idx.select, c(groupby)] ) )

	  	  for(k in 1:length(groups) )
	      {
	          idx.page <- idx.select[ which( tfbs@extra_info[idx.select, c(groupby)] == groups[k])];

			  grid.newpage();
			  vp1 <- viewport(x=0, y=0, width=1, height=1, just=c("left","bottom"))
			  pushViewport(vp1);

	  	      for(i in 1:length(idx.page) )
	          {
	              if(i>10 && i%%10==1)
	              {
	              	  popViewport();
	                  grid.newpage();
			          vp1 <- viewport(x=0, y=0, width=1, height=1, just=c("left","bottom"))
			          pushViewport(vp1);
                  }
                  
	              i.motif <- idx.page[i];
	      	      if( i.motif<=0 || i.motif>tfbs@ntfs )
	      	         next;
      	  
			      vp1 <- viewport(y = 0.8 - (((i-1)%%10)%/%2)/5, x=(i+1)%%2/2, width=0.5, height=0.2, just=c("left","bottom"))
			      pushViewport(vp1);
                  draw_viewport(i.motif, xaxis = FALSE, yaxis = FALSE, cex=0.6 );
			      popViewport();
			      
      	      }

			  popViewport();
          }
      }

	  if( !is.null(file.pdf) && !is.na(file.pdf) ) dev.off(); 
	
	})

## Draws all logos for each cluster.
setGeneric("tfbs.drawLogosForClusters", 
    def=function(tfbs, cluster.mat, file.pdf) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.drawLogosForClusters")
	})
setMethod("tfbs.drawLogosForClusters", c(tfbs="tfbs"),
    function(tfbs, cluster.mat, file.pdf=NULL) {

		if(  !is.null(file.pdf) && !is.na(file.pdf) )
	    	if( !check_folder_writable( file.pdf ) ) 
	  		    cat("! Can not create pdf file: ",file.pdf, "\n");

		if( !is.null(file.pdf) && !is.na(file.pdf) )
		{
			r.try <- try ( pdf( file.pdf ) );	
			if(class(r.try)=="try-error")
			{
				cat("! Failed to create PDF file for motif logos.\n");
				return;
			}	
		}
		
		draw_tf_name <- function(tfbs, idx, nmax.motifs)
		{
		  	tf_name <- tfbs@mgisymbols[idx];	
		  	if (!is.null(tfbs@extra_info))
		  		tf_name <- tfbs@extra_info[ idx, "TF_Name"];
			
			n.block  <- round(nmax.motifs/2);
		  	pushViewport( viewport(x=0.01, y=0, width=0.05, height=1, just=c("left","bottom")) );
		  	
	  		str.inches <- strwidth(tf_name, units = 'in');
	  	    y.scale <- convertHeight(unit(str.inches, "inches"), "native", valueOnly=T );
		
			cex <- 1;
			if( y.scale> 0.8 ) cex <- 0.8/y.scale;
			if( cex > 1 ) cex <- 1;

	  		grid.text( tf_name, rot=90, gp=gpar(cex=cex), check.overlap=T);  

		  	popViewport();
		}

		for(i in 1:max(cluster.mat[,2])) 
		{
			motifs_in_cluster <- cluster.mat[which(cluster.mat[,2]==i), 1];
			nmotifs <- length(motifs_in_cluster);
			
			nmax.motifs <- nmotifs;
			# at least4 motif space will be designed in each page.
			if( nmax.motifs <= 4 ) nmax.motifs <- 4;
			
			grid.newpage();
			for(j in seq(1, nmax.motifs, 2)) {
			
			  vp1 <- viewport(x=0, y=((j-1)/2)/ceiling(nmax.motifs/2), width=0.5, height=1/ceiling(nmax.motifs/2), just=c("left","bottom"));
			  pushViewport(vp1); 
			  #seqLogo(makePWM(exp(t(tfbs@pwm[[motifs_in_cluster[j]]]))), xaxis = FALSE, yaxis = FALSE);
			  if ( j <= nmotifs )
			  {
			  	draw_tf_name(tfbs, motifs_in_cluster[j], nmax.motifs);
		  	    pushViewport( viewport(x=0.04, y=0, width=0.96, height=1, just=c("left","bottom")) );
			  	seqLogo(exp(t(tfbs@pwm[[motifs_in_cluster[j]]])), xaxis = FALSE, yaxis = FALSE);
			  	popViewport();
			  }
			  popViewport();
			  
			  vp1 <- viewport(x=0.5, y=((j-1)/2)/ceiling(nmax.motifs/2), width=0.5, height=1/ceiling(nmax.motifs/2), just=c("left","bottom"));
			  pushViewport(vp1); 
			  #seqLogo(makePWM(exp(t(tfbs@pwm[[motifs_in_cluster[j]]]))), xaxis = FALSE, yaxis = FALSE);
			  if ( j+1 <= nmotifs )
			  {
			  	draw_tf_name(tfbs, motifs_in_cluster[j+1], nmax.motifs);
		  	    pushViewport( viewport(x=0.04, y=0, width=0.96, height=1, just=c("left","bottom")) );
			  	seqLogo( exp(t(tfbs@pwm[[motifs_in_cluster[j+1]]])), xaxis = FALSE, yaxis = FALSE);
			  	popViewport();
			  }
			  popViewport();
			}
		}

		if( !is.null(file.pdf) && !is.na(file.pdf) )
			dev.off();	
	})

#
# The following function has been moved to express.R 
#
##################################################
## Gets expression level of target TF.
## TODO: Add the MGI symbol to each TF.  Not 100% sure where to do this?!
setGeneric("tfbs.getExpression", 
    def=function(tfbs, file.bigwig.plus, file.bigwig.minus, file.bam=NA, file.twoBit=NA, file.gencode.gtf=NA, seq.datatype=NA, ncores = 3) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.getExpression")
	})

setMethod("tfbs.getExpression", c(tfbs="tfbs"), tfbs_getExpression );

#    function(tfbs, bed, file_plus, file_minus) {
#      bw_plus <- load.bigWig("file_plus")
#      bw_minus <- load.bigWig("file_minus")
#	  
#      indx <- match(tfbs@mgisymbols, bed[,4])
#      bw_plus <- load.bigWig(file_plus)
#      bw_minus <- load.bigWig(file_minus)
#      tfbs@expressionlevel <- bedQuery.bigWig(bed[indx,], bw_plus, bw_minus)
#      return(tfbs)
#})
####################################################

## find TF sites in the BED range from sequence data file(hg19/hg19.2bit);
## see codes in scan_sites.R
##
setGeneric("tfbs.scanTFsite", 
    def=function( tfbs, file.twoBit, tre.bed = NULL, return.type=c("matches", "posteriors", "maxposterior", "writedb"), file.prefix = NA,  usemotifs = NA, ncores = 3,  
                  fdr = NA, threshold = 6, gc.groups = NA, background.order = 2, background.length = 100000 ) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.scanTFsite")
	})

setMethod("tfbs.scanTFsite", c(tfbs="tfbs"), tfbs_scanTFsite );


## Comparative TFBS search between positive BED and negative BED
## see codes in comp_sites.R
##
setGeneric("tfbs.compareTFsite", 
    def=function( tfbs, file.twoBit, positive.bed, negative.bed, file.prefix = NA, use.cluster = NA, ncores = 3,
   	              gc.correction = TRUE, gc.correction.pdf=NA, fdr = NA, threshold = 6, gc.groups=1, background.order = 2, background.length = 100000, pv.adj=p.adjust.methods) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.compareTFsite")
	})

setMethod("tfbs.compareTFsite", c(tfbs="tfbs"), tfbs_compareTFsite );

setGeneric("tfbs.selectByGeneExp", 
    def=function(tfbs, cluster.mat) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.selectByGeneExp")
	})
setMethod("tfbs.selectByGeneExp", c(tfbs="tfbs"),
    function(tfbs, cluster.mat) {
      stopifnot(!is.null(tfbs@expressionlevel) );
	
	  cluster <- cluster.mat[,2]

      usemotifs <- sapply(1:max(cluster), function(x) {
	    a <- which(cluster == x)
		if(length(a) > 1) {
		  a.min <- which.min( tfbs@expressionlevel$prob[a] );
		  if(length(a.min)>0)
		  	return( a[ a.min[1] ] ) 
		  else
		    return( sample(a, 1) );
		} else {
	      return(a)
		}})
		
      return(cluster.mat[ usemotifs,1 ]);
})

setGeneric("tfbs.selectByRandom", 
    def=function(tfbs, cluster.mat) {
	  stopifnot(class(tfbs) == "tfbs")
	  standardGeneric("tfbs.selectByRandom")
	})
	
setMethod("tfbs.selectByRandom", c(tfbs="tfbs"),
    function(tfbs, cluster.mat) {
	  cluster <- cluster.mat[,2]

      usemotifs <- sapply(1:max(cluster), function(x) {
	    a <- which(cluster == x)
		if(length(a) > 1) {
		  return(sample(a, 1)) ## DANGEROUS!! If length(a) == 1, samples from 1:a[1].
		} else {
	      return(a)
		}})
		
      return(cluster.mat[ usemotifs, 1])
})

setMethod("show", "tfbs", function(object){

	cat("Species: ", object@species, "\n");
	cat("TF number: ", object@ntfs, "\n");

	if(all(which(object@distancematrix==0))) 
		cat( "Distance Matrix:  NULL\n" )
	else
		cat( "Distance Matrix:  [", NROW(object@distancematrix), ",", NCOL(object@distancematrix), "]\n" );
	
	if(NROW(object@expressionlevel)==0) 
		cat( "Expression:  NULL\n" )
	else
		cat( "Expression:  [", NROW(object@expressionlevel), ",", NCOL(object@expressionlevel), "]\n" );

	df <- NULL;	
	if(!is.null(object@extra_info))
	{
		df <- object@extra_info[,c("Motif_ID", "DBID", "TF_Name", "Family_Name", "Motif_Type", "MSource_Identifier")]
		df <- data.frame(df, filename=basename(object@filename));
	}
	else
		df <- data.frame(Motif_ID=object@mgisymbols, filename=basename(object@filename));

	if(NROW(object@expressionlevel)>0) 
		df <- data.frame(df, p.pois = object@expressionlevel[,c("p.pois")] );
	
	cat("\nPartial list of TFs\n");
	show(head(df, 20));
});


check_folder_writable<-function(file.prefix)
{
  if( dirname(file.prefix) != "." )
  	 dir.create( dirname(file.prefix), showWarnings = TRUE, recursive = FALSE );
  
  file.temp <- tempfile( tmpdir = dirname(file.prefix), fileext = "")
  r.try <- try( file.create(file.temp) )	
  if( class(r.try)=="try-error" || r.try==FALSE )
  	return(FALSE);
  	
  unlink(file.temp);
  return(TRUE);
}


