## Reads in a motif from a PWM using the TRANSFAC format.
## example: 
## setwd("/usr/data/GROseq.parser/pwm_data/teal")
## read.motif("FOXL1.RTAAAYA.pwm", header=TRUE)
##
## For loading JASPAR, Use Andre's: /usr/projects/GROseq.parser/jaspar.load.R
##

read.motif <- function(motif_path, pseudocount= -7, force_even= FALSE, ...) 
{
	## Read the pwm and sequence file.
	motif <- tryCatch({ 
			as.matrix(read.table(motif_path, ...))
		}, 
		error= function(e) 
		{ 
			show(e); 
			#read.pwm(motif_path);
			return( NULL ) ;
		})
	
	if( is.null(motif) || NROW(motif) == 0 ) return( NULL );

	if(sum(motif[1,])>0) {
		motif <- log(motif/rowSums(motif)) # Divide by counts.
	}
	
	## Set a pseudocount of ~0.1%.
	motif[motif==-Inf] <- pseudocount;
	
	## Off-by-one bug for odd motifs.  Not sure why?!  For now hack it.
	if((NROW(motif) %% 2) == 1 & force_even) 
		motif <- rbind(motif, log(c(0.25, 0.25, 0.25, 0.25)));

	return(motif);
}

## Reverse complement.
## Assume motif columns in order: A, C, G, T
reverse.complement <- function(motif) 
{
	colnames(motif) <- colnames(motif)[c(4:1)]

	return( motif[c(NROW(motif):1),c(4:1)] ) ## Reverse complement.
}

## Correlation between motifs w/ equal size.
cor.motif.eq.size <- function(motif1, motif2) 
{
	stopifnot(NROW(motif1) == NROW(motif2))

	return(cor(as.vector(motif1), as.vector(motif2)))
}

## Extends a shorter motif w/ a pre-defined background (BG)= {A, C, G, T}.
extend.motif <- function(motif2, left, right, BG) 
{
	leftBG  <- t(matrix(rep(BG, left), nrow=4))
	rightBG <- t(matrix(rep(BG, right), nrow=4))

	return(rbind(leftBG, motif2, rightBG))
}

## Compares two motifs.  Returns Pearson's R for the log-values.
compare.motifs <- function(motif1, motif2, BG=log(c(0.25, 0.25, 0.25, 0.25))) 
{
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
	}
	else 
	{
		max_pos <- which.max(max_minus)
		maxStr_mot1 <- motif1_rc
	}

	## One we have the optimal placement, extend the motifs with background.  
	motif_score <- cor.motif.eq.size(maxStr_mot1, extend.motif(motif2, (max_pos-1), (ld-max_pos), BG))

	return(motif_score)
}