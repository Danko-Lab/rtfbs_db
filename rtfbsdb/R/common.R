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


write.starchbed <- function(bed, file.starch) {
	# pipe bed into starch file
	r <- try( write.table(bed, file = pipe(paste("sort-bed - | starch - >", file.starch)), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t") );

	if(class(r)=="try-error") 
		return("ERROR");

	if( !file.exists(file.starch) )    
	{
		cat("! Failed to write starch file (", file.starch, ") using the sort-bed and starch commands.\n");
		return("ERROR");
	}		
	else
		return(file.starch);
}

# get chromosome size based on genemo data of 2 bit file.
#
# @file.twoBit, e.g. hh19.2bit, mm10.2bit
#
get_chromosome_size <- function(file.twoBit)
{
	file.tmp <- tempfile();

	err_code <- system(paste("twoBitInfo", file.twoBit, file.tmp, sep=" "));
	chromInfo <- read.table( file.tmp );
	unlink( file.tmp );

	if( err_code!=0 ) cat("! Failed to call the twoBitInfo command\n");

	return(chromInfo);
}

# save the chromosme size into one ftable file
#
# save_chromosome_size("mm10.2bit", "chrom_info_mm19.tab")
#
save_chromosome_size <- function( file.twoBit, file.tab)
{
	write.table(get_chromosome_size(file.twoBit), file=file.tab, row.names=F, col.names=F, quote=F);
}


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

check_bed<-function( df.bed )
{
	if (NCOL(df.bed)<3 )
	{
		warning("At least 3 columns in BED file.");
		return (FALSE);
	}	

	if (any(is.na(df.bed[,c(1:3)]))) 
	{
		warning("NA values in first three columns in BED file.");
		return (FALSE);
	}	

	if (  !(class(df.bed[,3]) == "numeric" || class(df.bed[,3]) == "integer")
	   || !(class(df.bed[,2]) == "numeric" || class(df.bed[,2]) == "integer") ) 
	{
		warning("The 2nd column or 3rd column are not numeric in BED file.");
		return (FALSE);
	}	
	
	if ( any(df.bed[,2] > df.bed[,3])) 
	{
		warning("The 2nd column is greater than the 3rd column in BED file.");
		return (FALSE);
	}	
	
	if ( NCOL(df.bed)==6 )
	{
		if (!(unique(df.bed[,6]) %in% c("+", "-", ".")))
		{
			warning("Three values are avalaible for the 6th column in BED file, '+', '-' or '.'.");
			return (FALSE);
		}	
	}	
	
	return(TRUE);
}