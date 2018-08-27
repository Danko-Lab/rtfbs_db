#' S4 class for storing TFBS.
#' sequence preferences of a set of TFs.

setClass("tfbs",
	representation(
		species        = "character",    ## such as, Homo_sapiens or Mus_musculus.
		ntfs           = "integer",      ## Number of motifs in matrix.
		mgisymbols     = "character",    ## Unique gene symbols for TF i.
		filename       = "character",    ## The filename of the PWM.
		pwm            = "list",         ## PWM for TF i.
		tf_info        = "data.frame",   ## TF information for PWMs, it maybe different with motif database
		tf_missing     = "data.frame",   ## Missing TF information for PWMs, it maybe different with motif database
		distancematrix = "matrix",       ## Distance matrix between motifs
		cluster        = "matrix",       ## The number of the cluster that this TF is included in.
		expressionlevel= "data.frame"    ## Expression level.
		#TFID           = "character",   ## A non-unique ID for TF i.
		#usemotifs     = "integer",      ## The indices of TFs to be used for analyses, such as scanning DNA sequences.
  ),
)

setMethod("show", "tfbs", function(object)
{
	cat("Species: ", object@species, "\n");
	cat("TF number: ", object@ntfs, "\n");

	if( NROW(object@distancematrix)==0 )
		cat( "Distance Matrix:  NULL\n" )
	else
		cat( "Distance Matrix:  [", NROW(object@distancematrix), ",", NCOL(object@distancematrix), "]\n" );

	if( NROW(object@cluster)==0 )
		cat( "Cluster Matrix:  NULL\n" )
	else
		cat( "Cluster Matrix:  [", NROW(object@cluster), ",", NCOL(object@cluster), "]\n" );

	if(NROW(object@expressionlevel)==0)
		cat( "Expression:  NULL\n" )
	else
		cat( "Expression:  [", NROW(object@expressionlevel), ",", NCOL(object@expressionlevel), "]\n" );

	df <- NULL;
	if( NROW(object@tf_info) >0 )
	{
		df.list <- lapply ( c("Motif_ID", "TF_Name", "Family_Name", "DBID", "Motif_Type"), function(x){
					if(x %in% colnames(object@tf_info) ) return(as.character(object@tf_info[,x])) else NA; } );

		df1 <- do.call( cbind, df.list );
		colnames(df1) <- c("Motif_ID", "TF_Name", "Family_Name", "DBID", "Motif_Type");

		df <- data.frame(df1, filename=basename(object@filename), stringsAsFactors=FALSE);
	}
	else
		df <- data.frame(Motif_ID=object@mgisymbols, filename=basename(object@filename), stringsAsFactors=FALSE);

	if( NROW(object@expressionlevel) > 0 )
		df <- data.frame(df, p.pois = object@expressionlevel[,c("p.pois")] );

	cat("\nPartial list of TFs\n");
	show(head(df, 20));
})

setGeneric("tfbs.importMotifs",
	def = function(tfbs, format, filenames, motif_ids=NULL, PPM.format=TRUE, skip.lines=0, pseudocount = -7, force_even = FALSE, ...){
		stopifnot(class(tfbs) == "tfbs");
		standardGeneric("tfbs.importMotifs");
	})

setMethod("tfbs.importMotifs", signature(tfbs="tfbs"), tfbs_importMotifs)


## Gets expression level of target TF.
## TODO: Add the MGI symbol to each TF.  Not 100% sure where to do this?!

setGeneric("tfbs.selectExpressedMotifs",
	def = function( tfbs,
			file.twoBit,
			file.gencode.gtf,
			file.bigwig.plus = NA,
			file.bigwig.minus = NA,
			file.bam = NA,
			seq.datatype=c( "GRO-seq", "PRO-seq", "RNA-seq" ),
			pvalue.threshold = 0.05,
			lowest.reads.RPKM = NA,
			include.DBID.missing = TRUE,
			use.strand=FALSE,
			ncores = 1 )
	{
		stopifnot(class(tfbs) == "tfbs");
		standardGeneric("tfbs.selectExpressedMotifs");
	})

setMethod("tfbs.selectExpressedMotifs", c(tfbs="tfbs"), tfbs_selectExpressedMotifs )

## Gets expression level of target TF.
## TODO: Add the MGI symbol to each TF.  Not 100% sure where to do this?!

setGeneric("tfbs.getExpression",
	def = function( tfbs,
			file.twoBit,
			file.gencode.gtf,
			file.bigwig.plus = NA,
			file.bigwig.minus = NA,
			file.bam = NA,
			seq.datatype=c( "GRO-seq", "PRO-seq", "RNA-seq" ),
			use.strand=FALSE,
			ncores = 1)
	{
		stopifnot(class(tfbs) == "tfbs")
		standardGeneric("tfbs.getExpression")
	})

setMethod("tfbs.getExpression", c(tfbs="tfbs"), tfbs_getExpression );

## Clusters TFs based on DNA sequence preferences.
setGeneric("tfbs.clusterMotifs",
	def = function(tfbs,
			method = c("agnes", "apcluster"),
			pdf.heatmap = NA,
			group.k = NA,
			apcluster.q = 0.99,
			ncores = 1,
			plot.style=c("rtfbsdb", "apcluster"),
			BG=log( c(0.25, 0.25, 0.25, 0.25)),
			... )
	{
		stopifnot(class(tfbs) == "tfbs")
		standardGeneric("tfbs.clusterMotifs")
	})

setMethod("tfbs.clusterMotifs", c(tfbs="tfbs"), tfbs_clusterMotifs);


## Draws the logo for a single tf.
setGeneric("tfbs.drawLogo",
	def = function(tfbs,
			file.pdf = NULL,
			index = NULL,
			tf_id = NULL,
			motif_id = NULL,
			tf_name = NULL,
			family_name = NULL,
			tf_status = NULL,
			groupby = NULL)
	{
		stopifnot(class(tfbs) == "tfbs")
		standardGeneric("tfbs.drawLogo")
	})

setMethod("tfbs.drawLogo", c(tfbs="tfbs"), tfbs_drawLogo );

## Draws all logos for each cluster.
setGeneric("tfbs.drawLogosForClusters",
	def = function(tfbs, file.pdf = NULL, nrow.per.page = 6, vertical=TRUE ) {
		stopifnot(class(tfbs) == "tfbs")
		standardGeneric("tfbs.drawLogosForClusters")
	})

setMethod("tfbs.drawLogosForClusters", c(tfbs="tfbs"), tfbs_drawLogosForClusters)


## find TF sites in the BED range from sequence data file(hg19/hg19.2bit);
## see codes in scan_sites.R
##
setGeneric("tfbs.scanTFsite",
	def = function( tfbs,
					file.genome,
					gen.bed = NULL,
					return.type = c("matches", "maxscore", "posteriors", "maxposterior", "writedb"),
					file.prefix = NA,
					usemotifs = NA,
					ncores = 1,
					threshold = 6,
					threshold.type = c("score", "fdr"),
					gc.groups = NA,
					background.order = 2,
					background.length = 100000,
					exclude_offset = 250,
					exclude_chromosome="_|chrM|chrY|chrX" )
	{
		stopifnot(class(tfbs) == "tfbs")
		standardGeneric("tfbs.scanTFsite")
	})

setMethod("tfbs.scanTFsite", c(tfbs="tfbs"), tfbs_scanTFsite )

## Comparative TFBS enrichment between positive BED and negative BED
## see codes in comp_sites.R
##
setGeneric("tfbs.enrichmentTest",
	def = function( tfbs,
					file.genome,
					positive.bed,
					negative.bed = NA,
					file.prefix  = NA,
					use.cluster  = FALSE,
					ncores = 1,
					gc.correction = TRUE,
					gc.correction.pdf = NA,
					gc.min.sample = 500,
					gc.robust.rep = NA,
					threshold = 6,
					threshold.type = c("score", "fdr"),
					gc.groups=1,
					background.order = 2,
					background.length = 100000,
					pv.adj=p.adjust.methods)
	{
		stopifnot(class(tfbs) == "tfbs");
		standardGeneric("tfbs.enrichmentTest");
	})

setMethod("tfbs.enrichmentTest", c(tfbs="tfbs"), tfbs_enrichmentTest)

setGeneric("tfbs.selectByGeneExp",
	def = function( tfbs ) {
		stopifnot(class(tfbs) == "tfbs")
		standardGeneric("tfbs.selectByGeneExp")
	})

setMethod("tfbs.selectByGeneExp", c(tfbs="tfbs"),
	function( tfbs ) {

		stopifnot( NROW(tfbs@expressionlevel)>0 && NROW(tfbs@cluster)>0 );

		cluster <- tfbs@cluster[,2]

		usemotifs <- sapply(1:max(cluster), function(x) {
			a <- which(cluster == x)
			if(length(a) > 1) {
				a.min <- which.min( tfbs@expressionlevel$p.pois[a] );
			if(length(a.min)>0)
				return( a[ a.min[1] ] )
			else
				return( sample(a, 1) );
			} else {
				return(a)
		}})

		select.col <- rep( 0, NROW(tfbs@cluster) );
		select.col[ usemotifs ]<-1;
		if(NCOL(tfbs@cluster)==2)
			tfbs@cluster <- cbind( tfbs@cluster, selected=select.col)
		else
			tfbs@cluster[,3] <- select.col;


		return( tfbs );
	})

setGeneric("tfbs.selectByRandom",
	def=function( tfbs ) {
		stopifnot(class(tfbs) == "tfbs")
		standardGeneric("tfbs.selectByRandom")
	})


setMethod("tfbs.selectByRandom", c(tfbs="tfbs"),
	function( tfbs ) {

		stopifnot( NROW(tfbs@cluster)>0 );

		cluster <- tfbs@cluster[,2]

		usemotifs <- sapply(1:max(cluster), function(x) {
			a <- which(cluster == x)
			if(length(a) > 1) {
				return(sample(a, 1)) ## DANGEROUS!! If length(a) == 1, samples from 1:a[1].
			} else {
				return(a)
		}})

		select.col <- rep( 0, NROW(tfbs@cluster) );
		select.col[ usemotifs ]<-1;

		if(NCOL(tfbs@cluster)==2)
			tfbs@cluster <- cbind( tfbs@cluster, selected=select.col )
		else
			tfbs@cluster[,3] <- select.col;

		return( tfbs );
	})
