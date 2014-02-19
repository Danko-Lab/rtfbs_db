## 
## Generates a tfbs object that stores
## PWMs for all motifs in R.
## 

setwd("/usr/data/GROseq.parser/pwm_data")

jolma <- paste("jolma/teal/", dir(path="jolma/teal/", pattern="pwm"), sep="")
jaspar<- paste("jaspar/", dir(path="jaspar/", pattern="pwm"), sep="")
neph  <- paste("neph/", dir(path="neph/", pattern="pwm"), sep="")

require(atif)
tfs <- tfbs(c(jolma, jaspar, neph), header=TRUE)
tfs <- tfbs.getDistanceMatrix(tfs, ncores=5)

save.image("tfbs.object.RData")

#############################################
##  Compute the optimal number of clusters.
##
load("tfbs.object.RData")
cluster_eval <- c(264,265) #seq(206, 299, 3) #seq(341,354,1) #seq(355,365,1) #seq(300,500,20) #seq(100, 1000, 100)
thresh <- 10

require(cluster)

d <- (1-tfs@distancematrix)^5
clu <- agnes(d, diss=TRUE)

###
## Minimize the within-cluster mean distance matrix

## Returns the minimum over clusters, of the mean within-cluster difference.
within_dif_max <- sapply(cluster_eval, function(i) {cn= cutree(clu, k=i); return(max(sapply(1:i, function(j) {mean(d[cn==j,cn==j])})))})
plot(c(cluster_eval), within_dif_max, type="b") ## Want to Minimize  this.
plot(c(cluster_eval), log(within_dif_max), type="b") ## Want to Minimize this.
## Smooth log-decay.

###
## Scan all DNAse-1 peaks in four cell types, with the set of all motifs.

k562_bed <- read.table("/usr/data/GROseq.parser/hg19/k562/dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz")
set.seed(200); k562_subset <- sample(1:NROW(k562_bed), 1000)
k562 <- scanDb_rtfbs(tfs, k562_bed[k562_subset,], ncores=10, threshold=thresh)

## Compuate a new difference object, with the fraction of sites that overlap between any two motifs.

system("mkdir ~/no-backup/tmp")
mean_overlap_between <- mclapply(cluster_eval, mc.cores=2, function(i) {
  ## Combine clusters into bed-like files.
  cn=cutree(clu, k=i)
  for(j in 1:NROW(k562)) {
    indx <- which(k562[[j]][,5] > thresh)
    write.table(k562[[j]][indx,], paste("~/no-backup/tmp/tmp.i", i, ".cn", cn[[j]], ".bed", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
  }
  
 n_bases <- integer(0)
 for(j in 1:i) {
  	system(paste("cat ~/no-backup/tmp/tmp.i", i, ".cn", j, ".bed | sort-bed - | bedops --merge - >  ~/no-backup/tmp/tmp.i", i, ".cn", j, ".sm.bed", sep=""))
	n_bases <- c(n_bases, sum(as.integer(tryCatch(system(paste("bedmap --bases ~/no-backup/tmp/tmp.i", i, ".cn", j, ".sm.bed", sep=""), intern=TRUE), error= function(x) {return(0)}))))
  }
  
  ## Now get the fraction.
  vals <- integer(0)
  for(j in 2:i) {
	  vals <- c(vals, sapply(1:(j-1), function(k) { 
	    	sum(as.integer(tryCatch( 
				system(paste("bedmap --bases-uniq ~/no-backup/tmp/tmp.i", i, ".cn", j, ".sm.bed ~/no-backup/tmp/tmp.i", i, ".cn", k, ".sm.bed", sep=""), intern=TRUE), 
				error= function(x) {return(0)}))) / (min(n_bases[j], n_bases[k])+1)
	  }))
  }
  
  ## Cleanup
  system(paste("rm ~/no-backup/tmp/tmp.i", i, ".*.bed", sep=""))
  
  ## And return...
  return(vals)
})
between_overlap_max <- sapply(mean_overlap_between, max)

save.image("nClusters.264-265.1.RData")

#############################################
##  Evaluate the optimal number of clusters.
##
## NOTE: Stoppped and started several times here.
#q("no")

norm <- function(x){
  return((x-min(x))/(max(x)-min(x)))
}

load("nClusters.264-265.1.RData")
between_overlap_max_264_265 <- sapply(mean_overlap_between, max)
within_dif_max_264_265 <- within_dif_max
cluster_eval_264_265 <- cluster_eval

load("nClusters.206-299.3.RData")
between_overlap_max_206_299 <- sapply(mean_overlap_between, max)
within_dif_max_206_299 <- within_dif_max
cluster_eval_206_299 <- cluster_eval

load("nClusters.341-354.1.RData")
between_overlap_max_341_354 <- sapply(mean_overlap_between, max)
within_dif_max_341_354 <- within_dif_max
cluster_eval_341_354 <- cluster_eval

load("nClusters.355-365.1.RData")
between_overlap_max_355_365 <- sapply(mean_overlap_between, max)
within_dif_max_355_365 <- within_dif_max
cluster_eval_355_365 <- cluster_eval

load("nClusters.300-500.20.RData")
between_overlap_max_300_500 <- sapply(mean_overlap_between, max)
within_dif_max_300_500 <- within_dif_max
cluster_eval_300_500 <- cluster_eval

load("nClusters.100-1000.100.RData")
between_overlap_max_100_1000 <- sapply(mean_overlap_between, max)
within_dif_max_100_1000 <- within_dif_max
cluster_eval_100_1000 <- cluster_eval

between_overlap_max <- c(between_overlap_max_100_1000, between_overlap_max_300_500, between_overlap_max_355_365, between_overlap_max_341_354, between_overlap_max_206_299, between_overlap_max_264_265)
within_dif_max <- c(within_dif_max_100_1000, within_dif_max_300_500, within_dif_max_355_365, within_dif_max_341_354, within_dif_max_206_299, within_dif_max_264_265)
cluster_eval <-  c(cluster_eval_100_1000, cluster_eval_300_500, cluster_eval_355_365, cluster_eval_341_354, cluster_eval_206_299, cluster_eval_264_265)

tab <- data.frame(cluster_eval, within_dif_max, between_overlap_max, scaled_product=norm(1/between_overlap_max) * norm(1/within_dif_max))
tab <- tab[order(cluster_eval),]
tab

pdf("clusterGraphEval.pdf")
 par(mar= c(5,5,2,5))
 plot(tab$cluster_eval, tab$within_dif_max, type="b", ylab="Within Cluster Max Diff", xlab="Number of Clusters") ## This subtracts the mean.
 par(new=TRUE)
 plot(tab$cluster_eval, tab$between_overlap_max, type="b", col="dark blue", axes=FALSE, xlab=NA, ylab=NA)
 axis(side = 4)
 mtext(side = 4, line = 3, "Between Cluster Overlap")
dev.off()

