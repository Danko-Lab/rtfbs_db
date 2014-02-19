## 
## Generates a tfbs object that stores
## PWMs for all motifs in R.
## 

thresh=6

require(atif)
setwd("/usr/data/GROseq.parser/pwm_data")

##### 263 #

#load("TF.db.263.clusters.RData")
#cd4_bed <- read.table("/usr/projects/GROseq/NHP/tss_caller/HCM-U-PI.dREG-tss-clusters.tsv")[,c(1:3)]
#cd4 <-  scanDb_rtfbs(tfs, cd4_bed, ncores=10, threshold=thresh, return_type="writedb", file_prefix=paste("cd4.dreg.tfs263.th", thresh, sep=""))

##### 341 #

#load("TF.db.341.clusters.RData")
#cd4_bed <- read.table("/usr/projects/GROseq/NHP/tss_caller/HCM-U-PI.dREG-tss-clusters.tsv")[,c(1:3)]
#cd4 <-  scanDb_rtfbs(tfs, cd4_bed, ncores=10, threshold=thresh, return_type="writedb", file_prefix=paste("cd4.dreg.tfs341.th", thresh, sep=""))

##### 600 #

load("TF.db.500.noNeph.clusters.RData")
cd4_bed <- read.table("/usr/projects/GROseq/NHP/tss_caller/HCM-U-PI.dREG-tss-clusters.tsv")[,c(1:3)]
cd4 <-  scanDb_rtfbs(tfs, cd4_bed, ncores=10, threshold=thresh, return_type="writedb", file_prefix=paste("cd4.dreg.tfs500.noNeph.th", thresh, sep=""))


load("TF.db.600.clusters.RData")
cd4_bed <- read.table("/usr/projects/GROseq/NHP/tss_caller/HCM-U-PI.dREG-tss-clusters.tsv")[,c(1:3)]
cd4 <-  scanDb_rtfbs(tfs, cd4_bed, ncores=10, threshold=thresh, return_type="writedb", file_prefix=paste("cd4.dreg.tfs600.th", thresh, sep=""))

