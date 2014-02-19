## 
## Generates a tfbs object that stores
## PWMs for all motifs in R.
## 

thresh=0

require(atif)
setwd("/usr/data/GROseq.parser/pwm_data")
load("TF.db.263.clusters.RData")

## k562
k562_bed <- read.table("/usr/data/GROseq.parser/hg19/k562/dnase/wgEncodeOpenChromDnaseK562PkV2.narrowPeak.gz")
k562 <- scanDb_rtfbs(tfs, k562_bed, ncores=10, threshold=thresh, return_type="writedb", file_prefix="k562.tfs263.")

## CD4
cd4_bed <- read.table("/usr/data/GROseq.parser/hg19/cd4/dnase1fp/dnase1.peaks_peaks.narrowPeak")
cd4 <-  scanDb_rtfbs(tfs, cd4_bed, ncores=10, threshold=thresh, return_type="writedb", file_prefix="cd4.tfs263.")

## GM12878
gm12878_bed <- read.table("/usr/data/GROseq.parser/hg19/gm12878/dnase/wgEncodeOpenChromDnaseGm12878Pk.narrowPeak.gz")
gm12878 <-  scanDb_rtfbs(tfs, gm12878_bed, ncores=10, threshold=thresh, return_type="writedb", file_prefix="gm12878.tfs263.")

