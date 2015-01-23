## 
## Generates a tfbs object that stores
## PWMs for all motifs in R.
## 

setwd("/usr/data/GROseq.parser/pwm_data")
require(rtfbs_db)

jolma <- paste("jolma/teal/", dir(path="jolma/teal/", pattern="pwm"), sep="")
jaspar<- paste("jaspar/", dir(path="jaspar/", pattern="pwm"), sep="")
neph  <- paste("neph/", dir(path="neph/", pattern="pwm"), sep="")

tf_names <- c(substr(jolma, 12, 16), substr(jaspar, 8, 12)) ##, paste("neph",substr(neph, 16, 18), sep=""))

tfs <- tfbs(c(jolma, jaspar), tf_names, header=TRUE) ##neph
tfs <- tfbs.getDistanceMatrix(tfs, ncores=5)


pdf(paste("correlationMatrix.500.pdf", sep=""))
tfs <- tfbs.clusterMotifs(tfs, 600, draw_heatmap= TRUE)
tfs <- tfbs.setUseMotifs.random(tfs)
dev.off()

pdf(paste("logos.500.pdf", sep=""))
tfbs.drawLogosForClusters(tfs)
dev.off()

save.image(paste("TF.db.500.noNeph.clusters.RData"))


q("no")


## 600
pdf(paste("correlationMatrix.600.pdf", sep=""))
tfs <- tfbs.clusterMotifs(tfs, 600, draw_heatmap= TRUE)
tfs <- tfbs.setUseMotifs.random(tfs)
dev.off()

pdf(paste("logos.600.pdf", sep=""))
tfbs.drawLogosForClusters(tfs)
dev.off()

save.image(paste("TF.db.600.clusters.RData"))

q("no")

## 263
pdf(paste("correlationMatrix.263.pdf", sep=""))
tfs <- tfbs.clusterMotifs(tfs, 263, draw_heatmap= TRUE)
tfs <- tfbs.setUseMotifs.random(tfs)
dev.off()

pdf(paste("logos.263.pdf", sep=""))
tfbs.drawLogosForClusters(tfs)
dev.off()

save.image(paste("TF.db.263.clusters.RData"))

## 341
pdf(paste("correlationMatrix.341.pdf", sep=""))
  tfs <- tfbs.clusterMotifs(tfs, 341, draw_heatmap= TRUE)
  tfs <- tfbs.setUseMotifs.random(tfs)
dev.off()

pdf(paste("logos.341.pdf", sep=""))
  tfbs.drawLogosForClusters(tfs)
dev.off()

save.image(paste("TF.db.341.clusters.RData"))
