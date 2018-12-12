#
# Test a tfbs.scanTSsite procedure of tfbs class
#

library(rtfbsdb)

detach("package:rtfbsdb", unload=TRUE);

library(rtfbsdb)


file.dREG.all.bed    <- "/workdir/zw355/proj/prj10-dreg/new-rf-201803/G1/G1.dREG.peak.score.bed.gz" 
file.twoBit_path     <-  "/fs/cbsudanko/storage/data/hg19/hg19.2bit"

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db);

dREG_bed <- read.table(file.dREG.all.bed, header=F);
dREG_bed <- dREG_bed[ dREG_bed[,1]=="chr1", ]

dScan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_bed, return.type="writedb", file.prefix="todo-rm", ncores=12);
    
dScan

dScan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_bed, ncores=12);

save(dScan, file="G1.dREG.peak.rdata");

for (i in 1:NROW(dScan$result)) write.table(dScan$result[[i]], file=paste("dREG-", dScan$summary[i,2], "-", dScan$summary[i,1], ".bed", sep=""), quote=F, row.names=F, col.names=F, sep="\t");
