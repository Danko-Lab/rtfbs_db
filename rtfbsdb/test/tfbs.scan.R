#
# Test a tfbs.scanTSsite procedure of tfbs class
#

library(rtfbsdb)


file.dREG.H.change.bed <- "/home/zw355/src/rtfbs_db/testdata/dREG.H.change.bed"
file.dREG.all.bed    <- "/home/zw355/src/rtfbs_db/testdata/dREG.all.bed" 
file.twoBit_path     <- "/home/zw355/src/rtfbs_db/testdata/hg19.2bit"


db <- CisBP.extdata("Homo_sapiens");
tfs <- CisBP.find(db, family_name="AP-2");

dREG_H_change_bed <- read.table(file.dREG.H.change.bed, header=F);

t1.bed <- tfbs.scanTFsite( tfs, 
   file.twoBit_path, 
   bed_dat=dREG_H_change_bed, 
   return_posteriors=FALSE,
   threshold = 7);
    
save(t1.bed, tfs, file="tfbs.scan.rdata");