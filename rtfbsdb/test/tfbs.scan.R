#
# Test a tfbs.scanTSsite procedure of tfbs class
#

library(rtfbsdb)


file.dREG.H.change.bed <- "/work/03350/tg826494/hg19/dREG.H.change.bed"
file.dREG.all.bed    <- "/work/03350/tg826494/hg19/dREG.all.bed" 
file.twoBit_path     <- "/work/03350/tg826494/hg19/hg19.2bit"
file.hg19.diff.bed   <- "/work/03350/tg826494/hg19/hg19.diff.bed"


db <- CisBP.extdata();
tfs <- tfbs.find(db, family_name="AP-2");

dREG_H_change_bed <- read.table(file.dREG.H.change.bed, header=F);

t1.bed <- tfbs.scanTFsite( tfs, 
   file.twoBit_path, 
   bed_dat=dREG_H_change_bed, 
   return_posteriors=FALSE,
   threshold = 7);
   
save(t1.bed, tfs, file="tfbs.scan.rdata");