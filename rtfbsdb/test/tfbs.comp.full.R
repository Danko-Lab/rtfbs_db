library(rtfbsdb)

file.dREG.H.change.bed <- "/home/zw355/src/rtfbs_db/testdata/dREG.H.change.bed"
file.dREG.all.bed      <- "/home/zw355/src/rtfbs_db/testdata/dREG.all.bed" 
file.twoBit_path       <- "/local/storage/data/hg19/hg19.2bit";

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db);

dREG_H_change_bed <- read.table(file.dREG.H.change.bed, header=F);
dREG_all_bed <- read.table(file.dREG.all.bed, header=F);

t <- tfbs.enrichmentTest( tfs, 
	file.twoBit_path, 
	dREG_H_change_bed, 
	dREG_all_bed, 
	gc.correction=TRUE,
	ncores = 21);

tfbs.reportEnrichment(tfs, t, file.pdf="tfbs.comp.full.pdf", sig.only=F, report.title="TEST FULL");

save.image(file="tfbs.comp.full.rdata");
