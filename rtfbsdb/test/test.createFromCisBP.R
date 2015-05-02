#
# Test a full procedure of CisBP.db class
#

library(rtfbsdb)

file.dREG.H.change.bed <- "/home/zw355/src/rtfbs_db/testdata/dREG.H.change.bed"
file.dREG.all.bed    <- "/home/zw355/src/rtfbs_db/testdata/dREG.all.bed" 
file.twoBit_path     <- "/local/storage/data/hg19/hg19.2bit";
file.bigwig.plus <- "/home/zw355/src/rtfbs_db/testdata/GSM1480327_K562_PROseq_plus.bw"
file.bigwig.minus <- "/home/zw355/src/rtfbs_db/testdata/GSM1480327_K562_PROseq_minus.bw"
file.gencode.gtf <- "/local/storage/data/gencode/gencode.v21.annotation.gtf"

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, family_name="AP-2", file.bigwig.plus=file.bigwig.plus, file.bigwig.minus=file.bigwig.minus,  file.gencode.gtf=file.gencode.gtf, seq.datatype="RNA-seq" );


save.image(file="cisbp.test.rdata");