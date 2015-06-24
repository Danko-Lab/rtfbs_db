#
# Test a full procedure of CisBP.db class
#

detach("package:rtfbsdb", unload=TRUE);

library(rtfbsdb)

file.dREG.H.change.bed <- "/home/zw355/src/rtfbs_db/testdata/dREG.H.change.bed"
file.dREG.all.bed    <- "/home/zw355/src/rtfbs_db/testdata/dREG.all.bed" 
file.twoBit_path     <- "/local/storage/data/hg19/hg19.2bit";
file.bigwig.plus <- "/home/zw355/src/rtfbs_db/testdata/GSM1480327_K562_PROseq_plus.bw"
file.bigwig.minus <- "/home/zw355/src/rtfbs_db/testdata/GSM1480327_K562_PROseq_minus.bw"
file.gencode.gtf <- "/local/storage/data/gencode/gencode.v21.annotation.gtf"

db <- CisBP.extdata("Homo_sapiens");

tfs1 <- tfbs.createFromCisBP(db, family_name="AP-2", file.bigwig.plus=file.bigwig.plus, file.bigwig.minus=file.bigwig.minus,  file.gencode.gtf=file.gencode.gtf, seq.datatype="PRO-seq", ncores=1);

tfs2 <- tfbs.createFromCisBP(db, file.bigwig.plus=file.bigwig.plus, file.bigwig.minus=file.bigwig.minus,  file.gencode.gtf=file.gencode.gtf, seq.datatype="PRO-seq", ncores=21);


file.bam.plus <- "/local/storage/projects/NHP/AllData/bams/H3_U.fastq.gz.sort.bam";

tfs3 <- tfbs.createFromCisBP(db, family_name="AP-2", file.bigwig.plus=file.bam.plus, file.bigwig.minus=NULL,  file.gencode.gtf=file.gencode.gtf, seq.datatype="RNA-seq", ncores=1);

tfs4 <- tfbs.createFromCisBP(db, file.bigwig.plus=file.bam.plus, file.bigwig.minus=NULL,  file.gencode.gtf=file.gencode.gtf, seq.datatype="RNA-seq", ncores=21);


save.image(file="tfbs.createFromCisBP.rdata");