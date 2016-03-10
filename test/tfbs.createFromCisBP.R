#
# Test a full procedure of CisBP.db class
#
if(require(rtfbsdb))
	detach("package:rtfbsdb", unload=TRUE);

library(rtfbsdb)

file.dREG.H.change.bed <- "/home/zw355/src/rtfbs_db/testdata/dREG.H.change.bed"
file.dREG.all.bed    <- "/home/zw355/src/rtfbs_db/testdata/dREG.all.bed" 
file.twoBit_path     <- "/local/storage/data/hg19/hg19.2bit";
file.bigwig.plus <- "/home/zw355/src/rtfbs_db/testdata/GSM1480327_K562_PROseq_plus.bw"
file.bigwig.minus <- "/home/zw355/src/rtfbs_db/testdata/GSM1480327_K562_PROseq_minus.bw"
file.gencode.gtf <- "/local/storage/data/gencode/gencode.v19.annotation.gtf"

db <- CisBP.extdata("Homo_sapiens");

tfs1 <- tfbs.createFromCisBP(db, family_name="AP-2" );
tfs1 <- tfbs.selectExpressedMotifs( tfs1, file.twoBit_path, 
			file.gencode.gtf, 
			file.bigwig.plus  = file.bigwig.plus, 
			file.bigwig.minus = file.bigwig.minus, 
			seq.datatype      = "PRO-seq", 
			ncores=1);

tfs2 <- tfbs.createFromCisBP(db);
tfs2 <- tfbs.selectExpressedMotifs( tfs2, file.twoBit_path, 
			file.gencode.gtf, 
			file.bigwig.plus  = file.bigwig.plus, 
			file.bigwig.minus = file.bigwig.minus,  
			seq.datatype      = "PRO-seq", 
			pvalue.threshold  = 0.005, 
      		include.DBID.missing=F, 
      		ncores            = 21);


file.bam <- "/local/storage/projects/NHP/AllData/bams/H3_U.fastq.gz.sort.bam";

tfs3 <- tfbs.createFromCisBP(db, family_name="AP-2");
tfs3 <- tfbs.selectExpressedMotifs(tfs3, file.twoBit_path, 
			file.gencode.gtf, 
			file.bam    = file.bam,
			seq.datatype="RNA-seq", 
			ncores      = 1);

tfs4 <- tfbs.createFromCisBP(db);
tfs4 <- tfbs.selectExpressedMotifs(tfs4, file.twoBit_path, 
			file.gencode.gtf, 
			file.bam         = file.bam,
			seq.datatype     = "RNA-seq", 
			pvalue.threshold = 0.005, 
      		include.DBID.missing=F, 
			ncores           = 21);

file.bam <- "/local/storage/projects/NHP/AllData/bams/H2_U.fastq.gz.sort.bam";
tfs5 <- tfbs.createFromCisBP(db);
tfs5 <- tfbs.selectExpressedMotifs(tfs5, file.twoBit_path, 
			file.gencode.gtf, 
			file.bam         = file.bam,
			seq.datatype     = "RNA-seq", 
			pvalue.threshold = 0.005, 
      		include.DBID.missing=F, 
			ncores           = 21);


save.image(file="tfbs.createFromCisBP.rdata");