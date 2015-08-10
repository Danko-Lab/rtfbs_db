#
# Test a tfbs.scanTSsite procedure of tfbs class
#

if(require(rtfbsdb))
detach("package:rtfbsdb", unload=TRUE);

library(rtfbsdb);

file.dREG.H.change.bed <- "/home/zw355/src/rtfbs_db/testdata/dREG.H.change.bed"
file.dREG.all.bed    <- "/home/zw355/src/rtfbs_db/testdata/dREG.all.bed" 
file.twoBit_path     <- "/local/storage/data/hg19/hg19.2bit";

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, family_name="AP-2");

dREG_H_change_bed <- read.table(file.dREG.H.change.bed, header=F);
dREG_all_bed <- read.table(file.dREG.all.bed, header=F);

tfs <- tfbs.getDistanceMatrix(tfs, ncores=5);

hcluster <- tfbs.clusterMotifs(tfs, pdf.heatmap="correlationMatrix.pdf")

t1.comp <- tfbs.compareTFsite( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, ncores = 7);

t2.comp <- tfbs.compareTFsite( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correct=TRUE, ncores = 7);

t3.comp <- tfbs.compareTFsite( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correct=TRUE, file.prefix="test.db", ncores = 7, fdr=0.05);

t4.comp <- tfbs.compareTFsite( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correct=TRUE, file.prefix="test.db", ncores = 7, threshold=8);

t5.comp <- tfbs.compareTFsite( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correct=TRUE, file.prefix="test.db", ncores = 7, fdr=0.05, background.order = 3 );

t6.comp <- tfbs.compareTFsite( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correct=TRUE, file.prefix="test.db", );

t7.comp <- tfbs.compareTFsite( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correct=FALSE, ncores = 7, fdr=0.05, background.order = 3 , pv.adj = "fdr");

t8.comp <- tfbs.compareTFsite( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correct=TRUE, ncores = 7, threshold=8, background.order = 3 , pv.adj = "fdr");

t9.comp <- tfbs.compareTFsite( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correct=TRUE, usecluster=hcluster, ncores = 7, fdr=0.05, background.order = 3 , pv.adj = "fdr");

save(t1.comp, t2.comp, t3.comp, t4.comp, t5.comp, t6.comp, t7.comp, t8.comp, t9.comp, tfs, file="tfbs.comp.rdata");

t1.comp;

t2.comp;

t3.comp;

t4.comp;

t5.comp;

t6.comp;

t7.comp;

t8.comp;

t9.comp;

tfbs.reportComparson(tfs, t1.comp, file.pdf="test-tfcomp1.pdf", sig.only=F, report.title="Test Report");

tfbs.reportComparson(tfs, t9.comp, file.pdf="test-tfcomp9.pdf", sig.only=T, report.title="Significant Report", pv.cutoff=0.1, pv.adj="fdr");
