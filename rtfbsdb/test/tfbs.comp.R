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

tfs <- tfbs.clusterMotifs(tfs, pdf.heatmap="test.tfbs.comp.pdf")

t1.comp <- tfbs.enrichmentTest( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correction=F, use.cluster=T, ncores = 7);

t2.comp <- tfbs.enrichmentTest( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correction=T, gc.robust.rep=5, gc.correction.pdf="gc.correction.pdf", ncores = 7);

t3.comp <- tfbs.enrichmentTest( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, file.prefix="test.db", ncores = 7, threshold=0.05, threshold.type="fdr");

t4.comp <- tfbs.enrichmentTest( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correction=TRUE, file.prefix="test.db", ncores = 7, threshold=8);

t5.comp <- tfbs.enrichmentTest( tfs, file.twoBit_path, dREG_H_change_bed, file.prefix="test.db", ncores = 7, threshold=0.05, threshold.type="fdr", background.order = 3 );

t6.comp <- tfbs.enrichmentTest( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correction=TRUE, gc.robust.rep=5, file.prefix="test.db", );

t7.comp <- tfbs.enrichmentTest( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correction=FALSE, ncores = 7, threshold=0.05, threshold.type="fdr", background.order = 3 , pv.adj = "fdr");

t8.comp <- tfbs.enrichmentTest( tfs, file.twoBit_path, dREG_H_change_bed, dREG_all_bed, gc.correction=TRUE, gc.robust.rep=5, ncores = 7, threshold=8, background.order = 3 , pv.adj = "fdr");

t9.comp <- tfbs.enrichmentTest( tfs, file.twoBit_path, dREG_H_change_bed, gc.correction=TRUE, use.cluster=TRUE, ncores = 7, threshold=0.05, threshold.type="fdr", background.order = 3 , pv.adj = "fdr");

save.image("tfbs.comp.rdata");

t1.comp;

t2.comp;

t3.comp;

t4.comp;

t5.comp;

t6.comp;

t7.comp;

t8.comp;

t9.comp;

tfbs.reportEnrichment(tfs, t1.comp, file.pdf="test.tfbs.comp.report1.pdf", sig.only=F, report.title="Test Report");

tfbs.reportEnrichment(tfs, t9.comp, file.pdf="test.tfbs.comp.report2.pdf", sig.only=T, enrichment.type="depleted", report.title="Significant Report", pv.threshold=0.1, pv.adj="fdr");

