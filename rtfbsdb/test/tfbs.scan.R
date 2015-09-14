#
# Test a tfbs.scanTSsite procedure of tfbs class
#


detach("package:rtfbsdb", unload=TRUE);

library(rtfbsdb)


file.dREG.H.change.bed <- "/home/zw355/src/rtfbs_db/testdata/dREG.H.change.bed"
file.dREG.all.bed    <- "/home/zw355/src/rtfbs_db/testdata/dREG.all.bed" 
file.twoBit_path     <- "/local/storage/data/hg19/hg19.2bit";


db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db, family_name="AP-2");

dREG_H_change_bed <- read.table(file.dREG.H.change.bed, header=F);


t1.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, ncores=3);
    
t2.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, return.type="writedb", ncores=3);

t3.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, return.type="posteriors", ncores=1, threshold=8);

t4.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, return.type="maxposterior", ncores=1, threshold=8);

t5.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, ncores=3, file.prefix = "test.tfscan", threshold=0.05, threshold.type="fdr" );

t6.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, ncores=3, file.prefix = "test.tfscan", threshold=0.05, threshold.type="fdr", gc.groups=5, background.order = 3 );

t7.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, ncores=3, threshold=8, gc.groups=4, background.order = 3 );

t8.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, ncores=3, return.type="writedb", file.prefix = "test.tfscan", threshold=0.05, threshold.type="fdr", gc.groups=5, background.order = 3 );

t9.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, ncores=3, return.type="matches", file.prefix = "test.tfscan", threshold=0.05, threshold.type="fdr", gc.groups=5, background.order = 3 );


t1.scan;

t2.scan;

t3.scan;

t4.scan;

t5.scan;

t6.scan;

t7.scan;

t8.scan;

t9.scan;

tfbs.reportFinding(tfs, t9.scan,  file.pdf="test-tfscan.pdf", report.title="AP-2 Results");

