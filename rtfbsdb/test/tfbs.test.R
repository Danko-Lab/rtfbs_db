#
# Test a full procedure of tfbs class
#

file.bigwig.plus  <- "/home/zw355/src/rtfbs_db/testdata/GSM1480327_K562_PROseq_plus.bw";
file.bigwig.minus <- "/home/zw355/src/rtfbs_db/testdata/GSM1480327_K562_PROseq_minus.bw";
file.hg19         <- "/local/storage/data/hg19/hg19.2bit";
file.gencode.v19  <- "/local/storage/data/gencode/gencode.v19.annotation.gtf";

library(rtfbsdb)

db <- CisBP.extdata("Homo_sapiens");
db.sum <- CisBP.group( db, group.by="family_name", tf.information.type=1);
tfs <- tfbs.createFromCisBP( db );

tfs <- tfbs.getExpression( tfs, file.hg19, file.gencode.v19, file.bigwig.plus, file.bigwig.minus );

cluster.mat <- tfbs.clusterMotifs(tfs, pdf.heatmap="test-correlationMatrix.pdf")

tfbs.drawLogo( tfs ,1 )

usemotifs1 <- tfbs.selectByRandom(tfs, cluster.mat);

usemotifs2 <- tfbs.selectByGeneExp( tfs, cluster.mat );

save.image(file="test.tfbs.rdata");

q("no")
