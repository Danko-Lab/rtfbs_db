#
# Test a full procedure of tfbs class
#

file.bigwig.plus  <- "/home/zw355/src/rtfbs_db/testdata/GSM1480327_K562_PROseq_plus.bw";
file.bigwig.minus <- "/home/zw355/src/rtfbs_db/testdata/GSM1480327_K562_PROseq_minus.bw";
file.hg19 <- "/home/zw355/src/rtfbs_db/testdata/hg19.2bit";

library(rtfbsdb)

db <- CisBP.extdata("Homo_sapiens");
db.sum <- CisBP.group( db, group_by="family_name", motif_info_type=1);
tfs <- CisBP.find( db );

save(tfs, file="tfbs.rdata");

tfs <- tfbs.getExpression( tfs, file.bigwig.plus, file.bigwig.minus, file.hg19 );

tfs <- tfbs.getDistanceMatrix(tfs, ncores=5);

save(tfs, file="tfbs.rdata");

cluster.mat <- tfbs.clusterMotifs(tfs, pdf.heatmap="correlationMatrix.pdf")

tfbs.drawLogo( tfs ,1 )

usemotifs1 <- tfbs.selectByRandom(tfs, cluster.mat);

usemotifs2 <- tfbs.selectByGeneExp( tfs, cluster.mat );

save(tfs, cluster.mat, usemotifs1, usemotifs2, file="tfbs.rdata");

q("no")
