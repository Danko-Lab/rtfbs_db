#
# Test a full procedure of tfbs class
#

file.bitwig.plus <- "/work/03350/tg826494/hg19/GSM1480327_K562_PROseq_plus.bw";
file.bigwig.minus <- "/work/03350/tg826494/hg19/GSM1480327_K562_PROseq_minus.bw";
file.hg19 <- "/work/03350/tg826494/hg19/hg19.2bit";

library(rtfbsdb)

db <- CisBP.extdata();
db.sum <- tfbs.group( db, group_by="family_name", motif_info_type=1);
tfs <- tfbs.find( db, motif_info_type=1);

save(tfs, file="tfbs.rdata");

tfs <- tfbs.getDistanceMatrix(tfs, ncores=5)

save(tfs, file="tfbs.rdata");

pdf(paste("correlationMatrix.500.pdf", sep=""))
tfs <- tfbs.clusterMotifs(tfs, 500, draw_heatmap= TRUE)
dev.off()

pdf(paste("logos.all.pdf", sep=""))
tfbs.drawLogo( tfs ,1 )
dev.off()

pdf(paste("logos.500.pdf", sep=""))
tfs <- tfbs.setUseMotifs.random(tfs);
tfbs.drawLogosForClusters( tfs )
dev.off()

tfs <- tfbs.getExpression( tfs, NULL, file.bitwig.plus, file.bigwig.minus, file.hg19 );

pdf(paste("logos.500.pdf", sep=""))
tfs <- tfbs.selectByGeneExp( tfs );
tfbs.drawLogosForClusters( tfs )
dev.off()

save(tfs, file="tfbs.rdata");

q("no")
