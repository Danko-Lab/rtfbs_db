detach("package:rtfbsdb", unload=TRUE);

library(rtfbsdb)

db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP( db );

motif_id   <- c( "M5604_1.01", "M6001_1.01" , "M5440_1.01", "M5441_1.01", "M5162_1.01", "M5352_1.01");
tf_id      <- c( "T093250_1.01", "T093251_1.01","T093252_1.01","T093253_1.01");
family_name<- c( "p53", "Homeodomain", "Paired box", "Pipsqueak");


#tfbs.drawLogo(tfbs, file.pdf, index=NA, tf_id=NA, motif_id=NA, tf_name=NA, family_name=NA, tf_status=NA, groupby=NA) 

tfbs.drawLogo(tfs, file.pdf="tfbs.drawLogo1.pdf", index=c(1:100), groupby="Family_Name");

tfbs.drawLogo(tfs, file.pdf="tfbs.drawLogo2.pdf", motif_id=motif_id, tf_id=tf_id, tf_name="AP-2", family_name=family_name, groupby="TF_Status");

tfbs.drawLogo(tfs, file.pdf="tfbs.drawLogo3.pdf", tf_status="D", groupby="TF_Status");

tfbs.drawLogo(tfs, file.pdf="tfbs.drawLogo4.pdf", index=c(1:100), groupby="TF_Status");

tfbs.drawLogo(tfs, file.pdf="tfbs.drawLogo5.pdf", groupby="TF_Status" );
 