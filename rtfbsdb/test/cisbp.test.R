#
# Test a full procedure of CisBP.db class
#

library(rtfbsdb)

#downloa human databse
db1 <- CisBP.download("Homo_sapiens");

#downloa mouse databse
db2 <- CisBP.download("Mus_musculus");

#reading data from zip file
db3 <- CisBP.zipload(db2@zip.file, species="Mus_musculus");

#reading data from inner file
db4 <- CisBP.extdata(species="Homo_sapiens");


#group method

gr1 <-CisBP.group(db1, group.by="family_name", tf.information.type=1 ); 

gr2 <-CisBP.group(db2, group.by="tf_status", tf.information.type=1 ); 

gr3 <-CisBP.group(db3, group.by="tf_status", tf.information.type=2); 

gr4 <-CisBP.group(db4, group.by="tf_status", tf.information.type=3); 

gr5 <-CisBP.group(db3, group.by="motif_type", tf.information.type=3); 

gr6 <-CisBP.group(db3, group.by="msource_id", tf.information.type=3); 

save(gr1, gr2, gr3, gr4, gr5, gr6, file="cisbp.test.rdata");

#find method
tfs0 <-tfbs.createFromCisBP(db4, family_name="Homeodomain", tf_status="D",  motif_type="ChIP-seq", msource_id= "MS01_1.01", tf.information.type=1 ); 

tfs1 <-tfbs.createFromCisBP(db4, family_name="Homeodomain", tf_status="D" ); 

tfs2 <-tfbs.createFromCisBP(db4, motif_type="ChIP-seq" ); 

tfs3 <-tfbs.createFromCisBP(db4, tf.information.type=2 ); 

save.image(file="cisbp.test.rdata");