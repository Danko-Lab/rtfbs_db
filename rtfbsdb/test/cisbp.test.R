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

gr1 <-CisBP.group(db1, group_by="family_name", motif_info_type=1 ); 

gr2 <-CisBP.group(db2, group_by="tf_status", motif_info_type=1 ); 

gr3 <-CisBP.group(db3, group_by="tf_status", motif_info_type=2); 

gr4 <-CisBP.group(db4, group_by="tf_status", motif_info_type=3); 

gr5 <-CisBP.group(db3, group_by="motif_type", motif_info_type=3); 

gr6 <-CisBP.group(db3, group_by="msource_id", motif_info_type=3); 

save(gr1, gr2, gr3, gr4, gr5, gr6, file="cisbp.test.rdata");

#find method
tfs0 <-CisBP.find(db4, family_name="Homeodomain", tf_status="D",  motif_type="ChIP-seq", msource_id= "MS01_1.01", motif_info_type=1 ); 

tfs1 <-CisBP.find(db4, family_name="Homeodomain", tf_status="D",  motif_info_type=1 ); 

tfs2 <-CisBP.find(db4, motif_type="ChIP-seq", motif_info_type=1 ); 

tfs3 <-CisBP.find(db4, motif_info_type=2); 

save.image(file="cisbp.test.rdata");