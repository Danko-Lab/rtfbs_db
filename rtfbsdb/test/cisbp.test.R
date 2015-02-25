#
# Test a full procedure of tfbs.db class
#

library(rtfbsdb)

#downloa human databse
db1 <- CisBP.download("Homo_sapiens");

#downloa mouse databse
db2 <- CisBP.download("Mus_musculus");

#reading data from zip file
db3 <- CisBP.zipload(db2@zip.file, species="Mus_musculus");

#reading data from inner file
db4 <- CisBP.extdata(species="Mus_musculus");


#group method

gr1 <-tfbs.group(db3, group_by="family_name", motif_info_type=1 ); 

gr2 <-tfbs.group(db3, group_by="tf_status", motif_info_type=1 ); 

gr3 <-tfbs.group(db3, group_by="tf_status", motif_info_type=2); 

gr4 <-tfbs.group(db3, group_by="tf_status", motif_info_type=3); 

gr5 <-tfbs.group(db3, group_by="motif_type", motif_info_type=3); 

gr6 <-tfbs.group(db3, group_by="msource_id", motif_info_type=3); 

save(gr1, gr2, gr3, gr4, gr5, gr6, file="cisbp.test.rdata");

#find method
tfs0 <-tfbs.find(db2, family_name="Homeodomain", tf_status="D",  motif_type="ChIP-seq", msource_id= "MS01_1.01", motif_info_type=1 ); 

tfs1 <-tfbs.find(db2, family_name="Homeodomain", tf_status="D",  motif_info_type=1 ); 

tfs2 <-tfbs.find(db2, motif_type="ChIP-seq", motif_info_type=1 ); 

tfs3 <-tfbs.find(db2, motif_info_type=2); 

save.image(file="cisbp.test.rdata");