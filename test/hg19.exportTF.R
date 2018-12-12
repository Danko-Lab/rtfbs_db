library(rtfbsdb)


file.dREG.all.bed    <- "/home/zw355/src/rtfbs_db/test/G1.dREG.peak.score.bed.gz"
file.twoBit_path     <- "/local/storage/data/hg19/hg19.2bit";


db <- CisBP.extdata("Homo_sapiens");
tfs <- tfbs.createFromCisBP(db);
df <- as.data.frame(table(tfs@tf_info$Motif_Type))
df <- df[df[,2]!=0,]

tf_info <- tfs@tf_info[, c("Motif_ID", "TF_Name")];
tf_info <- do.call("rbind", lapply(unique(tf_info$Motif_ID), function(x){ 
          idx <- which( as.character(tf_info$Motif_ID) == as.character(x) ); 
          if(NROW(idx)>4)
            return(data.frame(Motif_ID=x, TF_Name=paste(c(as.character(tf_info[idx, "TF_Name"])[1:4], "..."), collapse=",")))
          else if( NROW(idx)>1 )
             return(data.frame(Motif_ID=x, TF_Name=paste(as.character(tf_info[idx, "TF_Name"]), collapse=",")))
          else 
            return(data.frame(Motif_ID=x, TF_Name=as.character(tf_info[idx, "TF_Name"])))
        }));

write.table(tf_info, file="hg19.tf.info", quote=F, row.names=F, col.names=F, sep="\t");

for (Motif_Type in df$Var1)
{
  if (!file.exists(paste0("hg19_Motif_Type_", Motif_Type, ".rdata")))
  {
  tfs1 <- tfbs.createFromCisBP(db, motif_type=Motif_Type);
  dScan <- tfbs.scanTFsite( tfs1, file.twoBit_path, ncores=5);
  save(dScan, file=paste0("hg19_Motif_Type_", Motif_Type, ".rdata") );
  rm(dScan);
  gc();
  }
}

for (Motif_Type in df$Var1)
{
  if (file.exists(paste0("hg19_Motif_Type_", Motif_Type, ".rdata")))
  {
    load( paste0("hg19_Motif_Type_", Motif_Type, ".rdata") ) ;

    for (i in 1:NROW(dScan$summary))
    {
	    dScan.bed <- dScan$result[[i]];
	    out_file <- paste0("temp/", as.character(dScan$summary$Motif_ID[i]), ".bed");
	    if(!file.exists(paste0(out_file, ".gz") ) )
          make_bed( dScan.bed, out_file);
    }
  }
}