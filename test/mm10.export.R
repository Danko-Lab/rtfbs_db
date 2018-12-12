library(rtfbsdb)

file.twoBit_path <- "/local/storage/data/mm10/mm10.2bit";

db <- CisBP.extdata("Mus_musculus");
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

write.table(tf_info, file="mm10.tf.info", quote=F, row.names=F, col.names=F, sep="\t");


for (Motif_Type in df$Var1)
{
  if (!file.exists(paste0("Motif_Type_", Motif_Type, ".rdata")))
  {
cat("Motif_Type=", Motif_Type, "\n");
  tfs1 <- tfbs.createFromCisBP(db, motif_type=Motif_Type);
  dScan <- tfbs.scanTFsite( tfs1, file.twoBit_path, ncores=20);
  save(dScan, file=paste0("Motif_Type_", Motif_Type, ".rdata") );
  rm(dScan);
  gc();
  }
}

make_bed<-function( df_bed, out_file)
{
    options(scipen=99, digits=4);
	file.tmp <- tempfile(fileext=".bed");
	write.table( df_bed, file=file.tmp, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
	system(paste( "sort-bed ", file.tmp, " | bgzip -f > ", out_file, ".gz", sep="") );
	unlink(file.tmp)
}


for (Motif_Type in df$Var1)
{
  if (file.exists(paste0("Motif_Type_", Motif_Type, ".rdata")))
  {
cat("load =",  Motif_Type, "\n");
    load( paste0("Motif_Type_", Motif_Type, ".rdata") ) ;

    for (i in 1:NROW(dScan$summary))
    {
	    dScan.bed <- dScan$result[[i]];
	    out_file <- paste0("temp/", as.character(dScan$summary$Motif_ID[i]), ".bed");
	    if(!file.exists(paste0(out_file, ".gz") ) )
          make_bed( dScan.bed, out_file);
    }
  }
}


