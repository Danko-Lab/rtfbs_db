library(rtfbsdb)


## Reference http://rsat01.biologie.ens.fr/rsa-tools/help.convert-matrix.html

tbs <- tfbs();

## Data is copied from http://rsat01.biologie.ens.fr/rsa-tools/help.convert-matrix.html
##
data.transfac <- system.file("extdata", "pwm.example.transfac.txt", package="rtfbsdb");
tfs.transfac <- tfbs.importMotifs(tbs, "transfac", data.transfac, skip.lines=0);
show(tfs.transfac);

## Data from http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt
##
data.jaspar <- system.file("extdata", "pwm.example.jaspar.2015.txt", package="rtfbsdb");
tfs.jaspar <- tfbs.importMotifs(tbs, "jaspar", data.jaspar, skip.lines=0);
show(tfs.jaspar);

## Data from http://jaspar.genereg.net/html/DOWNLOAD/ARCHIVE/JASPAR2010/JASPAR_CORE/pfms/pfms_all.txt
##
data.mscan <- system.file("extdata", "pwm.example.mscan.jaspar2010.txt", package="rtfbsdb");
tfs.mscan <- tfbs.importMotifs(tbs, "mscan", data.mscan, skip.lines=0);
show(tfs.mscan);

## Data from http://meme-suite.org/doc/download.html?man_type=web
##
data.meme <- system.file("extdata", "pwm.example.meme.Homo_sapiens.txt", package="rtfbsdb");
tfs.meme <- tfbs.importMotifs(tbs, "meme", data.meme, skip.lines=5);
show(tfs.meme);

format.style <- c(">$Motif_ID $TF_Name",
		"$A+",
		"$C+",
		"$G+",
		"$T+");
data.jaspar2010 <- system.file("extdata", "pwm.example.mscan.jaspar2010.txt", package="rtfbsdb");
tfs.jaspar2010 <- tfbs.importMotifs(tbs, format.style, data.jaspar2010, skip.lines=0)
show(tfs.jaspar2010);

## Data from http://www.nature.com/nature/journal/v527/n7578/full/nature15518.html#supplementary-information "S2 PWM moldes"
##                  
format.style <- c("Base $Motif_ID $TF_Name $Experiment $Ligand_sequbatch $Seed $Multinomial $Cycle $Is_matrix $Comment $Genomic_pvalue $Enrichment_pvalue $Category $SKIP",
		"A $A+",
		"C $C+",
		"G $G+",
		"T $T+" );
data.file <- system.file("extdata", "pwm.example.nature15518.s1.txt", package="rtfbsdb");
tfs.nature15518 <- tfbs.importMotifs(tbs, format.style, data.file, skip.lines=19)
show(tfs.nature15518);

options(warn=2);

tbs <- tfbs();
tbs <- tfbs.importMotifs(tbs, "transfac", data.transfac );
tbs <- tfbs.importMotifs(tbs, "jaspar", data.jaspar, skip.lines=0);

db <- CisBP.extdata("human");
tfs <- tfbs.createFromCisBP(db, tf_name="ATF1")
tfs <- tfbs.importMotifs(tfs, "transfac", data.transfac );
tfs <- tfbs.importMotifs(tfs, "jaspar", data.jaspar, skip.lines=0);

file.dREG.H.change.bed <- "/home/zw355/src/rtfbs_db/testdata/dREG.H.change.bed"
file.dREG.all.bed    <- "/home/zw355/src/rtfbs_db/testdata/dREG.all.bed" 
file.twoBit_path     <- "/local/storage/data/hg19/hg19.2bit";

dREG_H_change_bed <- read.table(file.dREG.H.change.bed, header=F);
dREG_all_bed <- read.table(file.dREG.all.bed, header=F);

t1.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, ncores=7);
show(t1.scan);

t2.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, return.type="writedb", ncores=7);
show(t2.scan);

t3.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, return.type="posteriors", ncores=7, threshold=8);
show(t3.scan);

t4.scan <- tfbs.scanTFsite( tfs, file.twoBit_path, dREG_H_change_bed, return.type="maxposterior", ncores=7, threshold=8);
show(t4.scan);

t <- tfbs.enrichmentTest( tfs, 
	file.twoBit_path, 
	dREG_H_change_bed, 
	dREG_all_bed, 
	gc.correction=TRUE,
	ncores = 21);

tfbs.reportEnrichment(tfs, t, file.pdf="tfbs.importMotifs.pdf", sig.only=F, report.title="IMPORT TEST");

save.image(file="tfbs.importMotifs.rdata");

