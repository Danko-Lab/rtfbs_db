library(rtfbsdb)
fs.dir <- system.file("extdata","", package="rtfbsdb")
tfs <- tfbs.dirs( fs.dir, args.read.motif = list(pseudocount=-7, header=TRUE, sep="\t" , row.names=1) );
 
tfbs.drawLogo(tfsN, file.pdf="/home/zw355/temp/a/a.pdf");
 
file.twoBit <- "/home/zw355/src/rtfbs_db/testdata/hg19.2bit"
gen.bed <- data.frame(chr="chr19", start=round(runif(10,1000000, 2000000)),stop=0, name="", score=0, strand=".");
gen.bed$stop <- gen.bed$start + 3000;
t1 <- tfbs.scanTFsite( tfs, file.twoBit, gen.bed=NULL, file.prefix="test.db", ncores = 1);


library(rtfbsdb)
db <- CisBP.extdata("Mus_musculus")
tfRORb <- tfbs.createFromCisBP(db,tf_name = "Rorb",motif_id = NULL,msource_id = NULL,tf.information.type = 1)


library(rtfbsdb)
fs.dir <- "/home/zw355/temp"
tfsN <- tfbs.dirs(fs.dir, species = "Mus_musculus", args.read.motif = list(header=TRUE, sep="\t", force_even=TRUE, row.names=1), pattern = glob2rx("*.pwm"))

tfsN;

tfbs.drawLogo(tfsN, file.pdf="/home/zw355/temp/a/a.pdf");

file.twoBit <- "/home/zw355/src/rtfbs_db/testdata/hg19.2bit"
gen.bed <- data.frame(chr="chr19", start=round(runif(10,1000000, 2000000)),stop=0, name="", score=0, strand=".");
gen.bed$stop <- gen.bed$start + 3000;
t1 <- tfbs.scanTFsite( tfsN, file.twoBit, gen.bed=NULL, file.prefix="test.db", ncores = 1);



tfsN <- tfbs.dirs(fs.dir, species = "Mus_musculus", args.read.motif = list(header=TRUE, sep="\t", force_even=TRUE), pattern = glob2rx("*.pwm"))
tfsN;
tfsN <- tfbs.dirs(fs.dir, species = "Mus_musculus", args.read.motif = list(header=TRUE, sep="\t", row.names=1), pattern = glob2rx("*.pwm"))
tfsN;
tfsN <- tfbs.dirs(fs.dir, species = "Mus_musculus", args.read.motif = list(header=TRUE, sep="\t"), pattern = glob2rx("*.pwm"))
tfsN;


