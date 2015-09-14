Change Log
========
09/14/2015
A: Data:
1) Adding the CisBP data file for dm3(Drosophila_melanogaster)
2) Adding partial hg19.2bit, partial gencode GTF for chrmosome 19
3) Adding BED files from ChipSeq data for chrmosome 19
4) Adding PRO-seq data for chrmosome 19
5) some fake motif PWM files.

B: Source: 
1) Adding function: CisBP.getTFinformation()
2) Adding slots "zip.date" for CisBP.db class
3) Adding slots "cluster" and "tf_missing", renaming slots "tf_info", removing TF_ID for tfbs class
4) Adding function: tfbs.importMotifs
5) Adding function:tfbs.selectExpressedMotifs(), gene expression and motif filtering, separated from the constructor CisBP.createFromCisBP.
6) Renaming function: tfbs.enrichmentTest()
7) Adding new cluster (apcluster) for tfbs.clusterMotifs
8) Renaming expected field for the results returned by tfbs.enrichmentTest.
9) Adding parameter gc.robust.rep, use.cluster, threshold.rype for tfbs.enrichmentTest.
10) Adding parameter threshold.rype for tfbs.scanTFsite.
11) Adding parameter enrichment.type and pv.threshold for tfbs.reportEnrichment
12) The parameter orders are changed in some functions.

