load("motif_intersected_sorted_uniq_sortedByAbsLogRatio.rdata")

library(rtfbsdb);

#db <- CisBP.extdata("Homo_sapiens");
#tfs <- tfbs.createFromCisBP(db);

load("tfs.rdata");

tfbs.plotEnrichment(tfs, t.comp, file.pdf="enrich-qqplot2.pdf", plot.type="nonpolar", enrichment.type="both", options=list( width=4, zoom.motif.label=1.5, zoom.motif.logo=2, zoom.tick=1, zoom.label=1, zoom.legend.label=1 ))

tfbs.plotEnrichment(tfs, t.comp, file.pdf="enrich-qqplot1.pdf", plot.type="nonpolar", enrichment.type="enriched")

tfbs.plotEnrichment(tfs, t.comp, file.pdf="enrich-qqplot3.pdf", plot.type="polar", enrichment.type="depleted")

tfbs.plotEnrichment(tfs, t.comp, file.pdf="enrich-qqplot4.pdf", plot.type="polar", enrichment.type="both")

