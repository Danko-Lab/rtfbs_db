load("motif_intersected_sorted_uniq_sortedByAbsLogRatio.rdata")

library(rtfbs_db);

tfbs.plotEnrichment(tfbs, t.comp, file.pdf="enrich-qqplot.pdf", plot.type="nonpolar", enrichment.type="enriched")

tfbs.plotEnrichment(tfbs, t.comp, file.pdf="enrich-qqplot.pdf", plot.type="nonpolar", enrichment. type="both")

tfbs.plotEnrichment(tfbs, t.comp, file.pdf="enrich-qqplot.pdf", plot.type="polar", enrichment.type="depleted")

tfbs.plotEnrichment(tfbs, t.comp, file.pdf="enrich-qqplot.pdf", plot.type="polar", enrichment.type="both")

