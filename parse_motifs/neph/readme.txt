## Motifs derived from de novo calls in the Neph et. al. ENCODE paper.
wget ftp://ftp.ebi.ac.uk/pub/databases/ensembl/encode/supplementary/integration_data_jan2011/byDataType/footprints/jan2011/de.novo.pwm
mv de.novo.pwm neph.txt
cat neph.txt | perl makeNeph.pl
#
# From the README:
#
#./jan2011/de.novo.pwm shows the 683 discovered de novo motifs within footprints
#  : these are positional weight matrices
#  : columns are A,C,G,T
#

