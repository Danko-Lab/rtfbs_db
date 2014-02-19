rtfbs_db
========

Parse TF motifs from public databases, read into R, and scan using 'rtfbs'.

WARNING:
--------

This package is in no way intended for general distribution or use... It is UNSUPPORTED and will receive bugfixes and updates only if they are useful to me.

Sources of Position Weight Matrices (PWMs).
-------------------------------------------

* Jolma: Downloaded from the supplementary data of this paper: http://www.cell.com/retrieve/pii/S0092867412014961
	Jolma A, et. al.  DNA-binding specificities of human transcription factors. Cell. 2013 Jan 17;152(1-2):327-39. doi: 10.1016/j.cell.2012.12.009. 

* Jaspar: http://129.177.120.189/cgi-bin/jaspar2010/jaspar_db.pl
	Mathelier A, et. al. JASPAR 2014: an extensively expanded and updated open-access database of transcription factor binding profiles. Nucleic Acids Res. 2014 Jan 1;42(1):D142-7. doi: 10.1093/nar/gkt997. Epub 2013 Nov 4.

* Neph: URL: ftp://ftp.ebi.ac.uk/pub/databases/ensembl/encode/supplementary/integration_data_jan2011/byDataType/footprints/jan2011/de.novo.pwm 
	Neph S, et. al. An expansive human regulatory lexicon encoded in transcription factor footprints. Nature. 2012 Sep 6;489(7414):83-90. doi: 10.1038/nature11212.

Requires
--------

* rtfbs: 
	## In R
	install.package('rtfbs')

* bedops:
	* Get the latest version of the bedops binaries here: https://bedops.readthedocs.org/en/latest/
	* Install, and add them to your path.

* The twoBitToFa program from the Kent libraries.  Download it here: http://hgdownload.cse.ucsc.edu/admin/exe/

* 2bit files for your genome of interest.  Find links to these here: http://hgdownload.cse.ucsc.edu/downloads.html

Usage
-----

See files in create_db/ for examples.
