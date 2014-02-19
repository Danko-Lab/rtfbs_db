# Obtained from: http://129.177.120.189/cgi-bin/jaspar2010/jaspar_db.pl
# Downloaded entire DB.
# Options for extracting: 
rm ../*.pwm
cat matrix_list.txt | perl makeJASPAR.pl ## ALL MOTIFS.
grep vertebrates matrix_list.txt | perl makeJASPAR.pl ## Just vertebrate motifs.
grep insects matrix_list.txt | perl makeJASPAR.pl ## Just insect motifs.
