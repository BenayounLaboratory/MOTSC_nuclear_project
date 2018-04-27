cat ../DEseq2/2018-03-22_2018-03-22_EV_M3_RNAseq_HEK293_3hGR__FDR5_genes_statistics.txt | perl -lane 'print if ($F[2] < 0)' | cut -f 1 | sort -u  > 2018-03-22_EV_M3_RNAseq_HEK293_3hGR.results.FDR5_DOWN.txt
cat ../DEseq2/2018-03-22_2018-03-22_EV_M3_RNAseq_HEK293_3hGR__FDR5_genes_statistics.txt  | perl -lane 'print if ($F[2] > 0)' | cut -f 1 | sort -u > 2018-03-22_EV_M3_RNAseq_HEK293_3hGR.results.FDR5_UP.txt

findMotifs.pl 2018-03-22_EV_M3_RNAseq_HEK293_3hGR.results.FDR5_DOWN.txt human ./Motifs_MOTS-C_GR_DOWN/ -start -500 -end 100 -p 4 -nogo
findMotifs.pl 2018-03-22_EV_M3_RNAseq_HEK293_3hGR.results.FDR5_UP.txt human ./Motifs_MOTS-C_GR_UP/ -start -500 -end 100 -p 4 -nogo
