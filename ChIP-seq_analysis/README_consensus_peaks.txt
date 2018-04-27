# Getting consensus peaks for ATF1, JunD and NRF2
# annotatePeaks.pl is part of the HOMER software suite


# ATF1 consensus peak and annotation
multiIntersectBed -header -i K562_ATF1_peaks.bed K562_HepG2_peaks.bed LoVo_ATF1_peaks.bed | perl -lane 'print if $F[3]>= 3' > ATF1_peaks_in_at_least_2_cells_K562_HepG2_LoVo.bed
annotatePeaks.pl ATF1_peaks_in_at_least_2_cells_K562_HepG2_LoVo.bed hg19 > HOMER_ATF1_peaks_in_at_least_2_cells_K562_HepG2_LoVo.xls

# JunD consensus peak and annotation
multiIntersectBed -header -i HudsonAlpha_ChipSeq_HCT116_JunD_peaks.bed HudsonAlpha_ChipSeq_HepG2_JunD_peaks.bed HudsonAlpha_ChipSeq_MCF7_JunD_peaks.bed HudsonAlpha_ChipSeq_SKNSH_JunD_peaks.bed HudsonAlpha_ChipSeq_T47D_JunD_peaks.bed | perl -lane 'print if $F[3]>= 3' > JunD_peaks_in_at_least_3_cells_HCT116_HepG2_MCF7_SKNSH_T47D.bed
annotatePeaks.pl JunD_peaks_in_at_least_3_cells_HCT116_HepG2_MCF7_SKNSH_T47D.bed hg19 > HOMER_JunD_peaks_in_at_least_3_cells_HCT116_HepG2_MCF7_SKNSH_T47D.xls

# NRF2 consensus peak and annotation
multiIntersectBed -header -i A549_NRF2_peaks.bed Lymphoblastoid_NRF2_peaks.bed NRF2_ChIP-seq/HeLa-S3_NRF2_peaks.bed HepG2_NRF2_peaks.bed NRF2_ChIP-seq/IMR-90_NRF2_peaks.bed | perl -lane 'print if $F[3]>= 3' > NRF2_peaks_in_at_least_3_cells_A549_Lymphoblatoid_HeLa_HepG2_IMR90.bed
annotatePeaks.pl NRF2_peaks_in_at_least_3_cells_A549_Lymphoblatoid_HeLa_HepG2_IMR90.bed hg19 > HOMER_NRF2_peaks_in_at_least_3_cells_A549_Lymphoblatoid_HeLa_HepG2_IMR90.xls

