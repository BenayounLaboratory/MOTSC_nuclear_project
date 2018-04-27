setwd('/Volumes/BB_USC_1/Collaborations/David_Lee_collaboration/HEK293_MOTSC_RNAseq/Comparisons')
options(stringsAsFactors = F)
library(DESeq2)
library(Vennerable)

########################################################################################################################################################################
# ChIPseq

load('2018-03-22_2018-03-22_EV_M3_RNAseq_HEK293_3hGR__statistics.RData')
motsc.rna <- data.frame(res); rm(res)
motsc.rna$GeneName <- rownames(motsc.rna)

my.nrf2.chip <- read.csv('HOMER_NRF2_peaks_in_at_least_3_cells_A549_Lymphoblatoid_HeLa_HepG2_IMR90.xls', sep = "\t", header = T)
my.atf1.chip <- read.csv('HOMER_ATF1_peaks_in_at_least_2_cells_K562_HepG2_LoVo.xls', sep = "\t", header = T)
my.jund.chip <- read.csv('HOMER_JunD_peaks_in_at_least_3_cells_HCT116_HepG2_MCF7_SKNSH_T47D.xls', sep = "\t", header = T)

my.M3.rna.sig <- motsc.rna[motsc.rna$padj < 0.05,] # 802
my.M3.rna_dwn.sig <- my.M3.rna.sig[my.M3.rna.sig$log2FoldChange < 0,]# 412
my.M3.rna_up.sig <- my.M3.rna.sig[my.M3.rna.sig$log2FoldChange > 0,]#390


########################################################################################################################
# Comparison to NRF2 targets

# Up M3 target genes
my.RNA.up <- list("NRF2" = unique(my.nrf2.chip$Gene.Name),
                  "Mots-c" = unique(rownames(my.M3.rna_up.sig)))

write.table(intersect(my.RNA.up[[1]],my.RNA.up[[2]]), file = "NRF2_amd_motsc_up.txt", quote = F, row.names = F, col.names = F)

Venn.RNA.up <- Venn(my.RNA.up)

pdf(paste(Sys.Date(),"Venn_MOTSC_UP_NRF2.pdf", sep = "_"))
plot(Venn.RNA.up, doWeights=T)
dev.off()

# DOWN M3 target genes
my.RNA.dwn <- list("NRF2" = unique(my.nrf2.chip$Gene.Name),
                   "Mots-c" = unique(rownames(my.M3.rna_dwn.sig)))

Venn.RNA.dwn <- Venn(my.RNA.dwn)
write.table(intersect(my.RNA.dwn[[1]],my.RNA.dwn[[2]]), file = "NRF2_amd_motsc_DWN.txt", quote = F, row.names = F, col.names = F)

pdf(paste(Sys.Date(),"Venn_MOTSC_DOWN_NRF2.pdf", sep = "_"))
plot(Venn.RNA.dwn, doWeights=T)
dev.off()


########################################################################################################################
# Comparison to ATF1 targets

# Up M3 target genes
my.RNA.up <- list("ATF1" = unique(my.atf1.chip$Gene.Name),
                  "Mots-c" = unique(rownames(my.M3.rna_up.sig)))

Venn.RNA.up <- Venn(my.RNA.up)

pdf(paste(Sys.Date(),"Venn_MOTSC_UP_ATF1.pdf", sep = "_"))
plot(Venn.RNA.up, doWeights=T)
dev.off()

# DOWN M3 target genes
my.RNA.dwn <- list("ATF1" = unique(my.atf1.chip$Gene.Name),
                   "Mots-c" = unique(rownames(my.M3.rna_dwn.sig)))

Venn.RNA.dwn <- Venn(my.RNA.dwn)

pdf(paste(Sys.Date(),"Venn_MOTSC_DOWN_ATF1.pdf", sep = "_"))
plot(Venn.RNA.dwn, doWeights=T)
dev.off()


########################################################################################################################
# Comparison to JUND targets

# Up M3 target genes
my.RNA.up <- list("JUND" = unique(my.jund.chip$Gene.Name),
                  "Mots-c" = unique(rownames(my.M3.rna_up.sig)))

Venn.RNA.up <- Venn(my.RNA.up)

pdf(paste(Sys.Date(),"Venn_MOTSC_UP_JUND.pdf", sep = "_"))
plot(Venn.RNA.up, doWeights=T)
dev.off()

# DOWN M3 target genes
my.RNA.dwn <- list("JUND" = unique(my.jund.chip$Gene.Name),
                   "Mots-c" = unique(rownames(my.M3.rna_dwn.sig)))

Venn.RNA.dwn <- Venn(my.RNA.dwn)

pdf(paste(Sys.Date(),"Venn_MOTSC_DOWN_JUND.pdf", sep = "_"))
plot(Venn.RNA.dwn, doWeights=T)
dev.off()


########################################################################################################################
# Compare overlap of targets of the 3 ARE-binding TFs

my.chips <- list("JUND" = unique(my.jund.chip$Gene.Name),
                 "NRF2" = unique(my.nrf2.chip$Gene.Name),
                 "ATF1" = unique(my.atf1.chip$Gene.Name))

Venn.ChIPs <- Venn(my.chips)

pdf(paste(Sys.Date(),"Venn_NRF2_JUND_ATF1.pdf", sep = "_"))
plot(Venn.ChIPs, doWeights=T)
dev.off()
