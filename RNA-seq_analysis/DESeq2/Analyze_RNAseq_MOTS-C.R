setwd('/Volumes/BB_USC_1/Collaborations/David_Lee_collaboration/HEK293_MOTSC_RNAseq/DEseq2/')
options(stringsAsFactors = F)

library(DESeq2)
library(pheatmap)
library('pvclust')

# 2018-03-22
# analyze HEK293 MOTSC/Glucose restriction RNAseq

##############################################################################################################
# read in Kallisto mappings
my.data <- read.csv("../Kallisto/2018-03-22_HEK293_EV_M3_GR_kallisto_mapping.txt", sep = "\t", header = T,row.names = NULL)

# sum read over genes (to not have results over transcripts for DEseq2)
my.data.per.gene <- aggregate(my.data[,5:28],by=list(my.data$GeneSymbol),FUN=sum)

# round counts (DESeq needs integers)
my.data.per.gene[,2:25] <- round(my.data.per.gene[,2:25])
rownames(my.data.per.gene) <- my.data.per.gene$Group.1

# get the genes with no reads out
my.null <- which(apply(my.data.per.gene[,2:25], 1, sum) <= 5) # see deseq2 vignette, need to remove too low genes; 6774 removed
my.filtered.matrix <- my.data.per.gene[-my.null,2:25] # 21039 genes

###### Analyze DE *only* upon 3h of GR
my.outprefix <- paste(Sys.Date(),"EV_M3_RNAseq_HEK293_3hGR", sep = "_")

# get covariates
my.m3.expression <- c(rep("EV",6),rep("M3",6))

# design matrix
dataDesign = data.frame( row.names = colnames( my.filtered.matrix[,c(7:12,19:24)] ), 
                         motsc = my.m3.expression)

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix[,c(7:12,19:24)],
                              colData = dataDesign,
                              design = ~ motsc)

# run DESeq normalizations and export results
dds.deseq <- DESeq(dds)

res <- results(dds.deseq, contrast=c("motsc","M3","EV")) # added the name of the tested variable and levels

# plot dispersion
my.disp.out <- paste(Sys.Date(),my.outprefix,"_dispersion_plot.pdf")

pdf(my.disp.out)
plotDispEsts(dds.deseq)
dev.off()

# parse sample names
my.sample.names <- unlist(strsplit(colnames( my.filtered.matrix ), c(".kallisto_res.abundance.tsv")))

# normalized expression value
tissue.cts <- log2( counts(dds.deseq, normalize = TRUE) + 0.01)
colnames(tissue.cts) <- my.sample.names

# do MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.colors <- c(rep("dodgerblue2",6),
               rep("red4",6))

my.mds.out <- paste(Sys.Date(),my.outprefix,"MDS_plot.pdf", sep ="_")

pdf(my.mds.out)
plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2",main="Multi-dimensional Scaling",cex=2)
points(x, y, pch=16,col=my.colors,cex=2)
legend("bottomleft",c("EV_3h","M3_3h"),col=c("dodgerblue2","red4"),pch=16,bty='n',pt.cex=2)
dev.off()

### get the heatmap of aging changes at FDR5
## exclude NA
res <- res[!is.na(res$padj),]

genes.m3 <- rownames(res)[res$padj < 0.05]
my.num.m3 <- length(genes.m3) # 802

# heatmap drawing - only if there is at least one gene
my.heatmap.out <- paste(Sys.Date(),my.outprefix,"_Heatmap_significant_genes.pdf", sep = "_")

pdf(my.heatmap.out, width = 10, height = 5, onefile = F)
my.heatmap.title <- paste("M3 significant (FDR<5%), ",my.num.m3, " genes",sep="")
pheatmap(tissue.cts[genes.m3,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 20)
dev.off()

# output result tables to files
my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix.txt", sep = "_")
my.out.stats <- paste(my.outprefix,"_all_genes_statistics.txt", sep = "_")
my.out.fdr5 <- paste(my.outprefix,"_FDR5_genes_statistics.txt", sep = "_")
my.out.rdata <- paste(my.outprefix,"_statistics.RData", sep = "_")

write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
write.table(res, file = my.out.stats , sep = "\t" , row.names = T, quote=F)
write.table(res[genes.m3,], file = my.out.fdr5, sep = "\t" , row.names = T, quote=F)
save(res,file=my.out.rdata)