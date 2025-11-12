#load packages
library(DESeq2)
library(readr)
library(tximport)
library(vsn)
library(pheatmap)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(org.At.tair.db)

#read information from samplesheet
samples <- read.csv("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/skripts/samlpesheet_DESeq2.csv", header=TRUE)

samples$group <- as.character(samples$group)
samples$timepoint <- as.character(samples$timepoint)#0h -> control, 3h -> 3h cold treatment, 24h -> 24h cold treatment
samples$repetition <- as.character(samples$repetition)


#paths to quant.sf
quantFiles <- file.path("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/star_salmon", samples$sample, "quant.sf")
names(quantFiles) <- samples$sample

#read tx2gene file 
tx2gene <- read_tsv("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/star_salmon/tx2gene.tsv")

#import files with tximport
txi <- tximport(quantFiles, type = "salmon", tx2gene = tx2gene)

#create a DESeqDataSet
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ timepoint) 
dim(dds)

############## Visualization #################

#pre-filter
smallestGroupSize <- 2 #two independent repetitions of the experiment
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds)

#give factor levels
dds$timepoint <- relevel(dds$timepoint, ref = "0h")
dds$timepoint

#transform data for visualization
vsd <- vst(dds, blind=TRUE) #variance stabilizing transformations
rld <- rlog(dds, blind=TRUE) #regularized logarithm
head(assay(vsd), 3)
ntd <- normTransform(dds)

#comparing the transformations
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld))
meanSdPlot(assay(vsd)) # this one gave the best flat curve 

png("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/meanSdPlot_vsd_jia.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(vsd)) 
dev.off()

#create Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$timepoint)
colnames(sampleDistMatrix) <- NULL

png("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/heatmap_jia.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()


#create PCA plot
pcaData <- plotPCA(vsd, intgroup=c("timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/pca_jia.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData, aes(PC1, PC2, color=timepoint)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()


############## Differential expression analysis #############
##############      with no-normalized data     #############
dds <- DESeq(dds)
res <- results(dds)
res
write.table(res, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/res_table_jia.csv", sep="\t", quote=F, col.names=NA)

#exporting normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/normalized_counts_jia.txt", sep="\t", quote=F, col.names=NA)


#plot Dispersion Estimate 
png("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/dispersion_estimate_jia.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds)
dev.off()

#exploring names of the groups based on timepoint
resultsNames(dds)

#define groups to compare
res_24_vs_0 <- results(dds, name="timepoint_24h_vs_0h")
res_3_vs_0 <- results(dds, name="timepoint_3h_vs_0h")

#create MA plot to see differential expressed genes
png("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/MA_plots_jia.png", width = 18, height = 16, units = "cm", res=300)
par(mfrow = c(2, 2), mar=c(2,2,2,2))
plotMA(res_24_vs_0, ylim=c(-2.5,2.5), main= "24_vs_0")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
plotMA(res_3_vs_0, ylim=c(-2.5,2.5), main= "3_vs_0")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
dev.off()

#set values of log2FoldChange and p-value
l2fc = 1.0
alpha = 0.05

#########################################
# compare 24h with 0h of cold treatment #
#########################################

all_genes_24_vs_0<- sort(as.character(rownames(res_24_vs_0)), decreasing = TRUE)

#extract significant results
##upregulated
signif_res_24_vs_0_up <- res_24_vs_0[ which((res_24_vs_0$padj < alpha & !is.na(res_24_vs_0$padj)) & res_24_vs_0$log2FoldChange > l2fc), ]
signif_res_24_vs_0_up <- sort(as.character(rownames(signif_res_24_vs_0_up)), decreasing = TRUE)
#downregulated
signif_res_24_vs_0_down <- res_24_vs_0[ which((res_24_vs_0$padj < alpha & !is.na(res_24_vs_0$padj)) & res_24_vs_0$log2FoldChange < -l2fc), ]
signif_res_24_vs_0_down <- sort(as.character(rownames(signif_res_24_vs_0_down)), decreasing = TRUE)

#GO enrichment analysis for upreguated genes
ego_24_vs_0_up <- enrichGO(gene = signif_res_24_vs_0_up,
                          universe = all_genes_24_vs_0,
                          OrgDb = "org.At.tair.db",
                          keyType = "TAIR",
                          ont = "BP", #biological process
                          pAdjustMethod = "BH",
                          readable = TRUE)

dotplot(ego_24_vs_0_up, showCategory = 20)
upsetplot(ego_24_vs_0_up, n = 10)

png("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/enrich_24_vs_0_up_jia.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_24_vs_0_up, showCategory = 15, title = "group 24_vs_0 up LFC 1.0")
dev.off()

#for more detail
write.csv2(as.data.frame(ego_24_vs_0_up), file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/24_vs_0_up_jia.csv")

#GO enrichment analysis for downregulated genes
ego_24_vs_0_down <- enrichGO(gene = signif_res_24_vs_0_down,
                            universe = all_genes_24_vs_0,
                            OrgDb = "org.At.tair.db",
                            keyType = "TAIR",
                            ont = "BP", #biological process
                            pAdjustMethod = "BH",
                            readable = TRUE)

dotplot(ego_24_vs_0_down, showCategory = 20)
upsetplot(ego_24_vs_0_down, n = 10)

png("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/enrich_24_vs_0_down_jia.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_24_vs_0_down, showCategory = 15, title = "group 24_vs_0 down LFC 1.0")
dev.off()

#for more detail
write.csv2(as.data.frame(ego_24_vs_0_down), file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/24_vs_0_down_jia.csv")


########################################
# compare 3h with 0h of cold treatment #
########################################

all_genes_3_vs_0 <- sort(as.character(rownames(res_3_vs_0)), decreasing = TRUE)

#extract significant results
#upregulated
signif_res_3_vs_0_up <- res_3_vs_0[ which((res_3_vs_0$padj < alpha & !is.na(res_3_vs_0$padj)) & res_3_vs_0$log2FoldChange > l2fc), ]
signif_res_3_vs_0_up <- sort(as.character(rownames(signif_res_3_vs_0_up)), decreasing = TRUE)
#downregulated
signif_res_3_vs_0_down <- res_3_vs_0[ which((res_3_vs_0$padj < alpha & !is.na(res_3_vs_0$padj)) & res_3_vs_0$log2FoldChange < -l2fc), ]
signif_res_3_vs_0_down <- sort(as.character(rownames(signif_res_3_vs_0_down)), decreasing = TRUE)

#GO enrichment analysis for upregulated genes
ego_3_vs_0_up <- enrichGO(gene = signif_res_3_vs_0_up,
                          universe = all_genes_3_vs_0,
                          OrgDb = "org.At.tair.db",
                          keyType = "TAIR",
                          ont = "BP", #biological process
                          pAdjustMethod = "BH",
                          readable = TRUE)

dotplot(ego_3_vs_0_up, showCategory = 20)
upsetplot(ego_3_vs_0_up, n = 10)

png("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/enrich_3_vs_0_up_jia.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_3_vs_0_up, showCategory = 15, title = "3_vs_0 up LFC 1.0")
dev.off()

#locking for more detail
write.csv2(as.data.frame(ego_3_vs_0_up), file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/3_vs_0_up_jia.csv")

#GO enrichment analysis for downregulated genes
ego_3_vs_0_down <- enrichGO(gene = signif_res_3_vs_0_down,
                            universe = all_genes_3_vs_0,
                            OrgDb = "org.At.tair.db",
                            keyType = "TAIR",
                            ont = "BP", #biological process
                            pAdjustMethod = "BH",
                            readable = TRUE)

dotplot(ego_3_vs_0_down, showCategory = 20)
upsetplot(ego_3_vs_0_down, n = 10)

png("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/enrich_3_vs_0_down_jia.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_3_vs_0_down, showCategory = 15, title = "3_vs_0 down LFC 1.0")
dev.off()

#looking for more detail
write.csv2(as.data.frame(ego_3_vs_0_down), file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/3_vs_0_down_jia.csv")





#########################################################
# compare results with genes predicted by the ML models #
#########################################################

#extract all differential expressed genes found in the publication 
all_up_genes <- union(signif_res_24_vs_0_up, signif_res_3_vs_0_up)
all_down_genes <- union(signif_res_24_vs_0_down, signif_res_3_vs_0_down)
diff_genes <- union(all_up_genes, all_down_genes)
write.table(diff_genes, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/all_diff_reg_genes_jia.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)


##############
# RFE model  #
##############
#genes predicted by the RFE model
data1 <- read.csv2("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/rfe model - Important.csv", header = TRUE, sep = ",")

#extract the genes
rfe <- data1$Genes
length(rfe)

#intersection of the upregulated genes
rfe_up <- intersect(rfe, all_up_genes)
length(rfe_up)

print(rfe_up)
write.table(rfe_up, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/rfe_up_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#intersection of the downregulated genes
rfe_down <- intersect(rfe, all_down_genes)
length(rfe_down)

print(rfe_down)
write.table(rfe_down, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/rfe_down_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#general different regulated
rfe_general <- union(rfe_up, rfe_down)
length(rfe_general)

print(rfe_general)
write.table(rfe_general, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/rfe_general_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


##############
# XGB model  #
##############
#genes predicted by the XGB model
data2 <- read.csv2("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/xgb full - Important.csv", header = TRUE, sep = ",")

xgb <- data2$Genes
length(xgb)

#upregulated
xgb_up <- intersect(xgb, all_up_genes)
length(xgb_up)

print(xgb_up)
write.table(xgb_up, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/xgb_up_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#downregulated
xgb_down <- intersect(xgb, all_down_genes)
length(xgb_down)

print(xgb_down)
write.table(xgb_down, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/xgb_down_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#general different regulated
xgb_general <- union(xgb_up, xgb_down)
length(xgb_general)

print(xgb_general)
write.table(xgb_general, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/xgb_general_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


################
# Lasso model  #
################
#genes predicted by the Lasso model
data3 <- read.csv2("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/lasso model - Important.csv", header = TRUE, sep = ",")

lasso <- data3$Genes
length(lasso)

#upregulated
lasso_up <- intersect(lasso, all_up_genes)
length(lasso_up)

print(lasso_up)
write.table(lasso_up, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/lasso_up_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#downregulated
lasso_down <- intersect(lasso, all_down_genes)
length(lasso_down)

print(lasso_down)
write.table(lasso_down, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/lasso_down_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#general different regulated
lasso_general <- union(lasso_up, lasso_down)
length(lasso_general)

print(lasso_general)
write.table(lasso_general, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/lasso_general_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)



#####################
# Overlap Lasso RFE #
#####################
#genes predicted by Lasso and RFE model
data4 <- read.csv2("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/overlaps - Lasso RFE models.csv", header = TRUE, sep = ",")

over <- data4$Genes
length(over)

#upregulated
over_up <- intersect(over, all_up_genes)
length(over_up)

print(over_up)
write.table(over_up, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/overlap_up_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#downregulated
over_down <- intersect(over, all_down_genes)
length(over_down)

print(over_down)
write.table(over_down, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/overlap_down_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#general different regulated
over_general <- union(over_up, over_down)
length(over_general)

print(over_general)
write.table(over_general, file="C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/overlap_general_genes_jia.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


########################
# calculate percentage #
########################

#RFE
percentage_rfe_general <- (length(rfe_general)/length(rfe))*100
print(percentage_rfe_general)

#XGB
percentage_xgb_general <- (length(xgb_general)/length(xgb))*100
print(percentage_xgb_general)

#LASSO
percentage_lasso_general <- (length(lasso_general)/length(lasso))*100
print(percentage_lasso_general)

#OVERLAP
percentage_over_general <- (length(over_general)/length(over))*100
print(percentage_over_general)


####################
# plot percentages #
####################
#only general different regulated
names <- c("RFE", "XGB", "Lasso", "Overlap")
values <- c(percentage_rfe_general, percentage_xgb_general, percentage_lasso_general, percentage_over_general)


png("C:/Users/souku/job/eneza/cold_treatment/arabid_jia/results/all_models_diff_regulated_jia.png", width = 18, height = 16, units = "cm", res = 300)
h <- barplot(values, xlab = "Model", ylab = "Found Genes in %",
             main=paste("Comparison of the ML-models (p-value=",alpha," and l2fc=",l2fc,")\nNumber of genes found in Jia et al., 2016: ", length(diff_genes), sep=""), cex.main=1, col="lightblue", names.arg=names)
text(h, c(4), labels = paste(round(values, digits=1), "%"), pos=3)
dev.off()
