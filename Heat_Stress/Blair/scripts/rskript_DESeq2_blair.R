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

setwd("C:/Users/souku/job/eneza/heat_treatment/blair/")

#read information from samplesheet
samples <- read.csv("skripts/samplesheet_DESeq2_blair.csv", header=TRUE)

samples$group <- as.character(samples$group)#defines wether the sample is control or heat treatment
samples$replicate <- as.character(samples$replicate)
samples$timepoint <- as.character(samples$timepoint)

#read tx2gene file
tx2gene <- read_tsv("star_salmon/tx2gene.tsv")

######################## ANALYSING TIMEPOINT ZT1 ######################################
#create subset
samples1 <- subset(samples, timepoint == "ZT1")
samples1

#paths to quant.sf
quantFiles1  <- file.path("star_salmon",samples1$sample, "quant.sf")
names(quantFiles1) <- samples1$sample

#import files with tximport
txi1 <- tximport(quantFiles1, type = "salmon", tx2gene = tx2gene)

#create a DESeqDataSet
dds1 <- DESeqDataSetFromTximport(txi = txi1,
                                colData = samples1,
                                design = ~ group) 
dim(dds1)


############## Visualization #################
# pre-filter
smallestGroupSize <- 4 #four biological replicates
keep <- rowSums(counts(dds1) >= 10) >= smallestGroupSize
dds1 <- dds1[keep,]
dim(dds1)
#give factor levels
dds1$group <- relevel(dds1$group, ref = "control")
dds1$group

#transform data for visualization
vsd <- vst(dds1, blind=TRUE) #variance stabilizing transformations
rld <- rlog(dds1, blind=TRUE)
head(assay(vsd), 3)
ntd <- normTransform(dds1)

#comparing the transformations
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld))# this one gave the best flat curve 
meanSdPlot(assay(vsd)) 

png("results/ZT1/meanSdPlot_rld_ZT1_blair.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(rld)) 
dev.off()

#create Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$group)
colnames(sampleDistMatrix) <- NULL

png("results/ZT1/heatmap_ZT1_blair.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()


#create PCA plot
pcaData <- plotPCA(rld, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("results/ZT1/pca_ZT1_blair.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

############## Differential expression analysis #############
##############      with no-normalized data     #############
dds1 <- DESeq(dds1)
res1 <- results(dds1)
res1
write.table(res1, file="results/ZT1/res_table_ZT1_blair.csv", sep="\t", quote=F, col.names=NA)


#exporting normalized counts
normalized_counts <- counts(dds1, normalized=TRUE)
write.table(normalized_counts, file="results/ZT1/normalized_counts_ZT1_blair.txt", sep="\t", quote=F, col.names=NA)


#plot Dispersion Estimate 
png("results/ZT1/dispersion_estimate_ZT1_blair.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds1)
dev.off()

#exploring name of the group
resultsNames(dds1)

res_warm_vs_control_ZT1 <- results(dds1, name="group_warm_vs_control")

#create MA plot to see differential expressed genes
png("results/ZT1/MA_plots_ZT1_blair.png", width = 18, height = 16, units = "cm", res=300)
plotMA(res_warm_vs_control_ZT1, ylim=c(-2.5,2.5), main= "warm_vs_control")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
dev.off()

#set values for log2FoldChange and p-value
l2fc = 1.0
alpha = 0.05

#############################
# compare warm with control #
#############################

all_genes_ZT1_warm_vs_control <- sort(as.character(rownames(res_warm_vs_control_ZT1)), decreasing = TRUE)

#extract significant results
#upregulated
signif_res_warm_vs_control_up_ZT1 <- res_warm_vs_control_ZT1[ which((res_warm_vs_control_ZT1$padj < alpha & !is.na(res_warm_vs_control_ZT1$padj)) & res_warm_vs_control_ZT1$log2FoldChange > l2fc), ]
signif_res_warm_vs_control_up_ZT1 <- sort(as.character(rownames(signif_res_warm_vs_control_up_ZT1)), decreasing = TRUE)
#downregulated
signif_res_warm_vs_control_down_ZT1 <- res_warm_vs_control_ZT1[ which((res_warm_vs_control_ZT1$padj < alpha & !is.na(res_warm_vs_control_ZT1$padj)) & res_warm_vs_control_ZT1$log2FoldChange < -l2fc), ]
signif_res_warm_vs_control_down_ZT1 <- sort(as.character(rownames(signif_res_warm_vs_control_down_ZT1)), decreasing = TRUE)

#GO enrichment analysis for upregulated genes
ego_warm_vs_control_up_ZT1 <- enrichGO(gene = signif_res_warm_vs_control_up_ZT1,
                          universe = all_genes_ZT1_warm_vs_control,
                          OrgDb = "org.At.tair.db",
                          keyType = "TAIR",
                          ont = "BP", #biological process
                          pAdjustMethod = "BH",
                          readable = TRUE)

dotplot(ego_warm_vs_control_up_ZT1, showCategory = 20)
upsetplot(ego_warm_vs_control_up_ZT1, n = 10)

png("results/ZT1/enrich_warm_vs_control_up_ZT1_blair.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_up_ZT1, showCategory = 15, title = "ZT1 warm_vs_control up LFC 1.0")
dev.off()

#more details
write.csv2(as.data.frame(ego_warm_vs_control_up_ZT1), file="results/ZT1/warm_vs_control_up_ZT1_blair.csv")

#GO enrichment analysis for downregulated genes
ego_warm_vs_control_down_ZT1 <- enrichGO(gene = signif_res_warm_vs_control_down_ZT1,
                            universe = all_genes_ZT1_warm_vs_control,
                            OrgDb = "org.At.tair.db",
                            keyType = "TAIR",
                            ont = "BP", #biological process
                            pAdjustMethod = "BH",
                            readable = TRUE)

dotplot(ego_warm_vs_control_down_ZT1, showCategory = 20)
upsetplot(ego_warm_vs_control_down_ZT1, n = 10)

png("results/ZT1/enrich_warm_vs_control_down_ZT1_blair.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_down_ZT1, showCategory = 15, title = "ZT1 warm_vs_control down LFC 1.0")
dev.off()

#more details
write.csv2(as.data.frame(ego_warm_vs_control_down_ZT1), file="results/ZT1/warm_vs_control_down_ZT1_blair.csv")



######################## ANALYSING TIMEPOINT ZT6 ######################################
#create subset
samples2 <- subset(samples, timepoint == "ZT6")
head(samples2)

#paths to quant.sf
quantFiles2  <- file.path("star_salmon",samples2$sample, "quant.sf")
names(quantFiles2) <- samples2$sample

#import files with tximport
txi2 <- tximport(quantFiles2, type = "salmon", tx2gene = tx2gene)

#create a DESeqDataSet
dds2 <- DESeqDataSetFromTximport(txi = txi2,
                                 colData = samples2,
                                 design = ~ group) 
dim(dds2)

############## Visualization #################
# pre-filter 
smallestGroupSize <- 2 #two biological replicates
keep <- rowSums(counts(dds2) >= 10) >= smallestGroupSize
dds2 <- dds2[keep,]
dim(dds2)

#give factor levels
dds2$group <- relevel(dds2$group, ref = "control")
dds2$group

#transform data for visualization
vsd <- vst(dds2, blind=TRUE) #variance stabilizing transformations
rld <- rlog(dds2, blind=TRUE)
head(assay(vsd), 3)
ntd <- normTransform(dds2)

#comparing the transformations
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld)) 
meanSdPlot(assay(vsd))# this one gave the best flat curve 

png("results/ZT6/meanSdPlot_vsd_ZT6_blair.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(vsd)) 
dev.off()

#create Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group)
colnames(sampleDistMatrix) <- NULL

png("results/ZT6/heatmap_ZT6_bliar.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()


#create PCA plot
pcaData <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("results/ZT6/pca_ZT6_bliar.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

############## Differential expression analysis #############
##############      with no-normalized data     #############
dds2 <- DESeq(dds2)
res2 <- results(dds2)
res2
write.table(res2, file="results/ZT6/res_table_ZT6_blair.csv", sep="\t", quote=F, col.names=NA)


#export normalized counts
normalized_counts <- counts(dds2, normalized=TRUE)
write.table(normalized_counts, file="results/ZT6/normalized_counts_ZT6_blair.txt", sep="\t", quote=F, col.names=NA)


#plot Dispersion Estimate 
png("results/ZT6/dispersion_estimate_ZT6_blair.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds2)
dev.off()

#exploring name of group
resultsNames(dds2)

res_warm_vs_control_ZT6 <- results(dds2, name="group_warm_vs_control")

#create MA plot to see differential expressed genes
png("results/ZT6/MA_plots_ZT6_blair.png", width = 18, height = 16, units = "cm", res=300)
plotMA(res_warm_vs_control_ZT6, ylim=c(-2.5,2.5), main= "warm_vs_control")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
dev.off()

#set values for log2FoldChange and p-value (same as for ZT1)
l2fc = 1.0
alpha = 0.05

#############################
# compare warm with control #
#############################

all_genes_ZT6_warm_vs_control <- sort(as.character(rownames(res_warm_vs_control_ZT6)), decreasing = TRUE)

#extract significant results
#upregulated
signif_res_warm_vs_control_up_ZT6 <- res_warm_vs_control_ZT6[ which((res_warm_vs_control_ZT6$padj < alpha & !is.na(res_warm_vs_control_ZT6$padj)) & res_warm_vs_control_ZT6$log2FoldChange > l2fc), ]
signif_res_warm_vs_control_up_ZT6 <- sort(as.character(rownames(signif_res_warm_vs_control_up_ZT6)), decreasing = TRUE)
#downregulated
signif_res_warm_vs_control_down_ZT6 <- res_warm_vs_control_ZT6[ which((res_warm_vs_control_ZT6$padj < alpha & !is.na(res_warm_vs_control_ZT6$padj)) & res_warm_vs_control_ZT6$log2FoldChange < -l2fc), ]
signif_res_warm_vs_control_down_ZT6 <- sort(as.character(rownames(signif_res_warm_vs_control_down_ZT6)), decreasing = TRUE)

#GO enrichment analysis for upregulated genes
ego_warm_vs_control_up_ZT6 <- enrichGO(gene = signif_res_warm_vs_control_up_ZT6,
                                      universe = all_genes_ZT6_warm_vs_control,
                                      OrgDb = "org.At.tair.db",
                                      keyType = "TAIR",
                                      ont = "BP", #biological process
                                      pAdjustMethod = "BH",
                                      readable = TRUE)

dotplot(ego_warm_vs_control_up_ZT6, showCategory = 20)
upsetplot(ego_warm_vs_control_up_ZT6, n = 10)

png("results/ZT6/enrich_warm_vs_control_up_ZT6_blair.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_up_ZT6, showCategory = 15, title = "ZT6 warm_vs_control up LFC 1.0")
dev.off()

#more detail
write.csv2(as.data.frame(ego_warm_vs_control_up_ZT6), file="results/ZT6/warm_vs_control_up_ZT6_blair.csv")

#GO enrichment analysis for downregulated genes
ego_warm_vs_control_down_ZT6 <- enrichGO(gene = signif_res_warm_vs_control_down_ZT6,
                                        universe = all_genes_ZT6_warm_vs_control,
                                        OrgDb = "org.At.tair.db",
                                        keyType = "TAIR",
                                        ont = "BP", #biological process
                                        pAdjustMethod = "BH",
                                        readable = TRUE)

dotplot(ego_warm_vs_control_down_ZT6, showCategory = 20)
upsetplot(ego_warm_vs_control_down_ZT6, n = 10)

png("results/ZT6/enrich_warm_vs_control_down_ZT6_blair.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_down_ZT6, showCategory = 15, title = "ZT6 warm_vs_control down LFC 1.0")
dev.off()

#more detail
write.csv2(as.data.frame(ego_warm_vs_control_down_ZT6), file="results/ZT6/warm_vs_control_down_ZT6_blair.csv")


#########################################################
# compare results with genes predicted by the ML models #
#########################################################

#put all different expressed genes together
all_up_genes <- union(signif_res_warm_vs_control_up_ZT1, signif_res_warm_vs_control_up_ZT6)
all_down_genes <- union(signif_res_warm_vs_control_down_ZT1, signif_res_warm_vs_control_down_ZT6)
diff_genes <- union(all_up_genes, all_down_genes)
write.table(diff_genes, file="results/diff_reg_genes_blair.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

#genes predicted by the different ML models
data1 <- read.csv2("rfe_heat.csv", header = TRUE, sep = ";")
data2 <- read.csv2("lasso_heat.csv", header = TRUE, sep = ";")
data3 <- read.csv2("xgb-full_heat.csv", header = TRUE, sep = ";")
data4 <- read.csv2("overlaps-lasso-rfe_heat.csv", header = TRUE, sep = ";")

##############
# RFE model  #
##############
rfe <- data1$genes
length(rfe)

#upregulated
rfe_up <- intersect(rfe, all_up_genes)
length(rfe_up)

#downregulated
rfe_down <- intersect(rfe, all_down_genes)
length(rfe_down)

#general different regulated
rfe_general <- union(rfe_up, rfe_down)
length(rfe_general)

print(rfe_general)
write.table(rfe_general, file="results/rfe_genes_blair.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)


################
# LASSO model  #
################
lasso <- data2$genes
length(lasso)

#upregulated
lasso_up <- intersect(lasso, all_up_genes)
length(lasso_up)

#downregulated
lasso_down <- intersect(lasso, all_down_genes)
length(lasso_down)

#general different regulated
lasso_general <- union(lasso_up, lasso_down)
length(lasso_general)

print(lasso_general)
write.table(lasso_general, file="results/lasso_genes_blair.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

##############
# FULL model #
##############
full <- data3$genes
length(full)

#upregulated
full_up <- intersect(full, all_up_genes)
length(full_up)

#downregulated
full_down <- intersect(full, all_down_genes)
length(full_down)

#general different regulated
full_general <- union(full_up, full_down)
length(full_general)

print(full_general)
write.table(full_general, file="results/full_genes_blair.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)



##################
# OVERLAP model  #
##################
over <- data4$genes
length(over)

#upregulated
over_up <- intersect(over, all_up_genes)
length(over_up)

#downregulated
over_down <- intersect(over, all_down_genes)
length(over_down)

#general different regulated
over_general <- union(over_up, over_down)
length(over_general)

print(over_general)
write.table(over_general, file="results/overlap_genes_blair.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)



########################
# calculate percentage #
########################

#RFE
percentage_rfe_general <- (length(rfe_general)/length(rfe))*100
print(percentage_rfe_general)

#LASSO
percentage_lasso_general <- (length(lasso_general)/length(lasso))*100
print(percentage_lasso_general)

#FULL
percentage_full_general <- (length(full_general)/length(full))*100
print(percentage_full_general)

#OVERLAP
percentage_over_general <- (length(over_general)/length(over))*100
print(percentage_over_general)


###########
# Plotten #
###########
#only general different regulated
names <- c("Lasso", "RFE", "Overlap", "Full")
values <- c(percentage_lasso_general, percentage_rfe_general, percentage_over_general, percentage_full_general)


png("results/models_diff_regulated_blair.png", width = 18, height = 16, units = "cm", res = 300)
h <- barplot(values, xlab = "Model", ylab = "Found Genes in %",
             main=paste("Comparison of the ML-models (p-value=",alpha," and l2fc=",l2fc,")\nNumber of genes found in Blair et al., 2019: ", length(diff_genes), sep=""), cex.main=1, col="lightblue", names.arg=names)
text(h, c(4), labels = paste(round(values, digits=1), "%"), pos=3)
dev.off()
