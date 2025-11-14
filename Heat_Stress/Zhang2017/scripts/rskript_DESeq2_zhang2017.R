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

setwd("C:/Users/souku/job/eneza/heat_treatment/zhang/")

#read information from samplesheet
samples <- read.csv("skripts/samplesheet_DESeq2_zhang.csv", header=TRUE)

samples$group <- as.character(samples$group)#control or heat treatment
samples$repetition <- as.character(samples$repetition)
samples$type <- as.character(samples$type)#rosette leaves (RL), early flowers(EF) or late flowers(LF)

#read tx2gene file
tx2gene <- read_tsv("star_salmon/tx2gene.tsv")

######################## ANALYSING GROUP ROSETTE LEAVES ######################################
#create subset
samples1 <- subset(samples, type == "RL")
head(samples1)

#paths to quant.sf
quantFiles1  <- file.path("star_salmon",samples1$sample, "quant.sf")
names(quantFiles1) <- samples1$sample

#import files with tximport
txi1 <- tximport(quantFiles1, type = "salmon", tx2gene = tx2gene)

#create DESeqDataSet
dds1 <- DESeqDataSetFromTximport(txi = txi1,
                                colData = samples1,
                                design = ~ group) 
dim(dds1)

############## Visualization #################
# pre-filter 
smallestGroupSize <- 3 #three biological replicates
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

png("results/RL/meanSdPlot_rld_RL_zhang.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(rld)) 
dev.off()

#create Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$group)
colnames(sampleDistMatrix) <- NULL

png("results/RL/heatmap_RL_zhang.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()


#create PCA plot
pcaData <- plotPCA(rld, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("results/RL/pca_RL_zhang.png", width = 18, height = 16, units = "cm", res=300)
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
write.table(res1, file="results/RL/res_table_RL_zhang.csv", sep="\t", quote=F, col.names=NA)


#exporting normalized counts
normalized_counts <- counts(dds1, normalized=TRUE)
write.table(normalized_counts, file="results/RL/normalized_counts_RL_zhang.txt", sep="\t", quote=F, col.names=NA)


#plot Dispersion Estimate 
png("results/RL/dispersion_estimate_RL_zhang.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds1)
dev.off()

#exploring name of group
resultsNames(dds1)

res_warm_vs_control_RL <- results(dds1, name="group_warm_vs_control")

#create MA plot to see differential expressed genes
png("results/RL/MA_plots_RL_zhang.png", width = 18, height = 16, units = "cm", res=300)
plotMA(res_warm_vs_control_RL, ylim=c(-2.5,2.5), main= "warm_vs_control")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
dev.off()

#set values for log2FoldChange and p-value
l2fc = 1.0
alpha = 0.05

#############################
# compare warm with control #
#############################

all_genes_RL_warm_vs_control <- sort(as.character(rownames(res_warm_vs_control_RL)), decreasing = TRUE)

#extract significant results
#upregulated
signif_res_warm_vs_control_up_RL <- res_warm_vs_control_RL[ which((res_warm_vs_control_RL$padj < alpha & !is.na(res_warm_vs_control_RL$padj)) & res_warm_vs_control_RL$log2FoldChange > l2fc), ]
signif_res_warm_vs_control_up_RL <- sort(as.character(rownames(signif_res_warm_vs_control_up_RL)), decreasing = TRUE)
#downregulated
signif_res_warm_vs_control_down_RL <- res_warm_vs_control_RL[ which((res_warm_vs_control_RL$padj < alpha & !is.na(res_warm_vs_control_RL$padj)) & res_warm_vs_control_RL$log2FoldChange < -l2fc), ]
signif_res_warm_vs_control_down_RL <- sort(as.character(rownames(signif_res_warm_vs_control_down_RL)), decreasing = TRUE)

#GO enrichment analysis for upregulated genes
ego_warm_vs_control_up_RL <- enrichGO(gene = signif_res_warm_vs_control_up_RL,
                          universe = all_genes_RL_warm_vs_control,
                          OrgDb = "org.At.tair.db",
                          keyType = "TAIR",
                          ont = "BP", #biological process
                          pAdjustMethod = "BH",
                          readable = TRUE)

dotplot(ego_warm_vs_control_up_RL, showCategory = 20)
upsetplot(ego_warm_vs_control_up_RL, n = 10)

png("results/RL/enrich_warm_vs_control_up_RL_zhang.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_up_RL, showCategory = 15, title = "RL warm_vs_control up LFC 1.0")
dev.off()

#more detail
write.csv2(as.data.frame(ego_warm_vs_control_up_RL), file="results/RL/warm_vs_control_up_RL_zhang.csv")

#GO enrichment analysis for downregulated genes
ego_warm_vs_control_down_RL <- enrichGO(gene = signif_res_warm_vs_control_down_RL,
                            universe = all_genes_RL_warm_vs_control,
                            OrgDb = "org.At.tair.db",
                            keyType = "TAIR",
                            ont = "BP", #biological process
                            pAdjustMethod = "BH",
                            readable = TRUE)

dotplot(ego_warm_vs_control_down_RL, showCategory = 20)
upsetplot(ego_warm_vs_control_down_RL, n = 10)

png("results/RL/enrich_warm_vs_control_down_RL_zhang.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_down_RL, showCategory = 15, title = "RL warm_vs_control down LFC 1.0")
dev.off()

#more detail
write.csv2(as.data.frame(ego_warm_vs_control_down_RL), file="results/RL/warm_vs_control_down_RL_zhang.csv")


######################## ANALYSING GROUP EARLY FLOWERS ######################################
#create subset
samples2 <- subset(samples, type == "EF")
head(samples2)

#paths to quant.sf
quantFiles2  <- file.path("star_salmon",samples2$sample, "quant.sf")
names(quantFiles2) <- samples2$sample

#read files with tximport
txi2 <- tximport(quantFiles2, type = "salmon", tx2gene = tx2gene)

#create DESeqDataSet
dds2 <- DESeqDataSetFromTximport(txi = txi2,
                                 colData = samples2,
                                 design = ~ group) 
dim(dds2)

############## Visualization #################
# pre-filter
smallestGroupSize <- 3
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
meanSdPlot(assay(vsd)) # this one gave the best flat curve 

png("results/EF/meanSdPlot_vsd_EF_zhang.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(vsd)) 
dev.off()

#create Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group)
colnames(sampleDistMatrix) <- NULL

png("results/EF/heatmap_EF_zhang.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()


#create PCA plot
pcaData <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("results/EF/pca_EF_zhang.png", width = 18, height = 16, units = "cm", res=300)
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
write.table(res2, file="results/EF/res_table_EF_zhang.csv", sep="\t", quote=F, col.names=NA)


#exporting normalized counts
normalized_counts <- counts(dds2, normalized=TRUE)
write.table(normalized_counts, file="results/EF/normalized_counts_EF_zhang.txt", sep="\t", quote=F, col.names=NA)


#plot Dispersion Estimate 
png("results/EF/dispersion_estimate_EF_zhang.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds2)
dev.off()

#exploring name of the group
resultsNames(dds2)

res_warm_vs_control_EF <- results(dds2, name="group_warm_vs_control")

#create MA plot to see differential expressed genes
png("results/EF/MA_plots_EF_zhang.png", width = 18, height = 16, units = "cm", res=300)
plotMA(res_warm_vs_control_EF, ylim=c(-2.5,2.5), main= "warm_vs_control")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
dev.off()

#set values for log2FoldChange and p-value
l2fc = 1.0
alpha = 0.05

#############################
# compare warm with control #
#############################

all_genes_EF_warm_vs_control <- sort(as.character(rownames(res_warm_vs_control_EF)), decreasing = TRUE)

#extract significant results
#upregulated
signif_res_warm_vs_control_up_EF <- res_warm_vs_control_EF[ which((res_warm_vs_control_EF$padj < alpha & !is.na(res_warm_vs_control_EF$padj)) & res_warm_vs_control_EF$log2FoldChange > l2fc), ]
signif_res_warm_vs_control_up_EF <- sort(as.character(rownames(signif_res_warm_vs_control_up_EF)), decreasing = TRUE)
#downregulated
signif_res_warm_vs_control_down_EF <- res_warm_vs_control_EF[ which((res_warm_vs_control_EF$padj < alpha & !is.na(res_warm_vs_control_EF$padj)) & res_warm_vs_control_EF$log2FoldChange < -l2fc), ]
signif_res_warm_vs_control_down_EF <- sort(as.character(rownames(signif_res_warm_vs_control_down_EF)), decreasing = TRUE)

#GO enrichment analysis for upregulated genes
ego_warm_vs_control_up_EF <- enrichGO(gene = signif_res_warm_vs_control_up_EF,
                                      universe = all_genes_EF_warm_vs_control,
                                      OrgDb = "org.At.tair.db",
                                      keyType = "TAIR",
                                      ont = "BP", #biological process
                                      pAdjustMethod = "BH",
                                      readable = TRUE)

dotplot(ego_warm_vs_control_up_EF, showCategory = 20)
upsetplot(ego_warm_vs_control_up_EF, n = 10)

png("results/EF/enrich_warm_vs_control_up_EF_zhang.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_up_EF, showCategory = 15, title = "EF warm_vs_control up LFC 1.0")
dev.off()

#more detail
write.csv2(as.data.frame(ego_warm_vs_control_up_EF), file="results/EF/warm_vs_control_up_EF_zhang.csv")

#GO enrichment analysis downregulated genes
ego_warm_vs_control_down_EF <- enrichGO(gene = signif_res_warm_vs_control_down_EF,
                                        universe = all_genes_EF_warm_vs_control,
                                        OrgDb = "org.At.tair.db",
                                        keyType = "TAIR",
                                        ont = "BP", #biological process
                                        pAdjustMethod = "BH",
                                        readable = TRUE)

dotplot(ego_warm_vs_control_down_EF, showCategory = 20)
upsetplot(ego_warm_vs_control_down_EF, n = 10)

png("results/EF/enrich_warm_vs_control_down_EF_zhang.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_down_EF, showCategory = 15, title = "EF warm_vs_control down LFC 1.0")
dev.off()

#more detail
write.csv2(as.data.frame(ego_warm_vs_control_down_EF), file="results/EF/warm_vs_control_down_EF_zhang.csv")


######################## ANALYSING GROUP LATE FLOWERS ######################################
#create subset
samples3 <- subset(samples, type == "LF")
head(samples3)

#paths to quant.sf
quantFiles3  <- file.path("star_salmon",samples3$sample, "quant.sf")
names(quantFiles3) <- samples3$sample

#import files with tximport
txi3 <- tximport(quantFiles3, type = "salmon", tx2gene = tx2gene)

#create DESeqDataSet
dds3 <- DESeqDataSetFromTximport(txi = txi3,
                                 colData = samples3,
                                 design = ~ group) 
dim(dds3)

############## Visualization #################
# pre-filter
smallestGroupSize <- 3
keep <- rowSums(counts(dds3) >= 10) >= smallestGroupSize
dds3 <- dds3[keep,]
dim(dds3)
#give factor levels
dds3$group <- relevel(dds3$group, ref = "control")
dds3$group

#transform data for visualization
vsd <- vst(dds3, blind=TRUE) #variance stabilizing transformations
rld <- rlog(dds3, blind=TRUE)
head(assay(vsd), 3)
ntd <- normTransform(dds3)

#comparing the transformations
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld)) 
meanSdPlot(assay(vsd))# this one gave the best flat curve

png("results/LF/meanSdPlot_vsd_LF_zhang.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(vsd)) 
dev.off()

#create Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group)
colnames(sampleDistMatrix) <- NULL

png("results/LF/heatmap_LF_zhang.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()


#create PCA plot
pcaData <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("results/LF/pca_LF_zhang.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

############## Differential expression analysis #############
##############      with no-normalized data     #############
dds3 <- DESeq(dds3)
res3 <- results(dds3)
res3
write.table(res3, file="results/LF/res_table_LF_zhang.csv", sep="\t", quote=F, col.names=NA)


#exporting normalized counts
normalized_counts <- counts(dds3, normalized=TRUE)
write.table(normalized_counts, file="results/LF/normalized_counts_LF_zhang.txt", sep="\t", quote=F, col.names=NA)


#plot Dispersion Estimate 
png("results/LF/dispersion_estimate_LF_zhang.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds3)
dev.off()

#exploring name for group
resultsNames(dds3)

res_warm_vs_control_LF <- results(dds3, name="group_warm_vs_control")

#create MA plot to see different expressed genes
png("results/LF/MA_plots_LF_zhang.png", width = 18, height = 16, units = "cm", res=300)
plotMA(res_warm_vs_control_LF, ylim=c(-2.5,2.5), main= "warm_vs_control")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
dev.off()

#set values log2FoldChange and p-value
l2fc = 1.0
alpha = 0.05

#############################
# compare warm with control #
#############################

all_genes_LF_warm_vs_control <- sort(as.character(rownames(res_warm_vs_control_LF)), decreasing = TRUE)

#extract significant results
#upregulated
signif_res_warm_vs_control_up_LF <- res_warm_vs_control_LF[ which((res_warm_vs_control_LF$padj < alpha & !is.na(res_warm_vs_control_LF$padj)) & res_warm_vs_control_LF$log2FoldChange > l2fc), ]
signif_res_warm_vs_control_up_LF <- sort(as.character(rownames(signif_res_warm_vs_control_up_LF)), decreasing = TRUE)
#downregulated
signif_res_warm_vs_control_down_LF <- res_warm_vs_control_LF[ which((res_warm_vs_control_LF$padj < alpha & !is.na(res_warm_vs_control_LF$padj)) & res_warm_vs_control_LF$log2FoldChange < -l2fc), ]
signif_res_warm_vs_control_down_LF <- sort(as.character(rownames(signif_res_warm_vs_control_down_LF)), decreasing = TRUE)

#GO enrichment analysis for upregulated genes
ego_warm_vs_control_up_LF <- enrichGO(gene = signif_res_warm_vs_control_up_LF,
                                      universe = all_genes_LF_warm_vs_control,
                                      OrgDb = "org.At.tair.db",
                                      keyType = "TAIR",
                                      ont = "BP", #biological process
                                      pAdjustMethod = "BH",
                                      readable = TRUE)

dotplot(ego_warm_vs_control_up_LF, showCategory = 20)
upsetplot(ego_warm_vs_control_up_LF, n = 10)

png("results/LF/enrich_warm_vs_control_up_LF_zhang.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_up_LF, showCategory = 15, title = "LF warm_vs_control up LFC 1.0")
dev.off()

#more detail
write.csv2(as.data.frame(ego_warm_vs_control_up_LF), file="results/LF/warm_vs_control_up_LF_zhang.csv")

#GO enrichment analysis for downregulated genes
ego_warm_vs_control_down_LF <- enrichGO(gene = signif_res_warm_vs_control_down_LF,
                                        universe = all_genes_LF_warm_vs_control,
                                        OrgDb = "org.At.tair.db",
                                        keyType = "TAIR",
                                        ont = "BP", #biological process
                                        pAdjustMethod = "BH",
                                        readable = TRUE)

dotplot(ego_warm_vs_control_down_LF, showCategory = 20)
upsetplot(ego_warm_vs_control_down_LF, n = 10)

png("results/LF/enrich_warm_vs_control_down_LF_zhang.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_down_LF, showCategory = 15, title = "LF warm_vs_control down LFC 1.0")
dev.off()

#more detail
write.csv2(as.data.frame(ego_warm_vs_control_down_LF), file="results/LF/warm_vs_control_down_LF_zhang.csv")



#########################################################
# compare results with genes predicted by the ML models #
#########################################################

#put all different expressed genes together
all_up_genes <- union(signif_res_warm_vs_control_up_RL, union(signif_res_warm_vs_control_up_EF, signif_res_warm_vs_control_up_LF))
all_down_genes <- union(signif_res_warm_vs_control_down_RL, union(signif_res_warm_vs_control_down_EF, signif_res_warm_vs_control_down_LF))
diff_genes <- union(all_up_genes, all_down_genes)
write.table(diff_genes, file="results/all_diff_reg_genes_zhang.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

############ab hier updaten####################################
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
write.table(rfe_general, file="results/rfe_genes_zhang.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)


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
write.table(lasso_general, file="results/lasso_genes_zhang.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

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
write.table(full_general, file="results/full_genes_zhang.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)



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
write.table(over_general, file="results/overlap_genes_zhang.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)



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

png("results/models_diff_regulated_zhang.png", width = 18, height = 16, units = "cm", res = 300)
h <- barplot(values, xlab = "Model", ylab = "Found Genes in %",
             main=paste("Comparison of the ML-models (p-value=",alpha," and l2fc=",l2fc,")\nNumber of genes found in Zhang et al., 2017: ", length(diff_genes), sep=""), cex.main=1, col="lightblue", names.arg=names)
text(h, c(4), labels = paste(round(values, digits=1), "%"), pos=3)
dev.off()
