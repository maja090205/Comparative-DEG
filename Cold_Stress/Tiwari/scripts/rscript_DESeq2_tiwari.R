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
samples <- read.csv("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/skripts/samplesheet_DESeq2.csv", header=TRUE)

samples$timepoint <- as.character(samples$timepoint)#samples were extracted after 3h, 6h and 48h
samples$replicate <- as.character(samples$replicate)
samples$treatment <- as.character(samples$treatment)#treatment defines wether the sample is from the control or cold treatment


#read tx2gene file
tx2gene <- read_tsv("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/star_salmon/tx2gene.tsv")

######################## ANALYSING TIMEPOINT 3h ######################################
#create subset
samples1 <- subset(samples, timepoint == "3h")
samples1

#paths to quant.sf
quantFiles1 <- file.path("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/star_salmon", samples1$sample, "quant.sf")
names(quantFiles1) <- samples1$sample

#import files with tximport
txi1 <- tximport(quantFiles1, type = "salmon", tx2gene = tx2gene)

#create a DESeqDataSet
dds1 <- DESeqDataSetFromTximport(txi1,
                                colData = samples1,
                                design = ~ treatment) 
dim(dds1)

############## Visualization #################
# pre-filter
smallestGroupSize <- 3 #biological triplicates
keep <- rowSums(counts(dds1) >= 10) >= smallestGroupSize
dds1 <- dds1[keep,]
dim(dds1)

#give factor levels
dds1$treatment <- relevel(dds1$treatment, ref = "control")
dds1$treatment

#transform data for visualization
vsd <- vst(dds1, blind=TRUE) #variance stabilizing transformations
rld <- rlog(dds1, blind=TRUE) #regularized logarithm
head(assay(vsd), 3)
ntd <- normTransform(dds1)

#comparing the transformations
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld))
meanSdPlot(assay(vsd)) # this one gave the best flat curve  

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/3h/meanSdPlot_vsd_3h_tiwari.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(vsd)) 
dev.off()

#create Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment)
colnames(sampleDistMatrix) <- NULL

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/3h/heatmap_3h_tiwari.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()


#create PCA plot
pcaData <- plotPCA(vsd, intgroup=c("treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/3h/pca_3h_tiwari.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData, aes(PC1, PC2, color=treatment)) +
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
write.table(res1, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/3h/res_table_3h_tiwari.csv", sep="\t", quote=F, col.names=NA)

#exporting normalized counts
normalized_counts <- counts(dds1, normalized=TRUE)
write.table(normalized_counts, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/3h/normalized_counts_3h_tiwari.txt", sep="\t", quote=F, col.names=NA)


#plot Dispersion Estimate 
png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/3h/dispersion_estimate_3h_tiwari.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds1)
dev.off()

#exploring name of group 
resultsNames(dds1)

res_cold_vs_control_3h <- results(dds1, name="treatment_cold_vs_control")

#create MA plot to see differential expressed genes
png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/3h/MA_plots_3h_tiwari.png", width = 18, height = 16, units = "cm", res=300)
plotMA(res_cold_vs_control_3h, ylim=c(-2.5,2.5), main= "cold_vs_control")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
dev.off()

#set values of log2FoldChange and p-value
l2fc = 1.0
alpha = 0.05

########################################
# compare cold treatment with control  #
########################################

all_genes_cold_vs_control_3h<- sort(as.character(rownames(res_cold_vs_control_3h)), decreasing = TRUE)

#extract significant results
#upregulated
signif_res_cold_vs_control_up_3h <- res_cold_vs_control_3h[ which((res_cold_vs_control_3h$padj < alpha & !is.na(res_cold_vs_control_3h$padj)) & res_cold_vs_control_3h$log2FoldChange > l2fc), ]
signif_res_cold_vs_control_up_3h <- sort(as.character(rownames(signif_res_cold_vs_control_up_3h)), decreasing = TRUE)
#downregulated
signif_res_cold_vs_control_down_3h <- res_cold_vs_control_3h[ which((res_cold_vs_control_3h$padj < alpha & !is.na(res_cold_vs_control_3h$padj)) & res_cold_vs_control_3h$log2FoldChange < -l2fc), ]
signif_res_cold_vs_control_down_3h <- sort(as.character(rownames(signif_res_cold_vs_control_down_3h)), decreasing = TRUE)

#GO enrichment analysis for upregulated genes
ego_cold_vs_control_up_3h <- enrichGO(gene = signif_res_cold_vs_control_up_3h,
                          universe = all_genes_cold_vs_control_3h,
                          OrgDb = "org.At.tair.db",
                          keyType = "TAIR",
                          ont = "BP", #biological process
                          pAdjustMethod = "BH",
                          readable = TRUE)

dotplot(ego_cold_vs_control_up_3h, showCategory = 20)
upsetplot(ego_cold_vs_control_up_3h, n = 10)

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/3h/enrich_cold_vs_control_up_tiwari.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_cold_vs_control_up_3h, showCategory = 15, title = "cold_vs_control 3h up LFC 1.0")
dev.off()

#more details
write.csv2(as.data.frame(ego_cold_vs_control_up_3h), file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/3h/cold_vs_control_up_tiwari.csv")

#GO enrichment analysis for downregulated genes
ego_cold_vs_control_down_3h <- enrichGO(gene = signif_res_cold_vs_control_down_3h,
                            universe = all_genes_cold_vs_control_3h,
                            OrgDb = "org.At.tair.db",
                            keyType = "TAIR",
                            ont = "BP", #biological process
                            pAdjustMethod = "BH",
                            readable = TRUE)

dotplot(ego_cold_vs_control_down_3h, showCategory = 20)
upsetplot(ego_cold_vs_control_down_3h, n = 10)

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/3h/enrich_cold_vs_control_down_tiwari.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_cold_vs_control_down_3h, showCategory = 15, title = "cold_vs_control 3h down LFC 1.0")
dev.off()

#more details
write.csv2(as.data.frame(ego_cold_vs_control_down_3h), file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/3h/cold_vs_control_down_tiwari.csv")


######################## ANALYSING TIMEPOINT 6h ######################################
#create subset
samples2 <- subset(samples, timepoint == "6h")
samples2

#paths to quant.sf
quantFiles2 <- file.path("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/star_salmon", samples2$sample, "quant.sf")
names(quantFiles2) <- samples2$sample

#import files with tximport
txi2 <- tximport(quantFiles2, type = "salmon", tx2gene = tx2gene)

#create a DESeqDataSet
dds2 <- DESeqDataSetFromTximport(txi2,
                                 colData = samples2,
                                 design = ~ treatment) 
dim(dds2)

############## Visualization #################
# pre-filter
smallestGroupSize <- 3 #biological triplicates
keep <- rowSums(counts(dds2) >= 10) >= smallestGroupSize
dds2 <- dds2[keep,]
dim(dds2)

#give factor levels
dds2$treatment <- relevel(dds2$treatment, ref = "control")
dds2$treatment

#transform data for visualization
vsd <- vst(dds2, blind=TRUE) #variance stabilizing transformations
rld <- rlog(dds2, blind=TRUE) #regularized logarithm
head(assay(vsd), 3)
ntd <- normTransform(dds2)

#comparing the transformations
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld))# this one gave the best flat curve
meanSdPlot(assay(vsd))  

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/6h/meanSdPlot_rld_6h_tiwari.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(rld)) 
dev.off()

#create Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$treatment)
colnames(sampleDistMatrix) <- NULL

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/6h/heatmap_6h_tiwari.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()


#create PCA plot
pcaData <- plotPCA(rld, intgroup=c("treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/6h/pca_6h_tiwari.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData, aes(PC1, PC2, color=treatment)) +
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
write.table(res2, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/6h/res_table_6h_tiwari.csv", sep="\t", quote=F, col.names=NA)

#exporting normalized counts
normalized_counts <- counts(dds2, normalized=TRUE)
write.table(normalized_counts, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/6h/normalized_counts_6h_tiwari.txt", sep="\t", quote=F, col.names=NA)


#plot Dispersion Estimate 
png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/6h/dispersion_estimate_6h_tiwari.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds2)
dev.off()

#exploring name of group 
resultsNames(dds2)

res_cold_vs_control_6h <- results(dds2, name="treatment_cold_vs_control")

#create MA plot to see differential expressed genes
png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/6h/MA_plots_6h_tiwari.png", width = 18, height = 16, units = "cm", res=300)
plotMA(res_cold_vs_control_6h, ylim=c(-2.5,2.5), main= "cold_vs_control")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
dev.off()

#set values of log2FoldChange and p-value
l2fc = 1.0
alpha = 0.05

########################################
# compare cold treatment with control  #
########################################

all_genes_cold_vs_control_6h<- sort(as.character(rownames(res_cold_vs_control_6h)), decreasing = TRUE)

#extract significant results
#upregulated
signif_res_cold_vs_control_up_6h <- res_cold_vs_control_6h[ which((res_cold_vs_control_6h$padj < alpha & !is.na(res_cold_vs_control_6h$padj)) & res_cold_vs_control_6h$log2FoldChange > l2fc), ]
signif_res_cold_vs_control_up_6h <- sort(as.character(rownames(signif_res_cold_vs_control_up_6h)), decreasing = TRUE)
#downregulated
signif_res_cold_vs_control_down_6h <- res_cold_vs_control_6h[ which((res_cold_vs_control_6h$padj < alpha & !is.na(res_cold_vs_control_6h$padj)) & res_cold_vs_control_6h$log2FoldChange < -l2fc), ]
signif_res_cold_vs_control_down_6h <- sort(as.character(rownames(signif_res_cold_vs_control_down_6h)), decreasing = TRUE)

#GO enrichment analysis for upregulated genes
ego_cold_vs_control_up_6h <- enrichGO(gene = signif_res_cold_vs_control_up_6h,
                                   universe = all_genes_cold_vs_control_6h,
                                   OrgDb = "org.At.tair.db",
                                   keyType = "TAIR",
                                   ont = "BP", #biological process
                                   pAdjustMethod = "BH",
                                   readable = TRUE)

dotplot(ego_cold_vs_control_up_6h, showCategory = 20)
upsetplot(ego_cold_vs_control_up_6h, n = 10)

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/6h/enrich_cold_vs_control_up_tiwari.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_cold_vs_control_up_6h, showCategory = 15, title = "cold_vs_control 6h up LFC 1.0")
dev.off()

#more details
write.csv2(as.data.frame(ego_cold_vs_control_up_6h), file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/6h/cold_vs_control_up_tiwari.csv")

#GO enrichment analysis for downregulated genes
ego_cold_vs_control_down_6h <- enrichGO(gene = signif_res_cold_vs_control_down_6h,
                                     universe = all_genes_cold_vs_control_6h,
                                     OrgDb = "org.At.tair.db",
                                     keyType = "TAIR",
                                     ont = "BP", #biological process
                                     pAdjustMethod = "BH",
                                     readable = TRUE)

dotplot(ego_cold_vs_control_down_6h, showCategory = 20)
upsetplot(ego_cold_vs_control_down_6h, n = 10)

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/6h/enrich_cold_vs_control_down_tiwari.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_cold_vs_control_down_6h, showCategory = 15, title = "cold_vs_control 6h down LFC 1.0")
dev.off()

#more details
write.csv2(as.data.frame(ego_cold_vs_control_down_6h), file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/6h/cold_vs_control_down_tiwari.csv")


######################## ANALYSING TIMEPOINT 48h ######################################
#create subset
samples3 <- subset(samples, timepoint == "48h")
samples3

#paths to quant.sf
quantFiles3 <- file.path("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/star_salmon", samples3$sample, "quant.sf")
names(quantFiles3) <- samples3$sample

#import files with tximport
txi3 <- tximport(quantFiles3, type = "salmon", tx2gene = tx2gene)

#create a DESeqDataSet
dds3 <- DESeqDataSetFromTximport(txi3,
                                 colData = samples3,
                                 design = ~ treatment) 
dim(dds3)

############## Visualization #################
# pre-filter
smallestGroupSize <- 3 #biological triplicates
keep <- rowSums(counts(dds3) >= 10) >= smallestGroupSize
dds3 <- dds3[keep,]
dim(dds3)

#give factor levels
dds3$treatment <- relevel(dds3$treatment, ref = "control")
dds3$treatment

#transform data for visualization
vsd <- vst(dds3, blind=TRUE) #variance stabilizing transformations
rld <- rlog(dds3, blind=TRUE) #regularized logarithm
head(assay(vsd), 3)
ntd <- normTransform(dds3)

#comparing the transformations
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld))
meanSdPlot(assay(vsd)) # this one gave the best flat curve  

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/48h/meanSdPlot_vsd_48h_tiwari.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(vsd)) 
dev.off()

#create Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment)
colnames(sampleDistMatrix) <- NULL

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/48h/heatmap_48h_tiwari.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()


#create PCA plot
pcaData <- plotPCA(vsd, intgroup=c("treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/48h/pca_48h_tiwari.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData, aes(PC1, PC2, color=treatment)) +
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
write.table(res3, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/48h/res_table_48h_tiwari.csv", sep="\t", quote=F, col.names=NA)

#exporting normalized counts
normalized_counts <- counts(dds3, normalized=TRUE)
write.table(normalized_counts, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/48h/normalized_counts_48h_tiwari.txt", sep="\t", quote=F, col.names=NA)


#plot Dispersion Estimate 
png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/48h/dispersion_estimate_48h_tiwari.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds3)
dev.off()

#exploring name of group 
resultsNames(dds3)

res_cold_vs_control_48h <- results(dds3, name="treatment_cold_vs_control")

#create MA plot to see differential expressed genes
png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/48h/MA_plots_48h_tiwari.png", width = 18, height = 16, units = "cm", res=300)
plotMA(res_cold_vs_control_48h, ylim=c(-2.5,2.5), main= "cold_vs_control")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
dev.off()

#set values of log2FoldChange and p-value
l2fc = 1.0
alpha = 0.05

########################################
# compare cold treatment with control  #
########################################

all_genes_cold_vs_control_48h<- sort(as.character(rownames(res_cold_vs_control_48h)), decreasing = TRUE)

#extract significant results
#upregulated
signif_res_cold_vs_control_up_48h <- res_cold_vs_control_48h[ which((res_cold_vs_control_48h$padj < alpha & !is.na(res_cold_vs_control_48h$padj)) & res_cold_vs_control_48h$log2FoldChange > l2fc), ]
signif_res_cold_vs_control_up_48h <- sort(as.character(rownames(signif_res_cold_vs_control_up_48h)), decreasing = TRUE)
#downregulated
signif_res_cold_vs_control_down_48h <- res_cold_vs_control_48h[ which((res_cold_vs_control_48h$padj < alpha & !is.na(res_cold_vs_control_48h$padj)) & res_cold_vs_control_48h$log2FoldChange < -l2fc), ]
signif_res_cold_vs_control_down_48h <- sort(as.character(rownames(signif_res_cold_vs_control_down_48h)), decreasing = TRUE)

#GO enrichment analysis for upregulated genes
ego_cold_vs_control_up_48h <- enrichGO(gene = signif_res_cold_vs_control_up_48h,
                                   universe = all_genes_cold_vs_control_48h,
                                   OrgDb = "org.At.tair.db",
                                   keyType = "TAIR",
                                   ont = "BP", #biological process
                                   pAdjustMethod = "BH",
                                   readable = TRUE)

dotplot(ego_cold_vs_control_up_48h, showCategory = 20)
upsetplot(ego_cold_vs_control_up_48h, n = 10)

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/48h/enrich_cold_vs_control_up_tiwari.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_cold_vs_control_up_48h, showCategory = 15, title = "cold_vs_control 48h up LFC 1.0")
dev.off()

#more details
write.csv2(as.data.frame(ego_cold_vs_control_up_48h), file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/48h/cold_vs_control_up_tiwari.csv")

#GO enrichment analysis for downregulated genes
ego_cold_vs_control_down_48h <- enrichGO(gene = signif_res_cold_vs_control_down_48h,
                                     universe = all_genes_cold_vs_control_48h,
                                     OrgDb = "org.At.tair.db",
                                     keyType = "TAIR",
                                     ont = "BP", #biological process
                                     pAdjustMethod = "BH",
                                     readable = TRUE)

dotplot(ego_cold_vs_control_down_48h, showCategory = 20)
upsetplot(ego_cold_vs_control_down_48h, n = 10)

png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/48h/enrich_cold_vs_control_down_tiwari.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_cold_vs_control_down_48h, showCategory = 15, title = "cold_vs_control 48h down LFC 1.0")
dev.off()

#more details
write.csv2(as.data.frame(ego_cold_vs_control_down_48h), file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/48h/cold_vs_control_down_tiwari.csv")



#########################################################
# compare results with genes predicted by the ML models #
#########################################################

#extract all differential expressed genes found in the publication 
all_up_genes <- union(signif_res_cold_vs_control_up_3h, union(signif_res_cold_vs_control_up_6h, signif_res_cold_vs_control_up_48h))
all_down_genes <- union(signif_res_cold_vs_control_down_3h, union(signif_res_cold_vs_control_down_6h, signif_res_cold_vs_control_down_48h))
diff_genes <- union(all_up_genes, all_down_genes)
write.table(diff_genes, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/all_diff_reg_genes_tiwari.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)


##############
# RFE model  #
##############
#genes predicted by the RFE model
data1 <- read.csv2("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/rfe model - Important.csv", header = TRUE, sep = ",")

#extract the genes
rfe <- data1$Genes
length(rfe)

#intersection of the upregulated genes
rfe_up <- intersect(rfe, all_up_genes)
length(rfe_up)

print(rfe_up)
write.table(rfe_up, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/rfe_up_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#intersection of the downregulated genes
rfe_down <- intersect(rfe, all_down_genes)
length(rfe_down)

print(rfe_down)
write.table(rfe_down, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/rfe_down_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#general different regulated
rfe_general <- union(rfe_up, rfe_down)
length(rfe_general)

print(rfe_general)
write.table(rfe_general, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/rfe_general_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


##############
# XGB model  #
##############
#genes predicted by the XGB model
data2 <- read.csv2("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/xgb full - Important.csv", header = TRUE, sep = ",")

xgb <- data2$Genes
length(xgb)

#upregulated
xgb_up <- intersect(xgb, all_up_genes)
length(xgb_up)

print(xgb_up)
write.table(xgb_up, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/xgb_up_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#downregulated
xgb_down <- intersect(xgb, all_down_genes)
length(xgb_down)

print(xgb_down)
write.table(xgb_down, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/xgb_down_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#general different regulated
xgb_general <- union(xgb_up, xgb_down)
length(xgb_general)

print(xgb_general)
write.table(xgb_general, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/xgb_general_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


################
# Lasso model  #
################
#genes predicted by the Lasso model
data3 <- read.csv2("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/lasso model - Important.csv", header = TRUE, sep = ",")

lasso <- data3$Genes
length(lasso)

#upregulated
lasso_up <- intersect(lasso, all_up_genes)
length(lasso_up)

print(lasso_up)
write.table(lasso_up, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/lasso_up_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#downregulated
lasso_down <- intersect(lasso, all_down_genes)
length(lasso_down)

print(lasso_down)
write.table(lasso_down, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/lasso_down_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#general different regulated
lasso_general <- union(lasso_up, lasso_down)
length(lasso_general)

print(lasso_general)
write.table(lasso_general, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/lasso_general_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


#####################
# Overlap Lasso RFE #
#####################
#genes predicted by the Lasso and RFE model
data4 <- read.csv2("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/overlaps - Lasso RFE models.csv", header = TRUE, sep = ",")

over <- data4$Genes
length(over)

#upregulated
over_up <- intersect(over, all_up_genes)
length(over_up)

print(over_up)
write.table(over_up, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/overlap_up_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#downregulated
over_down <- intersect(over, all_down_genes)
length(over_down)

print(over_down)
write.table(over_down, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/overlap_down_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#general different regulated
over_general <- union(over_up, over_down)
length(over_general)

print(over_general)
write.table(over_general, file="C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/overlap_general_genes_tiwari.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


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


###########
# Plotten #
###########
#only general different regulated
names <- c("RFE", "XGB", "Lasso", "Overlap")
values <- c(percentage_rfe_general, percentage_xgb_general, percentage_lasso_general, percentage_over_general)


png("C:/Users/souku/job/eneza/cold_treatment/arabid_tiwari/results/all_models_diff_regulated_tiwari.png", width = 18, height = 16, units = "cm", res = 300)
h <- barplot(values, xlab = "Model", ylab = "Found Genes in %",
             main=paste("Comparison of the ML-models (p-value=",alpha," and l2fc=",l2fc,")\nNumber of genes found in Tiwari et al., 2020: ", length(diff_genes), sep=""), cex.main=1, col="lightblue", names.arg=names)
text(h, c(4), labels = paste(round(values, digits=1), "%"), pos=3)
dev.off()

