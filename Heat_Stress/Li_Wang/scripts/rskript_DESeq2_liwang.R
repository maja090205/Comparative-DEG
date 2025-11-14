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

setwd("C:/Users/souku/job/eneza/heat_treatment/li_wang/")

#read information from samplesheet
samples <- read.csv("skripts/samplesheet_DESeq2_liwang.csv", header=TRUE)

samples$group <- as.character(samples$group)
samples$repetition <- as.character(samples$repetition)

#paths to quant.sf
quantFiles <- file.path("star_salmon", samples$sample, "quant.sf")
names(quantFiles) <- samples$sample

#read tx2gene file 
tx2gene <- read_tsv("star_salmon/tx2gene.tsv")

#import files with tximport
txi <- tximport(quantFiles, type = "salmon", tx2gene = tx2gene)

#create DESeqDataSet
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ group) 
dim(dds)

############## Visualization #################
# pre-filter 
smallestGroupSize <- 3 #three biological replicates
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds)
#give factor levels
dds$group <- relevel(dds$group, ref = "control")
dds$group

#transform data for visualization
vsd <- vst(dds, blind=TRUE) #variance stabilizing transformations
rld <- rlog(dds, blind=TRUE)
head(assay(vsd), 3)
ntd <- normTransform(dds)

#comparing the transformations
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld))# this one gave the best flat curve 
meanSdPlot(assay(vsd)) 

png("results/meanSdPlot_rld_liwang.png", width = 14, height = 12, units = "cm", res=300)
meanSdPlot(assay(rld)) 
dev.off()

#create Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$group)
colnames(sampleDistMatrix) <- NULL

png("results/heatmap_liwang.png", width = 14, height = 12, units = "cm", res=300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
dev.off()


#create PCA plot
pcaData <- plotPCA(rld, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("results/pca_liwang.png", width = 18, height = 16, units = "cm", res=300)
ggplot(pcaData, aes(PC1, PC2, color=group)) +
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
write.table(res, file="results/res_table_liwang.csv", sep="\t", quote=F, col.names=NA)


#exporting normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="results/normalized_counts_liwang.txt", sep="\t", quote=F, col.names=NA)


#plot Dispersion Estimate 
png("results/dispersion_estimate_liwang.png", width = 18, height = 16, units = "cm", res=300)
plotDispEsts(dds)
dev.off()

#exploring names of the groups
resultsNames(dds)

res_heat_vs_control <- results(dds, name="group_heat_vs_control")
res_warm_vs_control <- results(dds, name="group_warm_vs_control")

#create MA plot to see differential expressed genes
png("results/MA_plots_liwang.png", width = 18, height = 16, units = "cm", res=300)
par(mfrow = c(2, 2), mar=c(2,2,2,2))
plotMA(res_heat_vs_control, ylim=c(-2.5,2.5), main= "heat_vs_control")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
plotMA(res_warm_vs_control, ylim=c(-2.5,2.5), main= "warm_vs_control")
abline(h=c(-1,1),col="dodgerblue",lwd=2)
dev.off()

#set values for log2FoldChange and p-value
l2fc = 1.0
alpha = 0.05

#############################
# compare heat with control #
#############################

all_genes_heat_vs_control<- sort(as.character(rownames(res_heat_vs_control)), decreasing = TRUE)

#extract significant results
#upregulated
signif_res_heat_vs_control_up <- res_heat_vs_control[ which((res_heat_vs_control$padj < alpha & !is.na(res_heat_vs_control$padj)) & res_heat_vs_control$log2FoldChange > l2fc), ]
signif_res_heat_vs_control_up <- sort(as.character(rownames(signif_res_heat_vs_control_up)), decreasing = TRUE)
#downregulated
signif_res_heat_vs_control_down <- res_heat_vs_control[ which((res_heat_vs_control$padj < alpha & !is.na(res_heat_vs_control$padj)) & res_heat_vs_control$log2FoldChange < -l2fc), ]
signif_res_heat_vs_control_down <- sort(as.character(rownames(signif_res_heat_vs_control_down)), decreasing = TRUE)

#GO enrichment analysis for upregulated genes
ego_heat_vs_control_up <- enrichGO(gene = signif_res_heat_vs_control_up,
                          universe = all_genes_heat_vs_control,
                          OrgDb = "org.At.tair.db",
                          keyType = "TAIR",
                          ont = "BP", #biological process
                          pAdjustMethod = "BH",
                          readable = TRUE)

dotplot(ego_heat_vs_control_up, showCategory = 20)
upsetplot(ego_heat_vs_control_up, n = 10)

png("results/enrich_heat_vs_control_up_liwang.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_heat_vs_control_up, showCategory = 15, title = "group heat_vs_control up LFC 1.0")
dev.off()

#for more detail
write.csv2(as.data.frame(ego_heat_vs_control_up), file="results/heat_vs_control_up_liwang.csv")

#GO enrichment analysis downregulated genes
ego_heat_vs_control_down <- enrichGO(gene = signif_res_heat_vs_control_down,
                            universe = all_genes_heat_vs_control,
                            OrgDb = "org.At.tair.db",
                            keyType = "TAIR",
                            ont = "BP", #biological process
                            pAdjustMethod = "BH",
                            readable = TRUE)

dotplot(ego_heat_vs_control_down, showCategory = 20)
upsetplot(ego_heat_vs_control_down, n = 10)

png("results/enrich_heat_vs_control_down_liwang.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_heat_vs_control_down, showCategory = 15, title = "group heat_vs_control down LFC 1.0")
dev.off()

#for more detail
write.csv2(as.data.frame(ego_heat_vs_control_down), file="results/heat_vs_control_down_liwang.csv")


#############################
# compare warm with control #
#############################

all_genes_warm_vs_control <- sort(as.character(rownames(res_warm_vs_control)), decreasing = TRUE)

#extract significant results
#upregulated
signif_res_warm_vs_control_up <- res_warm_vs_control[ which((res_warm_vs_control$padj < alpha & !is.na(res_warm_vs_control$padj)) & res_warm_vs_control$log2FoldChange > l2fc), ]
signif_res_warm_vs_control_up <- sort(as.character(rownames(signif_res_warm_vs_control_up)), decreasing = TRUE)
#downregulated
signif_res_warm_vs_control_down <- res_warm_vs_control[ which((res_warm_vs_control$padj < alpha & !is.na(res_warm_vs_control$padj)) & res_warm_vs_control$log2FoldChange < -l2fc), ]
signif_res_warm_vs_control_down <- sort(as.character(rownames(signif_res_warm_vs_control_down)), decreasing = TRUE)

#GO enrichment analysis for upregulated genes
ego_warm_vs_control_up <- enrichGO(gene = signif_res_warm_vs_control_up,
                          universe = all_genes_warm_vs_control,
                          OrgDb = "org.At.tair.db",
                          keyType = "TAIR",
                          ont = "BP", #biological process
                          pAdjustMethod = "BH",
                          readable = TRUE)

dotplot(ego_warm_vs_control_up, showCategory = 20)
upsetplot(ego_warm_vs_control_up, n = 10)

png("results/enrich_warm_vs_control_up_liwang.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_up, showCategory = 15, title = "warm_vs_control up LFC 1.0")
dev.off()

#for more detail
write.csv2(as.data.frame(ego_warm_vs_control_up), file="results/warm_vs_control_up_liwang.csv")

#GO enrichment analysis for downregulated genes
ego_warm_vs_control_down <- enrichGO(gene = signif_res_warm_vs_control_down,
                            universe = all_genes_warm_vs_control,
                            OrgDb = "org.At.tair.db",
                            keyType = "TAIR",
                            ont = "BP", #biological process
                            pAdjustMethod = "BH",
                            readable = TRUE)

dotplot(ego_warm_vs_control_down, showCategory = 20)
upsetplot(ego_warm_vs_control_down, n = 10)

png("results/enrich_warm_vs_control_down_liwang.png", width = 18, height = 16, units = "cm", res=300)
barplot(ego_warm_vs_control_down, showCategory = 15, title = "warm_vs_control down LFC 1.0")
dev.off()

#more detail
write.csv2(as.data.frame(ego_warm_vs_control_down), file="results/warm_vs_control_down_liwang.csv")



#########################################################
# compare results with genes predicted by the ML models #
#########################################################

all_up_genes <- union(signif_res_heat_vs_control_up, signif_res_warm_vs_control_up)
all_down_genes <- union(signif_res_heat_vs_control_down, signif_res_warm_vs_control_down)
diff_genes <- union(all_up_genes, all_down_genes)
write.table(diff_genes, file="results/diff_reg_genes_liwang.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

#genes predicted by the different models
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
write.table(rfe_general, file="results/rfe_genes_liwang.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)


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
write.table(lasso_general, file="results/lasso_genes_liwang.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

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
write.table(full_general, file="results/full_genes_liwang.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)



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
write.table(over_general, file="results/overlap_genes_liwang.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)




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

png("results/models_diff_regulated_liwang.png", width = 18, height = 16, units = "cm", res = 300)
h <- barplot(values, xlab = "Model", ylab = "Found Genes in %",
             main=paste("Comparison of the ML-models (p-value=",alpha," and l2fc=",l2fc,")\nNumber of genes found in Li Wang et al., 2020: ", length(diff_genes), sep=""), cex.main=1, col="lightblue", names.arg=names)
text(h, c(4), labels = paste(round(values, digits=1), "%"), pos=3)
dev.off()
