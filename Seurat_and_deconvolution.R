############################################################
#                                                          #
# Master's thesis: Single-cell RNA-Seq deconvolution       #
# of glioblastoma bulk transcriptomics                     #
# Seurat analysis and deconvolution                        #
#                                                          #
# Sofia Randelin                                           #
# 8.7.2019                                                 #
############################################################

library(dplyr)
library(Seurat)

# Single-cell aligned counts
singlecell_file <- "matrixcounts.txt"
count_mat <- read.table(singlecell_file, header=TRUE)

# Seurat: CreateSeuratObject to start. 
# Count matrix has UMIs as column names and ensembl gene IDs as row names
# Save to Seurat Object 'scrna'
scrna <- CreateSeuratObject(counts = count_mat, 
                            project = "reference", 
                            min.cells = 3, min.features = 200) 

# Mitochodrial percentage
scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
scrna <- subset(scrna, subset = percent.mt < 50)

# VlnPlot visualization of QC metrics
VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident")

# Visualize feature relationships
FeatureScatter(scrna, feature1 = "nFeature_RNA", 
               feature2 = "percent.mt")

# NORMALIZATION
scrna <- NormalizeData(scrna, 
                       normalization.method = "LogNormalize", 
                       scale.factor = 10000)

# Identification of 5000 highly variable features (feature selection)
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 5000)

#Scaling data
all.genes <- rownames(scrna)
scrna <- ScaleData(scrna, features = all.genes)

# Linear dimensional reduction with 5000 top variable genes
scrna <- RunPCA(scrna, 
                features = VariableFeatures(object = scrna)) 

# Examine and visualize PCA results a few different ways
print(scrna[["pca"]], dims = 1:5, nfeatures = 30)





#########################################################################
# Cell cycle scores where computed in one run for testing purposes
# and this part was excluded in the final analysis

# Cell cycle scores
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scrna <- CellCycleScoring(scrna, s.features = s.genes, 
                          g2m.features = g2m.genes, 
                          set.ident = TRUE)
scrna <- RunPCA(scrna, features = c(s.genes, g2m.genes))
DimPlot(scrna, group.by = "Phase")
scrna <- ScaleData(scrna, vars.to.regress = c("S.Score", "G2M.Score"), 
                   features = rownames(scrna))
scrna <- RunPCA(scrna, features = VariableFeatures(scrna), 
                nfeatures.print = 10)
DimPlot(scrna)

#########################################################################



# CLUSTERING
scrna <- FindNeighbors(scrna, dims = 1:10)
scrna <- FindClusters(scrna, resolution = 0.5)

table(scrna@active.ident)


# UMAP for visualizing clusters
scrna <- RunUMAP(scrna, dims = 1:10)
DimPlot(scrna, reduction = "umap", label = TRUE, label.size = 5)


#clusters 0, 1, 4, 8 (CD45 cells) --> New clustering
clusters <- scrna@meta.data$seurat_clusters
cell_idents <- scrna@active.ident
cd45 <- subset(scrna, idents = c("0", "1", "4", "8"), invert = FALSE)
cd45_counts <- cd45[["RNA"]]@counts

cd45_original_clusters <- cd45@active.ident
cd45 <- FindNeighbors(cd45, dims = 1:10)
cd45 <- FindClusters(cd45, resolution = 0.5)
cd45_new_clusters <-cd45@active.ident


# Update cd45 ident numbers
# so they do not mix with the scrna idents later
cd45 <- RenameIdents(cd45, `1` = "1a", `2` = "4", `3` = "0a", 
                     `0` = "0b", `4` = "1b", `5` = "8")
cd45_cell_idents <- cd45@active.ident

# Find matching cell names and create a vector of 
#new idents for clusters 0, 1, 4 and 8
new_clusters <- as.character(cell_idents)
indices <- match(names(cd45_cell_idents), names(cell_idents))
new_clusters[indices] <- as.character(cd45_cell_idents)
new_clusters <- as.character(new_clusters)

# Add new idents to scrna:
Idents(scrna) <- new_clusters

# scrna UMAP visualization
scrna <- RunUMAP(scrna, dims = 1:10)
DimPlot(scrna, reduction = "umap", label = TRUE, label.size = 10, pt.size = 3) + NoLegend()


# CD45 clusters only, visualization UMAP and tSNE
cd45 <- RunUMAP(cd45, dims = 1:10)
Dimplot(cd45, reduction = "umap", label = TRUE, label.size = 6)


# Test marker genes in article:
# "PTPRC", "EGFR", "MOG", "ETNPPL", "DCN", "GPR17", "STMN2"
# Plot marker genes in UMAP Featureplot and VlnPlot
FeaturePlot(scrna, 
            features = c("PTPRC", "EGFR", "MOG", "ETNPPL", 
                         "DCN", "GPR17", "STMN2"),
            pt.size = 1.5, combine = TRUE, label = TRUE, 
            label.size = 4)
VlnPlot(scrna, features = c("PTPRC", "EGFR", "MOG", "ETNPPL", 
                            "DCN", "GPR17", "STMN2"))



# Test different sets of marker genes
immunefeature <- c("CD4", "CD8A", "CD3E", "CD3G", "FOXP3", 
                   "ITGAM", "ITGAX", "CD68", "CD163", 
                   "CD14", "CX3CR1", "CSF1R", "MS4A1", 
                   "CEACAM8", "CD33", "CD24", "KIT")
leucocyte <- c("PTPRC")
neoplastic <- c("EGFR", "MDM2", "CDK4")
oligoden <- c("MOG")
opc <- c("GPR17", "PDGFRA")
astrocyte <- c("ETNPPL")
vascular <- c("DCN")
neuron <- c("STMN2")
bcell <- c("BLK","CD19","FCRL2","MS4A1","TNFRSF17",
           "TCL1A","SPIB", "PNOC")
cytotoxic <- c("PRF1","GZMA","GZMB","NKG7","KLRK1",
               "KLRB1","KLRD1","CTSW","GNLY")
DC <- c("ITGAM","ITGAX","CX3CR1","CSF1R","CCL13",
        "CD209","HSD11B1")
exhaustedcd8 <- c("LAG3","CD244","EOMES","PTGER4")
mac <- c("CD68","CD163","CD14","CD84","MS4A4A")
mast <- c("KIT","TPSB2","TPSAB1","CPA3","MS4A2","HDC")
neutro <- c("CEACAM3", "CEACAM8","CD33","CD24","FPR1",
            "SIGLEC5","CSF3R","FCAR","FCGR3B","S100A12")
nk <- c("IL21R")
tcell <- c("XCL1","XCL2")
th1 <- c("CD3E","CD3G", "CD4","CD6","SH2D1A","TRAT1")
treg <- c("FOXP3")
cd8 <- c("CD8A")

# Combine features for DotPlot plotting
all_features <- c(leucocyte,neoplastic,oligoden,
                  opc,astrocyte,vascular,neuron,bcell,
                  cytotoxic,DC,exhaustedcd8,mac,mast,
                  neutro,nk,tcell,th1,treg,cd8)
all_immune <- c(leucocyte, bcell,cytotoxic,DC,exhaustedcd8,
                mac,mast,neutro,nk,tcell,th1,treg,cd8)
all_other <- c(neoplastic,oligoden,opc,astrocyte,vascular,neuron)


# Load the "scales" package
require(scales)

PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

genegroups <- c("B cell", "Cytotoxic cell", "DC",
                "Exhausted CD8+ cell", "Macrophage",
                "Mast cell", "Neutrophil", "NK",
                "T cell", "Th1 cell", "Treg", "CD8+ T cell")
genegroups_all <- c("Leukocyte", "Neoplastic", "Oligodendrocyte",
                    "OPC", "Astrocyte", "Vascular cell", "Neuron",
                    "B cell", "Cytotoxic cell", "DC",
                    "Exhausted CD8+ cell", "Macrophage",
                    "Mast cell", "Neutrophil", "NK",
                    "T cell", "Th1 cell", "Treg", "CD8+ T cell")


inputgroups_all <- c(rep(genegroups_all[1], length(leucocyte)),
                     rep(genegroups_all[2], length(neoplastic)),
                     rep(genegroups_all[3], length(oligoden)),
                     rep(genegroups_all[4], length(opc)),
                     rep(genegroups_all[5], length(astrocyte)),
                     rep(genegroups_all[6], length(vascular)),
                     rep(genegroups_all[7], length(neuron)),
                     rep(genegroups_all[8], length(bcell)),
                     rep(genegroups_all[9], length(cytotoxic)), 
                     rep(genegroups_all[10], length(DC)),
                     rep(genegroups_all[11], length(exhaustedcd8)),
                     rep(genegroups_all[12], length(mac)),
                     rep(genegroups_all[13], length(mast)),
                     rep(genegroups_all[14], length(neutro)),
                     rep(genegroups_all[15], length(nk)),
                     rep(genegroups_all[16], length(tcell)),
                     rep(genegroups_all[17], length(th1)),
                     rep(genegroups_all[18], length(treg)),
                     rep(genegroups_all[19], length(cd8)))

inputgroups <- c(
  rep(genegroups[1], length(bcell)),
  rep(genegroups[2], length(cytotoxic)), 
  rep(genegroups[3], length(DC)),
  rep(genegroups[4], length(exhaustedcd8)),
  rep(genegroups[5], length(mac)),
  rep(genegroups[6], length(mast)),
  rep(genegroups[7], length(neutro)),
  rep(genegroups[8], length(nk)),
  rep(genegroups[9], length(tcell)),
  rep(genegroups[10], length(th1)),
  rep(genegroups[11], length(treg)),
  rep(genegroups[12], length(cd8)))

input_other <- c(
  rep(genegroups_all[2], length(neoplastic)),
  rep(genegroups_all[3], length(oligoden)),
  rep(genegroups_all[4], length(opc)),
  rep(genegroups_all[5], length(astrocyte)),
  rep(genegroups_all[6], length(vascular)),
  rep(genegroups_all[7], length(neuron)))

other <- subset(scrna, idents = c("0a","0b", "1a", "1b", "4", "8"), invert = TRUE)

# non-immune Dotplot
DotPlot(other, 
        features = all_other,
        #features = all_immune[1:length(all_immune)],
        gene.groups = input_other) +
  RotatedAxis()

# immune Dotplot
DotPlot(cd45, 
        features = all_immune[1:length(all_immune)],
        gene.groups = input_groups) +
  RotatedAxis()


scrna.markers <- FindAllMarkers(scrna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "t")
scrna.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- scrna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


# Optional: Find markers that separate immune clusters from each other
cluster0a.markers <- FindMarkers(object = scrna, only.pos = TRUE,
                                 ident.1 = '0a',
                                 ident.2 = c('0b','1a','1b','4'),
                                 min.pct = 0.25, logfc.threshold = 0.25)
cluster0b.markers <- FindMarkers(object = scrna, only.pos = TRUE,
                                 ident.1 = '0b',
                                 ident.2 = c('0a','1a','1b','4'),
                                 min.pct = 0.25, logfc.threshold = 0.25)
cluster1a.markers <- FindMarkers(object = scrna, only.pos = TRUE,
                                 ident.1 = '1a',
                                 ident.2 = c('0b','0a','1b','4'),
                                 min.pct = 0.25, logfc.threshold = 0.25)
cluster1b.markers <- FindMarkers(object = scrna, only.pos = TRUE,
                                 ident.1 = '1b',
                                 ident.2 = c('0b','1a','0a','4'),
                                 min.pct = 0.25, logfc.threshold = 0.25)
cluster4.markers <- FindMarkers(object = scrna, only.pos = TRUE,
                                ident.1 = '4',
                                ident.2 = c('0b','1a','1b','0a'),
                                min.pct = 0.25, logfc.threshold = 0.25)


##############################################################################################

# Deconvolution preparation of data
# Create deconvolution reference from scRNA-seq data


scrna <- RenameIdents(scrna, `9` = "Oligodendrocyte", `3` = "OPC", 
                      `0b` = "Neutrophil", `7` = "Neoplastic",
                      `1a` = "Dendritic cell", `4` = "Mast cell", 
                      `12` = "Vascular", `0a` = "Macrophage/Microglia", 
                      `1b` = "Dendritic cell", `11` = "Astrocyte", 
                      `13` = "Neuron", `12` = "Vascular", `2` = "Neoplastic", 
                      `5` = "Neoplastic",`6` = "Neoplastic", `10` = "Neoplastic", 
                      `8` = "G2M")

deconv_clusters <- subset(scrna, idents = c("OPC","Astrocyte","Vascular",
                                            "Neuron","Macrophage/Microglia",
                                            "Oligodendrocyte","Dendritic cell",
                                            "Mast cell","Neutrophil"), invert = FALSE)


# Update 'scrna' Object to deconvolution clusters only
scrna <- deconv_clusters
rm(deconv_clusters)


# TPM normalization 
clusters <- scrna@active.ident

tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

# Read gene lengths from featureCount output
featureCounts <- read.table("featureClength.txt", sep="\t", 
                            header=FALSE)

# TPM if scRNA-seq aligned counts
tpm_counts <- tpm(count_mat, featureCounts$V1)

# Limit to samples present in Seurat object
tpm_counts <- tpm_counts[, match(rownames(scrna@meta.data), 
                                 colnames(tpm_counts))]
# Limit to genes present in Seurat object
tpm_counts <- tpm_counts[match(all.genes, 
                               rownames(tpm_counts)),]


dim(tpm_counts)
# 39539  2344

# Quantile normalization with bulk data
# TPM normalized bulk counts (that have IDH mut / codel status):
bulk_counts <- read.table("bulk_counts.txt", header=TRUE)
bulk_counts <- bulk_counts[,-(1:5)]
sample_names <- colnames(bulk_counts)
sample_info <- c("ID", "Sampleinfo")

# Correct column names
for (i in 1:length(sample_names)) {
  name <- unlist(strsplit(sample_names[i], "[.]"))
  new <- paste(as.character(name[2]), as.character(name[3]), as.character(name[4]), sep=".")
  info <- name[5]
  sample_info <- rbind(sample_info, c(new, info))
}
colnames(sample_info) <- sample_info[1,]
sample_info <- sample_info[-1,]

dim(bulk_counts)
# 39539 170

# Quantile normalization
library(preprocessCore)

# Combine single-cell and bulk sample counts
all_counts <- cbind(tpm_counts, bulk_counts)

# Quantile normalize:
quantile_counts <- normalize.quantiles(as.matrix(all_counts), copy=FALSE)
# log2:
quantile_counts_log2 <- log2(quantile_counts+1)


# Quantile normalized matrix (and log2) of single cell reference:
quantile_log2_single <- quantile_counts_log2[,1:2344]

# Quantile normalized bulk counts:
quantile_log2_GBM <- quantile_counts_log2[,2345:length(colnames(quantile_counts_log2))]

rm(quantile_counts)
rm(quantile_counts_log2)


# Deconvolution MEDIAN SAMPLE from bulk data
median_samples <- c("G17188.TCGA.02.0047.01A.01R.1849.01.2.",
                    "G17189.TCGA.06.0132.01A.02R.1849.01.2.",
                    "G17190.TCGA.06.0174.01A.01R.1849.01.2.",
                    "G17192.TCGA.06.0645.01A.01R.1849.01.2.",
                    "G17193.TCGA.06.0743.01A.01R.1849.01.2.",
                    "G17194.TCGA.02.0055.01A.01R.1849.01.2.",
                    "G17195.TCGA.06.0138.01A.02R.1849.01.2.",
                    "G17196.TCGA.06.0178.01A.01R.1849.01.2.",
                    "G17198.TCGA.06.0646.01A.01R.1849.01.2.",
                    "G17199.TCGA.06.0744.01A.01R.1849.01.2.",
                    "G17801.TCGA.19.4065.02A.11R.2005.01.4.", # IDH codel NA
                    "G17638.TCGA.28.2499.01A.01R.1850.01.2.") # IDH codel NA

median_indices <- which(median_samples %in% colnames(quantile_log2_GBM))
median_indices
gbm_median <- apply(quantile_log2_GBM[,median_indices], 1, median)
gbm_median <- as.data.frame(gbm_median)

# Take medians of samples from same patient
sample_names <- strsplit(colnames(quantile_log2_GBM), "[.]")
patients <- character()
for ( i in 1:170) {
  patients[i] <- sample_names[[i]][4]  
}

indices <- which(duplicated(patients))
indices
for (i in indices) {
  j <- match(patients[i],patients) # first incidence of the gene symbol in j
  quantile_log2_GBM[,j] <- apply(rbind(quantile_log2_GBM[,i], quantile_log2_GBM[,j]), 2, FUN=mean) # change j:th row to median of duplicates
}

# remove extra columns (duplicates) aka samples from same patient
quantile_log2_GBM <- quantile_log2_GBM[,-indices]

dim(quantile_log2_GBM)
# 170 --> 162 samples



# 1 .Median of clusters
# 2. & 3  10% and 25% trimmed mean of clusters
median <- rep.int(0, 39539)
mean10 <- rep.int(0, 39539)
mean25 <- rep.int(0, 39539)

cluster_names <- as.character(clusters)
types <- levels(clusters)
for (i in as.character(types)) {
  Mcluster_i <- apply(quantile_log2_single[,cluster_names[]==i], 1, median)
  median <- cbind(median, Mcluster_i)
  Mcluster_i <- apply(quantile_log2_single[,cluster_names[]==i], 1, mean, trim=0.1)
  mean10 <- cbind(mean10, Mcluster_i)
  Mcluster_i <- apply(quantile_log2_single[,cluster_names[]==i], 1, mean, trim=0.25)
  mean25 <- cbind(mean25, Mcluster_i)
}

median <- median[,-1]
colnames(median) <- paste(as.character(types), "median", sep="_")
mean10 <- mean10[,-1]
colnames(mean10) <- paste(as.character(types), "mean10", sep="_")
mean25 <- mean25[,-1]
colnames(mean25) <- paste(as.character(types), "mean25", sep="_")


#refe <- m
refe <- cbind(mean10, mean25, median)

annotdf <- data.frame('Reference' = rep(c("10% trimmed mean", "25% trimmed mean", "Median"), c(9,9,9)),
                      'Cell type' = c(as.character(types),as.character(types),as.character(types)))
colnames(refe) <- c(paste(as.character(types), "(10% trimmed mean)", sep=" "),
                    paste(as.character(types), "(25% trimmed mean)", sep=" "),
                    paste(as.character(types), "(median)", sep=" "))
rownames(annotdf) <- colnames(refe)
refe <- refe[match(variable, rownames(refe)),] # limit to variable genes

annocol <- list('Reference' = c('10% trimmed mean' = "grey84",
                                '25% trimmed mean' = "gray60",
                                'Median' = "grey35"))

# decrease the number of genes plotted to suitable amount (210)
varrow <- apply(refe, 1, var)
refe2 <- refe
refe2 <- refe[varrow > 5,]
dim(refe2)

library(grid) 
## For pheatmap_1.0.8 and later, rotate axis names:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, 
                 hjust = 1, rot = 55, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

# Choose plotting colors
library(RColorBrewer)
coB <- brewer.pal(9, "Blues")
library(pheatmap)
# Plot pheatmap of reference sets
pheatmap(refe2, annotation_col = annotdf, col=coB, 
         show_rownames = F, show_colnames = F,
         annotation_colors = annocol[1],
         treeheight_row = 0, gaps_col = 3,
         border_color = F)
grid.ls(grid.force())
grid.gedit("col_annotation", gp=gpar(col="black"))




median <- cbind(median, gbm_median)
variable <- read.table("variable_genes.txt", header=FALSE) # 5000 variable genes from Seurat analysis
variable <- variable$V1

gbm <- gbm_median[match(variable, rownames(gbm_median)),]




# DECONVOLUTION (median reference)
# Run deconvolution and compare results between 3 reference sets
colnames(median)[10] <- c("GBM median")
deconvolution_reference <- as.matrix(median[match(variable, rownames(median)),])
deconvolution_mixture <- as.matrix(quantile_log2_GBM[match(variable, rownames(quantile_log2_GBM)),])


# Plot baseline correlation of median results
library(corrplot)
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
deconvolution_reference <- deconvolution_reference[,c("Neuron", "Oligodendrocyte", 
                                                      "OPC", "Astrocyte", "GBM median",
                                                      "Vascular", "Macrophage/Microglia",
                                                      "Mast cell", "Neutrophil", "Dendritic cell")]
corrplot(as.matrix(cor(deconvolution_reference)), method="color",
         type="lower", tl.col = "black", tl.srt = 35, is.corr = FALSE,
         tl.cex = 0.9, order="original", col=col2(50))



# Remove recurrent samples
indices <- c(which(grepl(".02A.",colnames(deconv_res))), which(grepl(".02B.",colnames(deconv_res))))
indices
deconv_res <- deconv_res[,-indices]
dim(deconv_res)
# 154 samples


# Reapeat deconvolution with other reference sets and use the following code for plotting
deconv_all <- deconv_res
deconv_all <- rbind(deconv_all, deconv_res)
deconv2 <- deconv_all

deconv_all <- deconv_all[-c(10,20,30),]
annot_deconv <- data.frame('Reference'= rep(c("10% trimmed mean", 
                                              "25% trimmed mean", 
                                              "Median"), c(9,9,9)),
                           'Cell type'=rep(c(as.character(types)), c(3)))
rownames(deconv_all)[c(10,20,30)] <- c("GBM median (10% trimmed mean)", 
                                       "GBM median (25% trimmed mean)", 
                                       "GBM median (median)")

rownames(annot_deconv) <- rownames(deconv_all)
head(annot_deconv)

annocol <- list('Reference' = c('10% trimmed mean' = "grey84",
                                '25% trimmed mean' = "gray60",
                                'Median' = "grey35"))
jpeg("Deconvolution_Referencesets_14042020.jpg", height = 550, width = 650)
pheatmap(t(deconv_all), annotation_col = annot_deconv, col=coB, 
         show_rownames = F, show_colnames = F,
         annotation_colors = annocol[1],
         treeheight_row = 0, gaps_col = 3,
         border_color = F)
grid.ls(grid.force())
grid.gedit("col_annotation", gp=gpar(col="black"))
dev.off()



deconv_all <- deconv_all[c("Neuron", "Oligodendrocyte", 
                           "OPC", "Astrocyte", "GBM median",
                           "Vascular", "Macrophage/Microglia",
                           "Mast cell", "Neutrophil", "Dendritic cell"),]
jpeg("Deconv_corr_MEAN10_20200319.jpeg", height = 500, width = 500)
corrplot(as.matrix(cor(t(deconv_all))), method="color",
         type="lower", tl.col = "black", tl.srt = 35, is.corr = TRUE,
         tl.cex = 0.9, order="original")
dev.off()


#Deconvolution matrices processing:

# Heatmap of deconvolution results
# Information about IDH codel status needed
# Change sample names to more clear ones

# IDH codel status
gbm_status <- read.table("GBM_IDHcodel2.txt", header=TRUE)

na <- which(is.na(gbm_status$IDH_codel_subtype))
gbm_status <- gbm_status[-na,]
statuses <- c(as.character(gbm_status$IDH_codel_subtype))
length(statuses)

labels <- colnames(deconv_all)
new_names <- rep("", length(labels))
for (i in 1:length(labels)) {
  name <- unlist(strsplit(labels[i], "[.]"))
  new <- paste(as.character(name[2]), as.character(name[3]), 
               as.character(name[4]), sep="-")
  new_names[i] <- as.character(new)
}

ids <- as.character(gbm_status$ID)
# Limit to IDs that have IDH/codel status
deconv_all2 <- deconv_all[,which(new_names %in% ids)]
new_names <- new_names[which(new_names %in% ids)]
colnames(deconv_all2) <- new_names

gbm_status <- gbm_status[which(ids %in% new_names),]
dim(gbm_status)
# 145 2

IDHwt <- gbm_status$ID[which(gbm_status$IDH_codel_subtype == "IDHwt")]
IDHmut_codel <- gbm_status$ID[which(gbm_status$IDH_codel_subtype == "IDHmut-codel")]
IDHmut_non_codel <- gbm_status$ID[which(gbm_status$IDH_codel_subtype == "IDHmut-non-codel")]



IDHwt_i <- match(as.character(IDHwt), colnames(deconv_all2))
#IDHmutcodel_i <- match(IDHmut_codel, colnames(deconv_all))
IDHmutnoncodel_i <- match(as.character(IDHmut_non_codel), colnames(deconv_all2))
c <- rep(1, 145)
c[IDHmutnoncodel_i] <- 2

# Define color for each of the IDHcodel type
colors <- c("#00AFBB", "#E7B800")
colors <- colors[as.numeric(c)]




# Do a heatmap of samples
categories <- c
categories[categories=="1"] <- "IDHwt"
categories[categories=="2"] <- "IDHmut, non-codel"
categories

# add immune group annotation
file_immune <- "..\\immune_groups.txt"
immune_data <- read.table(file_immune, header=FALSE)

immune_annotations <- immune_data$V2[match(colnames(deconv_all2), as.character(immune_data$V1))]
immune_samples <- immune_data$V1[match(colnames(deconv_all2), as.character(immune_data$V1))]

length(immune_annotations)
deconv_all_annot <- deconv_all2[,match(as.character(immune_samples), colnames(deconv_all2))]

head(immune_annotations)
#deconvolution_reference <- as.matrix(m[match(variable, rownames(m)),])

dim(deconv_all_annot)
categories_annot <- categories[match(as.character(immune_samples), colnames(deconv_all2))]



na <- which(is.na(categories_annot))
immune_annotations <- immune_annotations[-na]
deconv_all_annot <- deconv_all_annot[,-na]
categories_annot <- categories_annot[-na]
immune_samples <- immune_samples[-na]

annotdf <- data.frame(
  row.names = colnames(deconv_all_annot), 
  'IDH mutation status' = categories_annot,
  'Immune response' = immune_annotations)

annot_deconv <- data.frame('Reference'= rep(c("10% trimmed mean", "25% trimmed mean", "Median"), c(9,9,9)),
                           'Cell type'=rep(c(as.character(types)), c(3)))
rownames(annot_deconv) <- rownames(deconv_all2)
#head(annot_deconv)

annocol <- list('Cell.type' = c("Oligodendrocyte"="#e496fc",
                                "OPC"="#99c50c",
                                "Neutrophil"="#ff82c4",
                                "Dendritic cell"="#0bcaf7",
                                "Mast cell"="#f6a925",
                                "Vascular"="#fc9785",
                                "Macrophage/Microglia"= "#d5b802",
                                "Astrocyte"="#f884ff",
                                "Neuron"="#04d8e3"),
                'Reference' = c('10% trimmed mean' = "grey84",
                                '25% trimmed mean' = "gray60",
                                'Median' = "grey35"))

library(pheatmap)
library(RColorBrewer)
my_cols <- brewer.pal(9, "Blues")
jpeg("deconv_heatmap_COLannotations_immune_15042020.jpg", height = 900, width = 800)
pheatmap(t(deconv_all_annot), 
         annotation_col= annot_deconv,color=my_cols,
         annotation_colors = annocol,
         clustering_method = "complete",
         fontsize_row = 8, fontsize_col = 18, fontsize=18, show_rownames = F, 
         clustering_distance_rows="minkowski",
         show_colnames = T,
         angle_col = "315", cluster_cols = T, cluster_rows = T, treeheight_row = 0)
grid.ls(grid.force())
grid.gedit("col_annotation.3-3-3-3", gp=gpar(col="black"))
dev.off()


jpeg("deconv_all_pheatmap_median_20200212.jpg", width = 1500, height = 900)
pheatmap(as.matrix(deconv_all), main="")
dev.off()


# BASELINE
library(corrplot)

library(RColorBrewer)
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col1 <- colorRampPalette(c("white", "blue"))

baseline_cor <- cor(t(deconv_all[19:27,]))
dim(baseline_cor)
baseline_cor <- baseline_cor[c(2,1,6,3,4,5,7,8,9), ]
baseline_cor <- baseline_cor[, c(2,1,6,3,4,5,7,8,9)]
rownames(baseline_cor) <- c("OPC","Oligodendrocyte",
                            "Vascular","Neutrophil",          
                            "Dendritic cell",      
                            "Mast cell","Macrophage/Microglia",
                            "Astrocyte","Neuron")
colnames(baseline_cor) <- rownames(baseline_cor)
baseline_cor <- baseline_cor[c(3,6,4,7,1,5,8,2,9), ]
baseline_cor <- baseline_cor[, c(3,6,4,7,1,5,8,2,9)]
# Baseline correlation
library(corrplot)
jpeg("deconv_Baseline_corrplot_MEDIAN_20200415.jpg", height = 600, width = 600)
corrplot(as.matrix(baseline_cor), type= "lower", tl.col = "black", tl.srt = 35,
         method = "color", col=col2(50), is.corr = F, tl.cex = 1.2
)
dev.off()


corrplot(as.matrix(deconv_res))



# Results from Luoto et al. (2018) (suvi)
suvi_res <- read.table("results_regression_no_rec.txt")

# Find matching samples
new_names2 <- rep("", length(colnames(suvi_res)))
for (i in 1:length(colnames(suvi_res))) {
  name <- unlist(strsplit(colnames(suvi_res)[i], "[.]"))
  new <- paste(as.character(name[1]), as.character(name[2]), as.character(name[3]), sep="-")
  new_names2[i] <- new
}

length(which(new_names2 %in% gbm_status$ID))

length(new_names2)
dim(suvi_res)

deconv_res <- deconv_all

colnames(suvi_res) <- new_names2

common <- deconv_res[, which(colnames(deconv_res) %in% colnames(suvi_res))]

write.table(common, "deconv_commonsamples_median_20200415.txt", sep="\t")
common <- read.table("deconv_commonsamples_median_20200415.txt", header=TRUE)
# Correct colnames
new_names2 <- rep("", length(colnames(common)))
for (i in 1:length(colnames(common))) {
  name <- unlist(strsplit(colnames(common)[i], "[.]"))
  new <- paste(as.character(name[1]), as.character(name[2]), as.character(name[3]), sep="-")
  new_names2[i] <- new
}

length(which(new_names2 %in% gbm_status$ID))

length(new_names2)
dim(suvi_res)
# count correlations between both clusters
dim(common)
colnames(common) <- new_names2

suvi_res2 <- suvi_res[, which(colnames(suvi_res) %in% colnames(common))]

dim(suvi_res2)
suvi_res <- suvi_res2

suvi_res <- suvi_res[,colnames(common)]

# Correlation between bulk and scrna
correlations <- rep(0, 10)
library(stats)

for ( i in 1:10) {
  cor_cluster <- rep(0, 10)
  names(cor_cluster) <- rownames(common)
  for ( j in 1:10) {
    #m <- rbind(common[j,], suvi_res[i,])
    cor_cluster[j] <- cor(unlist(suvi_res[i,]), unlist(common[j,]), method = "pearson")
  }
  correlations <- cbind(correlations, cor_cluster)
}

dim(correlations)
correlations <- correlations[,-1]
rownames(correlations) <- rownames(common)
colnames(correlations) <- rownames(suvi_res)[1:10]


write.table(correlations, "deconv_cor_median_20200415.txt", sep = "\t")

correlations2 <- correlations[c(10,9,7,6,1,2,3,4,5,8), c(5,9,1,2,6,7,11,10,3,8,4)]
correlations2 <- correlations[c(5,6,10,2,3,8,4,7,9,1), c(3,4,10,2,6,11,5,9,1,8,7)]
n <- lapply(strsplit(rownames(correlations2), "[(]"), function(x) x[1])
pheatmap(as.matrix(correlations2),  main="scRNA deconvolution correlation to bulk deconvolution", 
         cluster_cols = T, cluster_rows = T, col=brewer.pal(11, "PRGn"), angle_col =315,
         fontsize = 12, display_numbers = T,fontsize_number = 8)

library(grid)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, 
                      height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(as.matrix(correlations2),  
         main="scRNA deconvolution correlation to bulk deconvolution", 
         cluster_cols = T, cluster_rows = T, col=brewer.pal(11, "RdBu"), angle_col =315,
         fontsize = 25, display_numbers = T, number_color = "gray0" ,
         fontsize_number = 20, xlab="Bulk deconvolution", 
         ylab="scRNA deconvolution")
setHook("grid.newpage", NULL, "replace")
grid.text("Bulk deconvolution", x=0.35, y=-0.01, gp=gpar(fontsize=30))
grid.text("scRNA deconvolution", x=0.97, rot=90, gp=gpar(fontsize=30))


library(corrplot)
corrplot(as.matrix(correlations), method = "square", is.corr = T, order="hclust")

rownames(common) <- n
# Compare suvi's results with cmdscale
mds <- cmdscale(dist(rbind(suvi_res[-5,], common[-10,]), method="euclidean"))
colnames(mds) <- c("var1", "var2")
# Define color for each of the 3 iris species
colors <- c("#00AFBB", "#E7B800")
display.brewer.all(colorblindFriendly = TRUE)
colors <- brewer.pal(10, "Paired")
colors <- colors[c(4,10)]
colors <- colors[c(rep(1,10),rep(2,9))]


# Define shapes
shapes = c(15, 17) 
shapes <- shapes[c(rep(1,10),rep(2,9))]
# Plot bulk and single-cell deconvolution MDS
plot(mds, col = colors, pch = shapes, xlab = "", ylab = "", asp = 1, #axes = FALSE,
     main="Bulk & scRNA-seq deconvolution dissimilarities", cex=2, cex.main=3, 
     ylim=c(-0.6, 0.8), cex.axis = 2)
text(mds, labels = row.names(mds), pos = 1, col = colors, pch = shapes, cex=2)
legend("topright", legend = c("Bulk RNA-seq reference", "scRNA-seq reference"),
       col =  c("#00AFBB", "#E7B800"),
       pch = c(16, 17), cex=2)


install.packages("ggplot2")
library(ggplot2)
library(gplots)
library(ggrepel)

reference <- c(rep("Bulk RNA-seq", 10), 
               rep("scRNA-seq", 9))
mds <- as.data.frame(mds)
mds$label <- rownames(mds)
mds$label[13] <- "Neutrophil"
mds$reference <- reference
#mds$other <- reference
mds$label[c(13, 19)] <- c("Macrophage/Microglia", "Neuron")

ggplot(mds, aes(x=var1, y=var2, label=label, color=factor(reference))) +
  scale_color_manual(values=c("#33A02C", "#6A3D9A")) +
  geom_point(size=5) +
  geom_text_repel(size=10, xlim = c(-0.55,1.25), show.legend = FALSE) +
  scale_y_continuous(breaks=seq(-0.5, 0.7, 0.1)) +         # Set tick every 0.05
  scale_x_continuous(breaks=seq(-0.65, 1.25, 0.1)) +         # Set tick every 0.05
  theme_bw(base_size = 20) +
  theme(legend.title=element_text(size=24, face = "bold"), 
        legend.text=element_text(size=22),
        axis.title = element_text(size=22, face = "bold"),
        axis.text = element_text(size=22),
        axis.text.x = element_text(angle = 90)) +
  
  labs(title="Metric MDS", y="Coordinate 1", x = "Coordinate 2", color="Deconvolution reference")




# Plot single-cell and bulk deconvolution results in heatmap
# 14.4.2020

suvi_median <- rbind(suvi_res[-5,], common[-10,])


# IDH codel status
gbm_status <- read.table("..\\..\\GBM_IDHcodel2.txt", header=TRUE)

na <- which(is.na(gbm_status$IDH_codel_subtype))
gbm_status <- gbm_status[-na,]
statuses <- c(as.character(gbm_status$IDH_codel_subtype))
length(statuses)


ids <- as.character(gbm_status$ID)
new_names <- colnames(suvi_median)
# Limit to IDs that have IDH/codel status
suvi_median <- suvi_median[,which(new_names %in% ids)]
new_names <- new_names[which(new_names %in% ids)]
colnames(suvi_median) <- new_names

gbm_status <- gbm_status[which(ids %in% new_names),]
dim(gbm_status)
# 134 2

IDHwt <- gbm_status$ID[which(gbm_status$IDH_codel_subtype == "IDHwt")]
IDHmut_codel <- gbm_status$ID[which(gbm_status$IDH_codel_subtype == "IDHmut-codel")]
IDHmut_non_codel <- gbm_status$ID[which(gbm_status$IDH_codel_subtype == "IDHmut-non-codel")]



IDHwt_i <- match(IDHwt, colnames(deconv_all))
#IDHmutcodel_i <- match(IDHmut_codel, colnames(deconv_all))
IDHmutnoncodel_i <- match(IDHmut_non_codel, colnames(deconv_all))
c <- rep(1, 134)
c[IDHmutnoncodel_i] <- 2

# Define color for each of the IDHcodel type
colors <- c("#00AFBB", "#E7B800")
colors <- colors[as.numeric(c)]




# Do a heatmap of samples
categories <- c
categories[categories=="1"] <- "IDHwt"
categories[categories=="2"] <- "IDHmut, non-codel"
categories

# add immune group annotation
file_immune <- "..\\immune_groups.txt"
immune_data <- read.table(file_immune, header=FALSE)

immune_annotations <- immune_data$V2[match(colnames(suvi_median), as.character(immune_data$V1))]
immune_samples <- immune_data$V1[match(colnames(suvi_median), as.character(immune_data$V1))]

length(immune_annotations)
deconv_all_annot <- suvi_median[,match(as.character(immune_samples), colnames(suvi_median))]

head(immune_annotations)
#deconvolution_reference <- as.matrix(m[match(variable, rownames(m)),])

dim(deconv_all_annot)
categories_annot <- categories[match(as.character(immune_samples), colnames(suvi_median))]



na <- which(is.na(categories_annot))
immune_annotations <- immune_annotations[-na]
deconv_all_annot <- deconv_all_annot[,-na]
categories_annot <- categories_annot[-na]
immune_samples <- immune_samples[-na]


deconv_sc <- deconv_all_annot


deconv2 <- deconv_all_annot
deconv_all_annot <- deconv2
rownames(deconv_all_annot)[c(11:19)] <- c("Oligodendrocyte","OPC",               
                                          "Neutrophil","Dendritic cell",     
                                          "Mast cell","Vascular",            
                                          "sc Macrophage/Microglia","Astrocyte", "Neuron (scRNA)")


annotdf <- data.frame(
  row.names = colnames(deconv_all_annot), 
  'IDH mutation status' = categories_annot,
  'Immune response' = immune_annotations)

annot_ref <- data.frame(row.names=rownames(deconv_all_annot),
                        'Reference'=rep(c("Bulk", "Single-cell"), c(10,9)))

annocol <- list('Immune.response' = c('Cellular-like' = "#ff928b",
                                      'Humoral' = "#00cbfd",
                                      'Negative' = "#c0c200"),
                'Reference'=c('Bulk'="gray84", 'Single-cell'="grey35"))

library(pheatmap)
library(RColorBrewer)
my_cols <- brewer.pal(9, "Blues")
pheatmap(t(deconv_all_annot), annotation_row = annotdf, 
         annotation_col = annot_ref, color=my_cols,
         annotation_colors = annocol,
         gaps_col = 10,
         clustering_method = "ward.D",
         fontsize_row = 8, fontsize_col = 18, fontsize=18, show_rownames = F, 
         clustering_distance_rows="minkowski",
         annotation_names_row = F,
         angle_col = "315", cluster_cols = F, cluster_rows = T, treeheight_row =0, treeheight_col = 15)
