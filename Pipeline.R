#####
# Load libraries
#####
library(BiocManager)
library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(readr)
library(patchwork)
library(scater)
library(SingleR)
library(celldex)
library(scuttle)
library(enrichR)


# ---
# Create Seurat Object

# Reads in matrix.mtx
count_matrix <- readMM('~/Desktop/GVHD-SC-Analysis/RawData\ /MWT-BM/matrix.mtx.gz')

# Reads in features.tsv
genes <- read_tsv('~/Desktop/GVHD-SC-Analysis/RawData\ /MWT-BM/features.tsv.gz', col_names = FALSE, show_col_types = FALSE)
gene_ids <- genes$X2

# Reads in barcodes.tsv
cell_ids <- read_tsv('~/Desktop/GVHD-SC-Analysis/RawData\ /MWT-BM/barcodes.tsv.gz', col_names = FALSE, show_col_types = FALSE)$X1

# Add row names (the gene_ids) and column names (the cell_ids) to the count data
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- cell_ids


seurat_object = CreateSeuratObject(counts = count_matrix, project = "scRNA-seq", assay = "RNA")
View(seurat_object@meta.data)

# ---
# QC and Cell Selection

# Mitochondrial Ratio
total_counts_per_cell <- colSums(seurat_object@assays$RNA@counts)
mito_genes <- rownames(seurat_object)[grep("^mt-", rownames(seurat_object))]
seurat_object$percent_mito <- colSums(seurat_object@assays$RNA@counts[mito_genes, ])/total_counts_per_cell
head(mito_genes, 20)

# Plot Unfiltered QC: Violin Plot
VlnPlot(seurat_object, features = "nFeature_RNA", assay = "RNA",log = TRUE, pt.size = 0) + geom_boxplot(width=0.1,fill="white")
VlnPlot(seurat_object, features = "nCount_RNA", assay = "RNA",log = TRUE, pt.size = 0) + geom_boxplot(width=0.1,fill="white")

# Plot Unfiltered QC: Scatter Plots
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Subset with thresholds
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent_mito < 5)
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mito < 5)
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filtered visuals: Violin
VlnPlot(seurat_object, features = "nFeature_RNA", assay = "RNA",log = TRUE, pt.size = 0) + geom_boxplot(width=0.1,fill="white")
VlnPlot(seurat_object, features = "nCount_RNA", assay = "RNA",log = TRUE, pt.size = 0) + geom_boxplot(width=0.1,fill="white")

# Filtered visuals: Scatter
FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# ---
# Normalization & Scaling

# Data Normalization
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Feature Selection
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# Identify top 20 most variable genes
top10 <- head(VariableFeatures(seurat_object), 10)

# Plot Variable Features
plot1 <- VariableFeaturePlot(seurat_object)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Scale Data
# (shifts exp of each gene so mean across cells is 0 and variance across cells is 1)
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)


# ---
# Dimensionality

# Linear Dimensional Reduction
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)

# Represent PCA
VizDimLoadings(seurat_object, dims = 1, reduction = "pca")
VizDimLoadings(seurat_object, dims = 2, reduction = "pca")
VizDimLoadings(seurat_object, dims = 3, reduction = "pca")
VizDimLoadings(seurat_object, dims = 4, reduction = "pca")
VizDimLoadings(seurat_object, dims = 1:4, reduction = "pca")

DimPlot(seurat_object, reduction = "pca")

DimHeatmap(seurat_object, dims = 1:2, cells = 500, balanced = TRUE)
DimHeatmap(seurat_object, dims = 3:4, cells = 500, balanced = TRUE)

# Determine "Dimensionality"
# (JackStraw procedure: permutes a subset of data and reruns PCA consturcting a 
# null distribution of feature scores; identifies significant PCs as those with 
# strong enrichment of low p-values features)
seurat_object <- JackStraw(seurat_object, num.replicate = 100)
seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20)

# Jack Straw Plot
# (compares distribution of p-values for each PC with a uniform distribution)
JackStrawPlot(seurat_object, dims = 1:20)

# Elbow Plot
ElbowPlot(seurat_object)


# ---
# Clustering

# Find Euclidean distance in PCA space
seurat_object <- FindNeighbors(seurat_object, dims = 1:15)

# Modularity Optimization Technique to iteratively group cells together
seurat_object <- FindClusters(seurat_object, resolution = 1)   ###standard resolution is 1 but can modify to change clusters 

# Cluster IDs of first 5 cells
head(Idents(seurat_object), 5)


# ---
# Non-Linear Dimensional Reduction

# UMAP non-linear dimensional reduction technique
seurat_object <- RunUMAP(seurat_object, dims = 1:15)

# Plot Clusters
DimPlot(seurat_object, reduction = "umap")
DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# ---
# Differentially Expressed Features (Cluster Biomarkers)

# Find markers of cluster 1
cluster1.markers <- FindMarkers(seurat_object, ident.1 = 1, min.pct = 0.25)

# Find markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(seurat_object, ident.1 = 5, ident.2 = c(0,3), 
                                min.pct = 0.25)

# Find markers for every cluster compared to all remaining cells, reporting only positive markers ##Slow Step 
# Use FindConservedMarkers() for two different conditions
seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, 
                                        min.pct = 0.25, logfc.threshold = 0.25)
seurat_object.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# To see the markers for each cluster:
seurat_object.markers
head(seurat_object.markers)

# ROC Test for Differential Expression
# (returns classification power for any individual marker, ranging from
# 0-random to 1-perfect)
cluster0.markers <- FindMarkers(seurat_object, ident.1 = 0, logfc.threshold = 0.25,
                                test.use = "roc", only.pos = TRUE)

# ---
# Visualizing Marker Expression

# Do Heatmap- expression heatmap for cells and features (plot the top 10 markers
# for each cluster)
seurat_object.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seurat_object, features = top10$gene) + NoLegend()

# ---
# Cell Type Labeling

# Obtain reference for annotation of cell types
mouse.se <- MouseRNAseqData()
reference <- mouse.se[,!is.na(mouse.se@colData@listData[["label.main"]])]

# convert gene names to uppercase
seurat_object@assays[["RNA"]]@counts@Dimnames[[1]] <- toupper(seurat_object@assays[["RNA"]]@counts@Dimnames[[1]])

# designate input data for SingleR()  (SingleR takes in a normalized count matrix with cells vs genes)
for_labels <- as.data.frame(seurat_object@assays[["RNA"]]@counts)

# annotate using SingleR()
annotatedcelltype <- SingleR(test = for_labels, ref = reference, labels = reference@colData@listData[["label.main"]])

# to look at labels for cell types:
View(table(annotatedcelltype@listData[["labels"]]))

# to add labels to Seurat object:
seurat_object[["SingleR.labels"]] <- annotatedcelltype@listData[["pruned.labels"]]

# (mainly for visualizing) rename clusters to their cell types
cell.type.freq <- table(annotatedcelltype@listData[["labels"]])
view(cell.type.freq)
Idents(seurat_object)

seurat_object.final <- SetIdent(seurat_object, value = 'SingleR.labels')

# Top 10 Genes for T cells 
seurat_object.final.markers <- FindAllMarkers(seurat_object.final, only.pos = TRUE, 
                                        min.pct = 0.25, logfc.threshold = 0.25)
seurat_object.final.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

top.10.genes.Tcells <- filter(seurat_object.final.markers, cluster == "T cells")

top.genes.plot <- ggplot(top.10.genes.Tcells, aes(x=reorder(gene, -avg_log2FC), y=avg_log2FC, fill = gene)) + 
  geom_bar(stat = "identity", show.legend = FALSE)+
  labs(title ="Top 10 Genes for T cells in MWT-BM", x= "Genes", y= "avg_log2FC")+ 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(l =50))

top.genes.plot

png(filename = "/Users/garciaduttonn2/Desktop/GVHD-SC-Analysis/FinalFigures\ /MWT-BM-topgenes.png", width = 1000, height = 720)
print(top.genes.plot)
dev.off()

# Produce UMAP dimension reduction plots labeled with cell types (instead of clusters)
DimPlot(seurat_object.final, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(seurat_object.final, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) 

png(filename = "/Users/garciaduttonn2/Desktop/GVHD-SC-Analysis/FinalFigures\ /MWT-SPL-UMAP.png", width = 1000, height = 720)
print(DimPlot(seurat_object.final, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE))
dev.off()

# Create Plots of Cell Types 
pt <- table(Idents(seurat_object.final), seurat_object.final$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
 cell.type.plot <- ggplot(pt, aes(x=Var1, y=Freq, fill = Var1)) + 
  geom_bar(stat = "identity", show.legend = FALSE)+
  labs(title ="Cell Type Frequency", x= "Cell Type", y= "Frequency")+ 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(l =50))

pt$Percentage <- pt$Freq/sum(pt$Freq)*100
pt$Round_off <- round(pt$Percentage, digits = 2)
 
pie.chart <- ggplot(pt, aes(x = "", y = -Round_off, 
                             fill = reorder(Var1, -Round_off))) + 
   geom_bar(stat = "identity", color = "black") + 
   labs(title = "Cell Type Frequency Share for MWT-SPL", fill = "Cell Type") +
   coord_polar("y") +
   theme_void()
 
pie.chart 
 
write.csv(pt, "/Users/garciaduttonn2/Desktop/GVHD-SC-Analysis/Tables/MWT-SPL-Percentages.csv")

png(filename = "/Users/garciaduttonn2/Desktop/GVHD-SC-Analysis/FinalFigures\ /MWT-SPL-PieChart.png", width = 1000, height = 720)
print(pie.chart)
dev.off()

png(filename = "/Users/garciaduttonn2/Desktop/GVHD-SC-Analysis/FinalFigures\ /FHET-SPL-celltypeplot.png", width = 1000, height = 720)
print(cell.type.plot)
dev.off()


# DE and EnrichR pathway visualization barplot 
gran.plot <- DEenrichRPlot(
  seurat_object.final,
  ident.1 = 'Granulocytes',
  ident.2 = 'T cells',
  balanced = TRUE,
  logfc.threshold = 0.25,
  assay = NULL,
  max.genes = 100,
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2021",
  num.pathway = 10,
  return.gene.list = FALSE)

png(filename = "/Users/garciaduttonn2/Desktop/GVHD-SC-Analysis/FinalFigures\ /FHET-SPL-Granulocytes-DEenrichPlot.png", width = 1000, height = 720)
print(gran.plot)
dev.off()

tcell.plot <- DEenrichRPlot(
  seurat_object.final,
  ident.1 = 'T cells',
  ident.2 = 'Granulocytes',
  balanced = TRUE,
  logfc.threshold = 0.25,
  assay = NULL,
  max.genes = 100,
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = "GO_Biological_Process_2021",
  num.pathway = 10,
  return.gene.list = FALSE)

png(filename = "/Users/garciaduttonn2/Desktop/GVHD-SC-Analysis/FinalFigures\ /FHET-SPL-Tcell-DEenrichPlot.png", width = 1000, height = 720)
print(tcell.plot)
dev.off()


# ---
# NOTES: EXPORTS

#       Exporting visuals:
# png(filename = "[insert filepath here/visualname.png]", width = [num], height = [num])
# print([insert call for the visual itself**])
# dev.off()

#       Exporting tables:
# write.csv([object to export], file = "[insert filepath here/tablename.csv])



# ** e.g. print(DimPlot(obj, red, lab, pt.size))
