# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)



setwd("/Volumes/Samsung_X5/Single Cell/GSE139829_MEL/GSE139829_RAW/")
counts <-NULL
counts <-readMM("GSM4147091_BSSR0022_matrix.mtx.gz")
genes <- read_tsv("GSM4147091_BSSR0022_genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv("GSM4147091_BSSR0022_barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts <- as(counts, "dgCMatrix")
counts <- CreateSeuratObject(counts = counts,
                                    min.features = 100)

counts.BSSR0022 <- counts



counts <-readMM("GSM4147092_UMM041L_matrix.mtx.gz")
genes <- read_tsv("GSM4147092_UMM041L_genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv("GSM4147092_UMM041L_barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts <- as(counts, "dgCMatrix")
counts <- CreateSeuratObject(counts = counts,
                             min.features = 100)
counts.UMM041L <- counts



counts <-readMM("GSM4147093_UMM059_matrix.mtx.gz")
genes <- read_tsv("GSM4147093_UMM059_genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv("GSM4147093_UMM059_barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts <- as(counts, "dgCMatrix")
counts <- CreateSeuratObject(counts = counts,
                             min.features = 100)
counts.UMM059 <- counts

  
counts <-readMM("GSM4147094_UMM061_matrix.mtx.gz")
genes <- read_tsv("GSM4147094_UMM061_genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv("GSM4147094_UMM061_barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts <- as(counts, "dgCMatrix")
counts <- CreateSeuratObject(counts = counts,
                             min.features = 100)
counts.UMM061 <- counts


counts <-readMM("GSM4147095_UMM062_matrix.mtx.gz")
genes <- read_tsv("GSM4147095_UMM062_genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv("GSM4147095_UMM062_barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts <- as(counts, "dgCMatrix")
counts <- CreateSeuratObject(counts = counts,
                             min.features = 100)
counts.UMM062 <- counts


counts <-readMM("GSM4147096_UMM063_matrix.mtx.gz")
genes <- read_tsv("GSM4147096_UMM063_genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv("GSM4147096_UMM063_barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts <- as(counts, "dgCMatrix")
counts <- CreateSeuratObject(counts = counts,
                             min.features = 100)
counts.UMM063 <- counts


counts <-readMM("GSM4147097_UMM064_matrix.mtx.gz")
genes <- read_tsv("GSM4147097_UMM064_genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv("GSM4147097_UMM064_barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts <- as(counts, "dgCMatrix")
counts <- CreateSeuratObject(counts = counts,
                             min.features = 100)
counts.UMM064 <- counts



counts <-readMM("GSM4147098_UMM065_matrix.mtx.gz")
genes <- read_tsv("GSM4147098_UMM065_genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv("GSM4147098_UMM065_barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts <- as(counts, "dgCMatrix")
counts <- CreateSeuratObject(counts = counts,
                             min.features = 100)
counts.UMM065 <- counts



counts <-readMM("GSM4147099_UMM066_matrix.mtx.gz")
genes <- read_tsv("GSM4147099_UMM066_genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv("GSM4147099_UMM066_barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts <- as(counts, "dgCMatrix")
counts <- CreateSeuratObject(counts = counts,
                             min.features = 100)
counts.UMM066 <- counts


counts <-readMM("GSM4147101_UMM069_matrix.mtx.gz")
genes <- read_tsv("GSM4147101_UMM069_genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv("GSM4147101_UMM069_barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts <- as(counts, "dgCMatrix")
counts <- CreateSeuratObject(counts = counts,
                             min.features = 100)
counts.UMM069 <- counts



counts <-readMM("GSM4147100_UMM067L_matrix.mtx.gz")
genes <- read_tsv("GSM4147100_UMM067L_genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2
cell_ids <- read_tsv("GSM4147100_UMM067L_barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
counts <- as(counts, "dgCMatrix")
counts <- CreateSeuratObject(counts = counts,
                             min.features = 100)
counts.UMM067L <- counts

##

merged_seurat <- merge(counts.BSSR0022, y = c(counts.UMM041L,counts.UMM059, counts.UMM061,counts.UMM062,counts.UMM063,counts.UMM064,counts.UMM065,
                                         counts.UMM066,counts.UMM069,counts.UMM067L), 
                  add.cell.ids = c("BSSR0022","UMM041L","UMM059","UMM061","UMM062","UMM063","UMM064","UMM065","UMM066","UMM069","UMM067L"), project = "GSE139829_MEL")

rm(counts.BSSR0022)
rm(counts.UMM041L)
rm(counts.UMM059)
rm(counts.UMM061)
rm(counts.UMM062)
rm(counts.UMM063)
rm(counts.UMM064)
rm(counts.UMM065)
rm(counts.UMM066)
rm(counts.UMM069)
rm(counts.UMM067L)

library(SeuratObject)



# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

metadata$sample <- NA
metadata$sample <- metadata$cells
metadata$sample <- gsub('_.*', '', metadata$sample)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)



# Single-cell RNA-seq - normalization
# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)
# Load cell cycle markers
load("cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
#View(seurat_phase@meta.data)  

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

# Plot the PCA colored by cell location
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "sample",
        split.by = "Phase")

# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))

# Plot the PCA colored by mitoFr
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

# normalize and identify variable features for each dataset independently
split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat)


anchors <- FindIntegrationAnchors(object.list = split_seurat, anchor.features = integ_features)
# this command creates an 'integrated' data assay
seurat_integrated <- IntegrateData(anchorset = anchors)

immune.combined <- seurat_integrated

# unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)



# Create .RData object to load at any time
setwd("/Volumes/Samsung_X5/Single Cell/GSE139829_MEL")
saveRDS(immune.combined, file="GSE139829_MEL_immune.combined.rds")