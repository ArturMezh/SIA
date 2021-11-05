library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(SeuratObject)
library(GEOquery)
library(Biobase)

counts <-NULL

setwd("/Volumes/Samsung_X5/Single Cell/MEL_ICI/data/")


counts <- read.table(file = "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt",  sep="\t", header=F, fill = TRUE) 

counts<-counts[-c(2), ]
names(counts) <- counts[1,]
counts <- counts[-1,]
rownames(counts)<- counts[,1]
counts<- counts[,-1]
counts  <- CreateSeuratObject(counts = counts, project = "MEL_ICI")

meta <- data.table::fread("GSE120575_patient_ID_single_cells.txt")
meta <- meta[20:nrow(meta),]
meta <- meta[1:16293,]
colnames(meta) <- as.character(meta[1,])
meta <- meta[-1, ]

merged_seurat <- counts

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

metadata <- merged_seurat@meta.data
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# add meta data
meta <- meta[,c(1:11)]
metadata <- cbind(metadata , meta, by="title")


# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
merged_seurat <- subset(merged_seurat, subset = orig.ident != "NA")

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata <- merged_seurat@meta.data
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = `characteristics: response`)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.78)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=`characteristics: response`, x=mitoRatio, fill=`characteristics: response`)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 100) & 
                            (log10GenesPerUMI > 0.78) & 
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

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.0085996, 0.0105469, 0.0126909, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))


# Plot the PCA colored by mitoFr
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")


head(seurat_phase)
split_seurat <- SplitObject(seurat_phase, split.by = "characteristics..response")


# normalize and identify variable features for each dataset independently
split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

rm(counts)
rm(filtered_counts)
rm(nonzero)
rm(filtered_seurat)
rm(seurat_phase)
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
setwd("/Volumes/Samsung_X5/Single Cell/MEL_ICI")
saveRDS(immune.combined, file="MEL_ICI_immune.combined.rds")









