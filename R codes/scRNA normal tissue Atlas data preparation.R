# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library("remotes")
library(SeuratObject)

counts <-NULL
setwd("/Volumes/Samsung_X5/Single Cell/ATLAS/GSE159929_RAW/")
#CN1

counts <- read.csv(file = "GSM4850577_Bladder_Counts.csv.gz", header = TRUE, row.names = 1)
Bladder  <- CreateSeuratObject(counts = counts, project = "Bladder", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850578_Blood_Counts.csv.gz", header = TRUE, row.names = 1)
Blood  <- CreateSeuratObject(counts = counts, project = "Blood", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850579_Common.bile.duct_Counts.csv.gz", header = TRUE, row.names = 1)
bile.duct  <- CreateSeuratObject(counts = counts, project = "bile.duct", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850580_Esophagus_Counts.csv.gz", header = TRUE, row.names = 1)
Esophagus  <- CreateSeuratObject(counts = counts, project = "Esophagus", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850581_Heart_Counts.csv.gz", header = TRUE, row.names = 1)
Heart  <- CreateSeuratObject(counts = counts, project = "Heart", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850582_Liver_Counts.csv.gz", header = TRUE, row.names = 1)
Liver  <- CreateSeuratObject(counts = counts, project = "Liver", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850583_Lymph.node_Counts.csv.gz", header = TRUE, row.names = 1)
Lymph.node  <- CreateSeuratObject(counts = counts, project = "Lymph.node", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850584_Marrow_Counts.csv.gz", header = TRUE, row.names = 1)
Marrow_Counts  <- CreateSeuratObject(counts = counts, project = "Marrow_Counts", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850585_Muscle_Counts.csv.gz", header = TRUE, row.names = 1)
Muscle  <- CreateSeuratObject(counts = counts, project = "Muscle", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850586_Rectum_Counts.csv.gz", header = TRUE, row.names = 1)
Rectum  <- CreateSeuratObject(counts = counts, project = "Rectum", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850587_Skin_Counts.csv.gz", header = TRUE, row.names = 1)
Skin_Counts  <- CreateSeuratObject(counts = counts, project = "Skin_Counts.duct", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850588_Small.intestine_Counts.csv.gz", header = TRUE, row.names = 1)
Small.intestine  <- CreateSeuratObject(counts = counts, project = "Small.intestine", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850589_Spleen_Counts.csv.gz", header = TRUE, row.names = 1)
Spleen  <- CreateSeuratObject(counts = counts, project = "Spleen", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850590_Stomach_Counts.csv.gz", header = TRUE, row.names = 1)
Stomach  <- CreateSeuratObject(counts = counts, project = "Stomach", min.cells = 3, min.features = 100)

counts <- read.csv(file = "GSM4850591_Trachea_Counts.csv.gz", header = TRUE, row.names = 1)
Trachea  <- CreateSeuratObject(counts = counts, project = "Trachea", min.cells = 3, min.features = 100)

##

merged_seurat <- merge(Bladder, y = c(Blood,
                                      bile.duct, 
                                      Esophagus,
                                      Heart,
                                      Liver,
                                      Lymph.node,
                                      Marrow_Counts,
                                      Muscle,
                                      Rectum,
                                      Skin_Counts,
                                      Small.intestine,
                                      Spleen,
                                      Stomach,
                                      Trachea), project = "ATLAS")


rm(Bladder)
rm(Blood)
rm(bile.duct)
rm(Esophagus)
rm(Heart)
rm(Liver)
rm(Lymph.node)
rm(Marrow_Counts)
rm(Muscle)
rm(Rectum)
rm(Skin_Counts)
rm(Small.intestine)
rm(Spleen)
rm(Stomach)
rm(Trachea)



# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

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
metadata$sample <- gsub('_.*', '', metadata$orig.ident)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
head(merged_seurat)


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

#split_seurat <- split_seurat[c("tumour core", "tumour border", "normal tissue adjacent to neoplasm")]

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
setwd("/Volumes/Samsung_X5/Single Cell/ATLAS")
saveRDS(immune.combined, file="ATLAS_immune.combined.rds")