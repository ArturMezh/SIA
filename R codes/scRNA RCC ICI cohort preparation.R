setwd("/Volumes/Samsung_X5/Single Cell/RCC_ICI/")

## Load Seurat library

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(SeuratObject)

counts <-NULL

setwd("/Volumes/Samsung_X5/Single Cell/RCC_ICI/data/")


counts <- read.delim("ccRCC_scRNASeq_NormalizedCounts.txt.gz", sep = "\t", check.names = T)
metadata  <- data.table::fread("Final_SCP_Metadata.txt")

metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$NAME

head(counts)
rownames(counts) <- counts$GENE
counts.2 <- counts[,-1]


obj <- CreateSeuratObject(counts = counts.2, min.cells = 3, min.features = 200, meta.data = metadata)

                        
                        merged_seurat <- obj
# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

merged_seurat <- subset(merged_seurat, subset = orig.ident != "NA")

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data
head(metadata)
# Add cell IDs to metadata

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = biosample_id, fill=biosample_id)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=biosample_id, x=mitoRatio, fill=biosample_id)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

str(merged_seurat@meta.data)
# Filter out low quality cells using selected thresholds - these will change with experiment

filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 100) & 
                            (nGene >= 100) & 
                            (log10GenesPerUMI > 0.08) & 
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

seurat_phase <- filtered_seurat




# Load cell cycle markers
setwd("/Volumes/Samsung_X5/Single Cell/RCC_ICI/")
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
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))



# Split seurat object by condition to perform cell cycle scoring and SCT on all samples

head(seurat_phase)
split_seurat <- SplitObject(seurat_phase, split.by = "biosample_id")
head(split_seurat)


setwd("/Volumes/Samsung_X5/Single Cell/RCC_ICI/")
saveRDS(seurat_phase,  "immune.combined RCC_ICI.rds")










