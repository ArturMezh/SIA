########################
# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)



setwd("/Volumes/Samsung_X5/Single Cell/E-MTAB-8410/non normalised/")



# Read in `matrix.mtx`
counts <-NULL
counts <-readMM("E-MTAB-8410.aggregated_filtered_counts.mtx")
# Read in `genes.tsv`
genes <- read_tsv("E-MTAB-8410.aggregated_filtered_counts.mtx_rows", col_names = FALSE)
gene_ids <- genes$X1
# Read in `barcodes.tsv`
cell_ids <- read_tsv("E-MTAB-8410.aggregated_filtered_counts.mtx_cols", col_names = FALSE)$X1
# Create a sparse matrix for more efficient computation
counts <- as(counts, "dgCMatrix")

# Make the column names as the cell IDs and the row names as the gene IDs
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids

                                    # Connect to AnnotationHub
                                    library(AnnotationHub)
                                    ah <- AnnotationHub()
                                    # Access the Ensembl database for organism
                                    ahDb <- query(ah, 
                                                  pattern = c("Homo sapiens", "EnsDb"), 
                                                  ignore.case = TRUE)
     
                                    # Check versions of databases available
                                    ahDb %>% 
                                      mcols()
                                    # Acquire the latest annotation files
                                    id <- ahDb %>%
                                      mcols() %>%
                                      rownames() %>%
                                      tail(n = 1)                
                                    # Download the appropriate Ensembldb database
                                    edb <- ah[[id]]                
                                    
                                    # Extract gene-level information from database
                                    annotations <- genes(edb, 
                                                         return.type = "data.frame")
                  

# replace cosed by genes
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name)



gID <- as.data.frame(gene_ids)
gID$gene <- NA
gID$gene <- gID$gene_ids
gID$gene <- as.character(annotations$gene_name)[match(unlist(gID$gene), annotations$gene_id)]
n_occur <- data.frame(table(gID$gene))
gID[gID$gene %in% n_occur$Var1[n_occur$Freq > 1],]
gID<- gID[!duplicated(gID$gene), ] 
rownames(counts) <- as.character(gID$gene)[match(unlist(rownames(counts)), gID$gene_ids)]
rownames(counts) <- make.names(rownames(counts), unique = TRUE)                        


# Turn count matrix into a Seurat object (output is a Seurat object)
merged_seurat <- CreateSeuratObject(counts = counts,
                           min.features = 100)

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

# read metadata 
meta <- data.table::fread("E-MTAB-8410.sdrf.txt")
# Create sample column
metadata$sample <- NA
metadata$sample <- metadata$cells
metadata$sample <- gsub('-.*', '', metadata$sample)
metadata$cells <- gsub('-.*', '', metadata$sample)
metadata$region <- gsub('-.*', '', metadata$sample)
metadata$region <- as.character(meta$`Characteristics[sampling site]`)[match(unlist(metadata$sample), meta$`Comment[BioSD_SAMPLE]`)]

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

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)
# Load cell cycle markers
setwd("/Volumes/Samsung_X5/Single Cell/E-MTAB-8410/")
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

head(seurat_phase)
split_seurat <- SplitObject(seurat_phase, split.by = "region")
head(split_seurat)

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

# save the data

setwd("/Volumes/Samsung_X5/Single Cell/E-MTAB-8410/")
saveRDS(immune.combined, "/Volumes/Samsung_X5/Single Cell/E-MTAB-8410/immune.combined CRC.rds")


