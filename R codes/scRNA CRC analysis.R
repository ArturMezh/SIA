########################

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(metap)
library(patchwork)
library(ggplot2)
library(data.table)


setwd("/Volumes/Samsung_X5/Single Cell/E-MTAB-8410/")


# Load the split seurat object into the environment
immune.combined <- readRDS("immune.combined CRC region_split.rds")

DefaultAssay(immune.combined) <- "RNA"

# cell number
cell.num <- table(immune.combined@active.ident)
cell.num <- as.data.frame(cell.num)
sum(cell.num$Freq)

# for supplementary figure
FeaturePlot(immune.combined, features = c( "CD68", "CD163",  "C1QA", "C1QB","C1QC","APOE"), split.by = "region",  min.cutoff = "q5", cols = c( "gray90","red"))& 
  theme(text = element_text(face = "bold"),
        axis.text.x=element_text(angle=45, hjust=1, size=8),
        axis.text.y=element_text( size=8),
        axis.title =   element_blank(),
        axis.title.y.right = element_text(size = 3),
        legend.text=element_text(size=0),
        legend.title=element_text(size=0),
        axis.line = element_line(size=0.15),
        legend.position = "none")








                                  

              
              
              
              
              
              
######  compute macrophages subclasses in tumor samples         
              
MF<- immune.combined 
MF<-subset(MF, subset = region != "normal tissue adjacent to neoplasm")

M1 <- subset(MF, subset = CD68 > 0.1)
M2 <- subset(M1, subset = CD163 > 0.1)

Myeloid <- subset(MF, subset = CD68 < 0.1)
Myeloid <- subset(Myeloid, subset = CD163 > 0.1)

MF@meta.data$MFclass <-
  ifelse(
    rownames(MF@meta.data) %in% colnames(M1),
    "M1",
    "xxx"
  )


MF@meta.data$MFclass <-
  ifelse(
    rownames(MF@meta.data) %in% colnames(M2),
    "M2",
    MF@meta.data$MFclass
  )

MF@meta.data$MFclass <-
  ifelse(
    rownames(MF@meta.data) %in% colnames(Myeloid),
    "Myeloid",
    MF@meta.data$MFclass
  )


# Set identity classes
Idents(object = MF) <- "MFclass"

MF<- subset(MF, idents = c("M1","M2", "Myeloid", "xxx"))
MFx<-MF
MF<- subset(MF, idents = c("M1","M2", "Myeloid"))


              # Find differentially expressed features between M1 and M2
              monocyte.de.markers <- FindMarkers(MF, ident.1 = "M1", ident.2 = "M2")
              
              # sort by fold change
              de.markers <- as.data.frame(monocyte.de.markers)
              de.markers <- de.markers[order(de.markers$avg_log2FC),] 

              setDT(de.markers, keep.rownames = TRUE)[]
              write.table(de.markers, file ="de.markers_M1_M2.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
              
              monocyte.de.markers <- FindMarkers(MF, ident.1 = "M1", ident.2 = "Myeloid")
              # sort by fold change
              de.markers <- as.data.frame(monocyte.de.markers)
              de.markers <- de.markers[order(de.markers$avg_log2FC),] 

              setDT(de.markers, keep.rownames = TRUE)[]
              write.table(de.markers, file ="de.markers_M1_Myeloid.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
              
              
              monocyte.de.markers <- FindMarkers(MF, ident.1 = "Myeloid", ident.2 = "M2")
              # sort by fold change
              de.markers <- as.data.frame(monocyte.de.markers)
              de.markers <- de.markers[order(de.markers$avg_log2FC),] 

              setDT(de.markers, keep.rownames = TRUE)[]
              write.table(de.markers, file ="de.markers_Myeloid_M2.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)


              
# plot heatmaps           
topDifferent <- c("CD68", "CD163", "C1QB" ,"C1QA" ,"C1QC","APOE" ,"APOC1", "CD14" ,  "CCL18" , "HLA.DRB1", "RNASE1" , "HLA.DRA",
                  "TFF3" ,"KRT19"  ,"KRT18"  ,"KRT8" ,"DCN" ,"MGP" ,"AGR2" ,"FXYD3" ,"PHGR1"   ,"LUM"  ,
                  "S100A8"   ,"S100A9"  ,"IGHA2"   ,"MS4A7"    ,"SELENOP"     ,"PDK4"    ,"JCHAIN"    ,"RGS2"    ,"GLUL"  )

DefaultAssay(MF) <- "integrated"
DefaultAssay(MFx) <- "integrated"
# plot only tumor tissue
MF@active.ident <- factor(x = MF@active.ident, levels = c( "M1", "Myeloid", "M2"))
DoHeatmap(MF, features = topDifferent) + NoLegend()
# plot with normal tissue
MFx@active.ident <- factor(x = MFx@active.ident, levels = c( "M1", "Myeloid", "M2", "xxx"))
DoHeatmap(MFx, features = topDifferent) + NoLegend()

# plot with normal tissue
genes <- c("C1QB" ,"C1QA" ,"C1QC","APOE" )
MFx@active.ident <- factor(x = MFx@active.ident, levels = c( "M1", "Myeloid", "M2", "xxx"))
DoHeatmap(MFx, features = genes) + NoLegend()


