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

####
####    stacked violin plots function:   https://rpubs.com/DarrenVan/628853
####

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0.00, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.title.x = element_text(size = 0.4, angle = 40), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(0.6), angle = 0), 
          axis.text.y = element_text(size = rel(0.5), angle = 0), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1), 
          axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

####
####    END of stacked violin plots function
####

setwd("/Volumes/Samsung_X5/Single Cell/E-MTAB-6653/")


# Load the split seurat object into the environment
immune.combined <- readRDS("immune.combined Lung.rds")
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

unique(immune.combined@meta.data$region)





######  compute macrophages subclasses in tumor samples   


MF<- immune.combined
MF<-subset(MF, subset = sample != "normal tissue adjacent to tumour")

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




# plot with normal tissue
DefaultAssay(MFx) <- "integrated"
genes <- c("C1QB" ,"C1QA" ,"C1QC","APOE" )
MFx@active.ident <- factor(x = MFx@active.ident, levels = c( "M1", "Myeloid", "M2", "xxx"))
DoHeatmap(MFx, features = genes) + NoLegend()

DefaultAssay(MF) <- "RNA"

features<- c("CD68", "CD163",  "C1QA", "C1QB", "C1QC", "APOE")
StackedVlnPlot(obj = MF, features = features)


