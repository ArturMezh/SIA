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
library(magrittr)

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


setwd("/Volumes/Samsung_X5/Single Cell/ATLAS/")


# Load the split seurat object into the environment
immune.combined <- readRDS("ATLAS_immune.combined.rds")
DefaultAssay(immune.combined) <- "RNA"

Idents(object = immune.combined) <- "orig.ident"



      

#####
#####                 percentage of cells
#####

          MF<- immune.combined
          M1 <- subset(MF, subset = CD68 > 0.1)
          M2 <- subset(M1, subset = CD163 > 0.1)
          
          Myeloid <- subset(MF, subset = CD68 < 0.1)
          Myeloid <- subset(Myeloid, subset = CD163 > 0.1)
          
          MF@meta.data$MFclass <-
            ifelse(
              rownames(MF@meta.data) %in% colnames(M1),
              "MF",
              "other cells"
            )
          
          
          MF@meta.data$MFclass <-
            ifelse(
              rownames(MF@meta.data) %in% colnames(M2),
              "MF",
              MF@meta.data$MFclass
            )
          
          MF@meta.data$MFclass <-
            ifelse(
              rownames(MF@meta.data) %in% colnames(Myeloid),
              "MF",
              MF@meta.data$MFclass
            )
          
          # Set identity classes
          Idents(object = MF) <- "MFclass"
          
          MF<- subset(MF, idents = c("MF", "other cells"))
          
                    ## extract meta data
                    md <- MF@meta.data %>% as.data.table
                    ## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
                    md <- md[, .N, by = c("sample", "MFclass")]
                    write.table(md, file ="MF fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)

                    
                  ####   for C1Qx-positive cell fraction
                    
                    FeatureScatter(immune.combined, feature1 = "CD68", feature2 = "C1QA")
                                        
                                        MF<- immune.combined
                                        head(MF)
                                        #MF<-subset(MF, subset = seq_folder == "Rectum")
                                        
                                        MFs <- subset(MF, subset = CD68 > 0.1 |CD163 > 0.1)
                                        C1QA <- subset(MF, subset = C1QA > 0.1)
                                        C1QB <- subset(MF, subset = C1QB > 0.1)
                                        C1QC <- subset(MF, subset = C1QC > 0.1)
                                        APOE <- subset(MF, subset = APOE > 0.1)
                                        
                                                              MF@meta.data$MFs <-
                                                                ifelse(
                                                                  rownames(MF@meta.data) %in% colnames(MFs),
                                                                  "MFs",
                                                                  "other cells"
                                                                )
                                                              
                                                              Idents(object = MF) <- "MFs"
                                                              MF<- subset(MF, idents = c("MFs", "other cells"))
                                                              md <- MF@meta.data %>% as.data.table
                                                              md <- md[, .N, by = c("sample", "MFs")]
                                                              write.table(md, file ="MFs fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                                                              
                                                              
                                                              
                                                              MF@meta.data$C1QA <-
                                                                ifelse(
                                                                  rownames(MF@meta.data) %in% colnames(C1QA),
                                                                  "C1QA",
                                                                  "other cells"
                                                                )
                                                              
                                                              Idents(object = MF) <- "C1QA"
                                                              MF<- subset(MF, idents = c("C1QA", "other cells"))
                                                              md <- MF@meta.data %>% as.data.table
                                                              md <- md[, .N, by = c("sample", "C1QA")]
                                                              write.table(md, file ="C1QA fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                                                              
                                                              
                                                              MF@meta.data$C1QB <-
                                                                ifelse(
                                                                  rownames(MF@meta.data) %in% colnames(C1QB),
                                                                  "C1QB",
                                                                  "other cells"
                                                                )
                                                              
                                                              Idents(object = MF) <- "C1QB"
                                                              MF<- subset(MF, idents = c("C1QB", "other cells"))
                                                              md <- MF@meta.data %>% as.data.table
                                                              md <- md[, .N, by = c("sample", "C1QB")]
                                                              write.table(md, file ="C1QB fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                                                              
                                                              
                                                              MF@meta.data$C1QC <-
                                                                ifelse(
                                                                  rownames(MF@meta.data) %in% colnames(C1QC),
                                                                  "C1QC",
                                                                  "other cells"
                                                                )
                                                              
                                                              Idents(object = MF) <- "C1QC"
                                                              MF<- subset(MF, idents = c("C1QC", "other cells"))
                                                              md <- MF@meta.data %>% as.data.table
                                                              md <- md[, .N, by = c("sample", "C1QC")]
                                                              write.table(md, file ="C1QC fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                                                              
                                                              
                                                              MF@meta.data$APOE <-
                                                                ifelse(
                                                                  rownames(MF@meta.data) %in% colnames(APOE),
                                                                  "APOE",
                                                                  "other cells"
                                                                )
                                                              
                                                              Idents(object = MF) <- "APOE"
                                                              MF<- subset(MF, idents = c("APOE", "other cells"))
                                                              md <- MF@meta.data %>% as.data.table
                                                              md <- md[, .N, by = c("sample", "APOE")]
                                                              write.table(md, file ="APOE fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                                        
                                        MF_C1QA <- subset(MFs, subset = C1QA > 0.1)
                                        MF_C1QB <- subset(MFs, subset = C1QB > 0.1)
                                        MF_C1QC <- subset(MFs, subset = C1QC > 0.1)
                                        MF_APOE <- subset(MFs, subset = APOE > 0.1)
                                                
                                                MF@meta.data$MF_C1QA <-
                                                  ifelse(
                                                    rownames(MF@meta.data) %in% colnames(MF_C1QA),
                                                    "MF_C1QA",
                                                    "other cells"
                                                  )
                                                
                                                Idents(object = MF) <- "MF_C1QA"
                                                MF<- subset(MF, idents = c("MF_C1QA", "other cells"))
                                                md <- MF@meta.data %>% as.data.table
                                                md <- md[, .N, by = c("sample", "MF_C1QA")]
                                                write.table(md, file ="MF_C1QA fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                                                
                                                MF@meta.data$MF_C1QB <-
                                                  ifelse(
                                                    rownames(MF@meta.data) %in% colnames(MF_C1QB),
                                                    "MF_C1QB",
                                                    "other cells"
                                                  )
                                                
                                                Idents(object = MF) <- "MF_C1QB"
                                                MF<- subset(MF, idents = c("MF_C1QB", "other cells"))
                                                md <- MF@meta.data %>% as.data.table
                                                md <- md[, .N, by = c("sample", "MF_C1QB")]
                                                write.table(md, file ="MF_C1QB fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                                                
                                                MF@meta.data$MF_C1QC <-
                                                  ifelse(
                                                    rownames(MF@meta.data) %in% colnames(MF_C1QC),
                                                    "MF_C1QC",
                                                    "other cells"
                                                  )
                                                
                                                Idents(object = MF) <- "MF_C1QC"
                                                MF<- subset(MF, idents = c("MF_C1QC", "other cells"))
                                                md <- MF@meta.data %>% as.data.table
                                                md <- md[, .N, by = c("sample", "MF_C1QC")]
                                                write.table(md, file ="MF_C1QC fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)

                                                
                                                MF@meta.data$MF_APOE <-
                                                  ifelse(
                                                    rownames(MF@meta.data) %in% colnames(MF_APOE),
                                                    "MF_APOE",
                                                    "other cells"
                                                  )
                                                
                                                Idents(object = MF) <- "MF_APOE"
                                                MF<- subset(MF, idents = c("MF_APOE", "other cells"))
                                                md <- MF@meta.data %>% as.data.table
                                                md <- md[, .N, by = c("sample", "MF_APOE")]
                                                write.table(md, file ="MF_APOE fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                                                

                                                
                                                #   for counting fraction of M2-like MFs which are C1Q POSITIVE
                                                
                                                MF<- immune.combined

                                                MF2 <- subset(MF, subset = CD68 > 0.1 & CD163 > 0.1)
                                                            
                                                            MF2_C1QA <- subset(MF2, subset = C1QA > 0.1)
                                                            MF2_C1QB <- subset(MF2, subset = C1QB > 0.1)
                                                            MF2_C1QC <- subset(MF2, subset = C1QC > 0.1)
                                                            MF2_APOE <- subset(MF2, subset = APOE > 0.1)
                                                            
                                                            MF2@meta.data$MF2_C1QA <-
                                                              ifelse(
                                                                rownames(MF2@meta.data) %in% colnames(MF2_C1QA),
                                                                "MF2_C1QA",
                                                                "other cells"
                                                              )
                                                            
                                                            Idents(object = MF2) <- "MF2_C1QA"
                                                            MF2<- subset(MF2, idents = c("MF2_C1QA", "other cells"))
                                                            md <- MF2@meta.data %>% as.data.table
                                                            md <- md[, .N, by = c("sample", "MF2_C1QA")]
                                                            write.table(md, file ="MF2_C1QA fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                                                            
                                                            MF2@meta.data$MF2_C1QB <-
                                                              ifelse(
                                                                rownames(MF2@meta.data) %in% colnames(MF2_C1QB),
                                                                "MF2_C1QB",
                                                                "other cells"
                                                              )
                                                            
                                                            Idents(object = MF2) <- "MF2_C1QB"
                                                            MF2<- subset(MF2, idents = c("MF2_C1QB", "other cells"))
                                                            md <- MF2@meta.data %>% as.data.table
                                                            md <- md[, .N, by = c("sample", "MF2_C1QB")]
                                                            write.table(md, file ="MF2_C1QB fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                                                            
                                                            MF2@meta.data$MF2_C1QC <-
                                                              ifelse(
                                                                rownames(MF2@meta.data) %in% colnames(MF2_C1QC),
                                                                "MF2_C1QC",
                                                                "other cells"
                                                              )
                                                            
                                                            Idents(object = MF2) <- "MF2_C1QC"
                                                            MF2<- subset(MF2, idents = c("MF2_C1QC", "other cells"))
                                                            md <- MF2@meta.data %>% as.data.table
                                                            md <- md[, .N, by = c("sample", "MF2_C1QC")]
                                                            write.table(md, file ="MF2_C1QC fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                                                            
                                                            MF2@meta.data$MF2_APOE <-
                                                              ifelse(
                                                                rownames(MF2@meta.data) %in% colnames(MF2_APOE),
                                                                "MF2_APOE",
                                                                "other cells"
                                                              )
                                                            
                                                            Idents(object = MF2) <- "MF2_APOE"
                                                            MF2<- subset(MF2, idents = c("MF2_APOE", "other cells"))
                                                            md <- MF2@meta.data %>% as.data.table
                                                            md <- md[, .N, by = c("sample", "MF2_APOE")]
                                                            write.table(md, file ="MF2_APOE fraction.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)


#####
#####                 percentage of cells END
#####

    

###

MF<- immune.combined
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
MF<- subset(MF, idents = c("M1","M2", "Myeloid"))
# list options for groups to perform differential expression on
cell.num <- table(MF@active.ident)

                    ## extract meta data
                    md <- MF@meta.data %>% as.data.table
                    
                    ## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
                    md <- md[, .N, by = c("sample", "MFclass")]
                    write.table(md, file ="md.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
                    
      
      

# Find differentially expressed features between M1 and M2
monocyte.de.markers <- FindMarkers(MF, ident.1 = "M1", ident.2 = "M2")
monocyte.de.markers <- FindMarkers(MF, ident.1 = "Myeloid", ident.2 = "M2")
# sort by fold change
de.markers <- as.data.frame(monocyte.de.markers)
de.markers <- de.markers[order(de.markers$avg_log2FC),] 

setDT(de.markers, keep.rownames = TRUE)[]
write.table(de.markers, file ="de.markers.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)






###
Idents(object = MF) <- "MFclass"
features<- c("CD68", "CD163",  "C1QA", "C1QB", "C1QC", "APOE")
DefaultAssay(MF) <- "RNA"



plots <- VlnPlot(object = MF, features =c("CD68", "CD163",  "C1QA", "C1QB", "C1QC", "APOE"), group.by = "MFclass", 
                 split.by = "sample",pt.size = 0, ncol = 1,  same.y.lims=T, combine = FALSE)

library(ggforce)
for(i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_bar(position="dodge", stat = "summary")  + theme(legend.position = 'none')
}
CombinePlots(plots)
CombinePlots(plots, ncol = 1)
