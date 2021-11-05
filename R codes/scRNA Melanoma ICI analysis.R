########################
# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)


setwd("/Volumes/Samsung_X5/Single Cell/MEL_ICI")

# Load the split seurat object into the environment
immune.combined <- readRDS("MEL_ICI_immune.combined.rds")
DefaultAssay(immune.combined) <- "RNA"
metadata <- immune.combined@meta.data
unique(metadata$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)
metadata$sample <-  metadata$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.
metadata$PrePost <-  gsub('_.*', '', metadata$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)
metadata$patient <- gsub('.*_P', 'P', metadata$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)
immune.combined@meta.data <- metadata


###
 metadata <- immune.combined@meta.data
          immune.combined.2 <- immune.combined
          
          immune.combined.3 <- subset(immune.combined.2, subset = CD8A > 0.0)
          CD8A <- table(immune.combined.3@meta.data$sample)
          CD8A <- as.data.frame(CD8A)

          immune.combined.3 <- subset(immune.combined.2, subset = CD68 > 0.0 &C1QA > 0.0)
          C1QA <- table(immune.combined.3@meta.data$sample)
          C1QA <- as.data.frame(C1QA)
          
          immune.combined.3 <- subset(immune.combined.2, subset = CD68 > 0.0 &C1QB > 0.0)
          C1QB <- table(immune.combined.3@meta.data$sample)
          C1QB <- as.data.frame(C1QB)
          
          immune.combined.3 <- subset(immune.combined.2, subset = CD68 > 0.0 &C1QC > 0.0)
          C1QC <- table(immune.combined.3@meta.data$sample)
          C1QC <- as.data.frame(C1QC)
          
          immune.combined.3 <- immune.combined.2
          total <- table(immune.combined.3@meta.data$sample)
          total <- as.data.frame(total)
          
          immune.combined.3 <- subset(immune.combined.2, subset = CD68 > 0.0 &  CD163 > 0.0)
          M2 <- table(immune.combined.3@meta.data$sample)
          M2 <- as.data.frame(M2)

          
          joint <- full_join(CD8A,M2, by="Var1")
          joint <- full_join(joint,C1QA, by="Var1")
          joint <- full_join(joint,C1QB, by="Var1")
          joint <- full_join(joint,C1QC, by="Var1")
          joint <- full_join(joint,total, by="Var1")
          
          colnames(joint) <- c('sample','CD8A','M2', 'C1QA', 'C1QB', 'C1QC','total')
          
          joint[is.na(joint)] <- 0
          
          joint$SIAold <- joint$CD8A/(joint$M2+joint$CD8A)
          joint$SIA_C1QA <- joint$CD8A/(joint$C1QA+joint$CD8A)
          joint$SIA_C1QB <- joint$CD8A/(joint$C1QB+joint$CD8A)
          joint$SIA_C1QC <- joint$CD8A/(joint$C1QC+joint$CD8A)
          

          short_meta <- metadata[,c(26,13)]
          short_meta <- short_meta[!duplicated(short_meta[,c(1)]), ]
          short_meta <- unique(short_meta[,c(1,2)])
          
          SIA_groups <- left_join(short_meta, joint, by="sample")
          SIA_groups <- na.omit(SIA_groups)

          wilcox.test(SIAold ~ characteristics..response, data=SIA_groups, paired = FALSE)
          wilcox.test(SIA_C1QA ~ characteristics..response, data=SIA_groups, paired = FALSE)
          wilcox.test(SIA_C1QB ~ characteristics..response, data=SIA_groups, paired = FALSE)
          wilcox.test(SIA_C1QC ~ characteristics..response, data=SIA_groups, paired = FALSE)

          XXXX<-SIA_groups[,c(2,9:12)]
          XXXX = reshape2::melt(XXXX, id = c("characteristics..response"))

          
          
          
          k <- ggplot(XXXX, aes(variable, value, fill=`characteristics..response`))
          
          k+ theme_bw()+ theme(panel.border = element_blank())+ theme(axis.line = element_line(colour = "black"))+
            geom_boxplot(outlier.shape=NA, stat = "boxplot",   alpha=0.5, coef=0)+
            theme(text = element_text(size=10),
                  axis.text.x = element_text(angle=0, vjust=1,size=10)) +
            coord_cartesian(ylim = c(0.8, 1)) 
          
          
##### AUC          

          
          SIA_groups$characteristics..response <- ifelse((SIA_groups$characteristics..response == "Responder"), 1, 0)
          
          
          library(precrec)

          precrec_obj <- evalmod(scores = SIA_groups$SIAold, labels = SIA_groups$characteristics..response)
          autoplot(precrec_obj)
          precrec_obj <- evalmod(scores = SIA_groups$SIA_C1QA, labels = SIA_groups$characteristics..response)
          autoplot(precrec_obj)
          precrec_obj <- evalmod(scores = SIA_groups$SIA_C1QB, labels = SIA_groups$characteristics..response)
          autoplot(precrec_obj)
          precrec_obj <- evalmod(scores = SIA_groups$SIA_C1QC, labels = SIA_groups$characteristics..response)
          autoplot(precrec_obj)
          
          