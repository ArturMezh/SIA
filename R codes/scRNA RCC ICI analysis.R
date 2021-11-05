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


setwd("/Volumes/Samsung_X5/Single Cell/RCC_ICI/")

# Load the split seurat object into the environment
immune.combined <- readRDS("immune.combined RCC_ICI.rds")
DefaultAssay(immune.combined) <- "RNA"
## 

metadata <- immune.combined@meta.data
          immune.combined@meta.data
          immune.combined.2 <- immune.combined
          
          immune.combined.3 <- subset(immune.combined.2, subset = CD8A > 0.0)
          CD8A <- table(immune.combined.3@meta.data$donor_id)
          CD8A <- as.data.frame(CD8A)

          immune.combined.3 <- subset(immune.combined.2, subset = CD68 > 0.0 &C1QA > 0.0)
          C1QA <- table(immune.combined.3@meta.data$donor_id)
          C1QA <- as.data.frame(C1QA)
          
          immune.combined.3 <- subset(immune.combined.2, subset = CD68 > 0.0 &C1QB > 0.0)
          C1QB <- table(immune.combined.3@meta.data$donor_id)
          C1QB <- as.data.frame(C1QB)
          
          immune.combined.3 <- subset(immune.combined.2, subset = CD68 > 0.0 &C1QC > 0.0)
          C1QC <- table(immune.combined.3@meta.data$donor_id)
          C1QC <- as.data.frame(C1QC)
          
          immune.combined.3 <- immune.combined.2
          total <- table(immune.combined.3@meta.data$donor_id)
          total <- as.data.frame(total)

          immune.combined.3 <- subset(immune.combined.2, subset = CD68 > 0.0 &  CD163 > 0.0)
          M2 <- table(immune.combined.3@meta.data$donor_id)
          M2 <- as.data.frame(M2)
          
          
          
          joint <- full_join(CD8A,M2, by="Var1")
          joint <- full_join(joint,C1QA, by="Var1")
          joint <- full_join(joint,C1QB, by="Var1")
          joint <- full_join(joint,C1QC, by="Var1")
          joint <- full_join(joint,total, by="Var1")
          
          colnames(joint) <- c('donor_id','CD8A','M2', 'C1QA', 'C1QB', 'C1QC', 'total')
          
          joint[is.na(joint)] <- 0

          joint$SIAold <- joint$CD8A/(joint$M2+joint$CD8A)
          joint$SIA_C1QA <- joint$CD8A/(joint$C1QA+joint$CD8A)
          joint$SIA_C1QB <- joint$CD8A/(joint$C1QB+joint$CD8A)
          joint$SIA_C1QC <- joint$CD8A/(joint$C1QC+joint$CD8A)
          
          short_meta <- metadata[,c(9,19)]
          unique(short_meta$donor_id)
          short_meta <- short_meta[!duplicated(short_meta[,c(1)]), ]
          short_meta <- unique(short_meta[,c(1,2)])
          
          SIA_groups <- left_join(short_meta, joint, by="donor_id")
          SIA_groups <- na.omit(SIA_groups)
          
          colnames(SIA_groups)
          XXXX<-SIA_groups[,c(2,9:12)]
          XXXX<-XXXX[which(XXXX$ICB_Response != "ICB_NE" & XXXX$ICB_Response != "NoICB"),]
          XXXX = reshape2::melt(XXXX, id = c("ICB_Response"))
          XXXX
          library(dplyr)

          k <- ggplot(XXXX, aes(variable, value, fill=`ICB_Response`))
          
          
          k+ theme_bw()+ theme(panel.border = element_blank())+ theme(axis.line = element_line(colour = "black"))+
            geom_point(aes(colour = `ICB_Response`, shape= `ICB_Response`), position=position_dodge(0.3), size=5.5)+
            theme(text = element_text(size=10),
                  axis.text.x = element_text(angle=0, vjust=1,size=10)) +
      coord_cartesian(ylim = c(0.6, 1)) 
      