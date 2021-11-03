########################
########################
require(gdata)
library(plyr)
require(xlsx)

                                          ##################################################
                                          ############ ***     TIL panel    *** ############
                                          ##################################################



##################################################
######### *** reading original files *** #########

setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub/datasets/validation cohort/")

##
### Endometrium
##
mergedTIL1 <-data.table::fread("Endometrium_TIL_1_merged_cell_seg_data.txt")
mergedTIL2 <-data.table::fread("Endometrium_TIL_2_merged_cell_seg_data.txt")

TIL1 <- mergedTIL1[,c(2,3, 8,9,   72, 77,  82, 41,  97)]  
TIL2 <- mergedTIL2[,c(2,3, 8,9,   72, 77,  82, 41,  97)]  

TIL <- rbind(TIL1, TIL2)
TIL <- na.omit(TIL)
names(TIL)[c(5:9)]<-c("CD4","CD20","CD8","FoxP3","CK")
TIL$`Sample Name` <- gsub("]_.*", "]", TIL$`Sample Name`)
TIL$`Sample Name` <- gsub('Core...', '[', TIL$`Sample Name`)
TIL$CD4 <- as.numeric(gsub(",", ".", TIL$CD4, fixed = T))
TIL$CD20 <- as.numeric(gsub(",", ".", TIL$CD20, fixed = T))
TIL$CD8 <- as.numeric(gsub(",", ".", TIL$CD8, fixed = T))
TIL$FoxP3 <- as.numeric(gsub(",", ".", TIL$FoxP3, fixed = T))
TIL$CK <- as.numeric(gsub(",", ".", TIL$CK, fixed = T))
TIL$`Cell X Position` <- as.numeric(gsub(",", ".", TIL$`Cell X Position`, fixed = T))
TIL$`Cell Y Position` <- as.numeric(gsub(",", ".", TIL$`Cell Y Position`, fixed = T))
TIL <- TIL[ which(TIL$`Tissue Category` != 'Blanc'), ]                     
TIL <- na.omit(TIL)
Endometrium <-TIL

##
### END Endometrium
##




##
### Esoph
##

mergedTIL <-data.table::fread("Esoph_TIL_merged_cell_seg_data.txt")

TIL <- mergedTIL[,c(2,3, 8,9,    72, 77, 82,  41,  97)]  
TIL <- na.omit(TIL)
names(TIL)[c(5:9)]<-c("CD4","CD20","CD8","FoxP3","CK")
TIL$`Sample Name` <- gsub("]_.*", "]", TIL$`Sample Name`)
TIL$`Sample Name` <- gsub('Core...', '[', TIL$`Sample Name`)
TIL$CD4 <- as.numeric(gsub(",", ".", TIL$CD4, fixed = T))
TIL$CD20 <- as.numeric(gsub(",", ".", TIL$CD20, fixed = T))
TIL$CD8 <- as.numeric(gsub(",", ".", TIL$CD8, fixed = T))
TIL$FoxP3 <- as.numeric(gsub(",", ".", TIL$FoxP3, fixed = T))
TIL$CK <- as.numeric(gsub(",", ".", TIL$CK, fixed = T))
TIL$`Cell X Position` <- as.numeric(gsub(",", ".", TIL$`Cell X Position`, fixed = T))
TIL$`Cell Y Position` <- as.numeric(gsub(",", ".", TIL$`Cell Y Position`, fixed = T))
TIL <- TIL[ which(TIL$`Tissue Category` != 'Blanc'), ]                     
Esoph <- na.omit(TIL)

##
### END Esoph
##





##
###  Lung
##

mergedTIL <-data.table::fread("Lung_TIL_merged_cell_seg_data.txt")

TIL <- mergedTIL[,c(2,3, 8,9,   72, 77, 82,  41,  97)]  
TIL <- na.omit(TIL)
names(TIL)[c(5:9)]<-c("CD4","CD20","CD8","FoxP3","CK")
TIL$`Sample Name` <- gsub("]_.*", "]", TIL$`Sample Name`)
TIL$`Sample Name` <- gsub('Core...', '[', TIL$`Sample Name`)
TIL$CD4 <- as.numeric(gsub(",", ".", TIL$CD4, fixed = T))
TIL$CD20 <- as.numeric(gsub(",", ".", TIL$CD20, fixed = T))
TIL$CD8 <- as.numeric(gsub(",", ".", TIL$CD8, fixed = T))
TIL$FoxP3 <- as.numeric(gsub(",", ".", TIL$FoxP3, fixed = T))
TIL$CK <- as.numeric(gsub(",", ".", TIL$CK, fixed = T))
TIL$`Cell X Position` <- as.numeric(gsub(",", ".", TIL$`Cell X Position`, fixed = T))
TIL$`Cell Y Position` <- as.numeric(gsub(",", ".", TIL$`Cell Y Position`, fixed = T))
TIL <- TIL[ which(TIL$`Tissue Category` != 'Blanc'), ]                     
Lung <- na.omit(TIL)
##
### END Lung
##





##
###  MEL
##

mergedTIL <-data.table::fread("Melanoma_TIL_merged_cell_seg_data.txt")

TIL <- mergedTIL[,c(2,3, 8,9,   67, 72, 77,  41,  87)]  
TIL <- na.omit(TIL)
names(TIL)[c(5:9)]<-c("CD4","CD20","CD8","FoxP3","CK")
TIL$`Sample Name` <- gsub("]_.*", "]", TIL$`Sample Name`)
TIL$`Sample Name` <- gsub('Core...', '[', TIL$`Sample Name`)
TIL$CD4 <- as.numeric(gsub(",", ".", TIL$CD4, fixed = T))
TIL$CD20 <- as.numeric(gsub(",", ".", TIL$CD20, fixed = T))
TIL$CD8 <- as.numeric(gsub(",", ".", TIL$CD8, fixed = T))
TIL$FoxP3 <- as.numeric(gsub(",", ".", TIL$FoxP3, fixed = T))
TIL$CK <- as.numeric(gsub(",", ".", TIL$CK, fixed = T))
TIL$`Cell X Position` <- as.numeric(gsub(",", ".", TIL$`Cell X Position`, fixed = T))
TIL$`Cell Y Position` <- as.numeric(gsub(",", ".", TIL$`Cell Y Position`, fixed = T))
TIL <- TIL[ which(TIL$`Tissue Category` != 'Blanc'), ]                     
MEL <- na.omit(TIL)

##
### END MEL
##


##
###  Ovca_Lund
##

mergedTIL <-data.table::fread("OvcaLund_TIL_merged_cell_seg_data.txt")

TIL <- mergedTIL[,c(2,3, 8,9,   67, 72, 77,  41,  87)]  
TIL <- na.omit(TIL)
names(TIL)[c(5:9)]<-c("CD4","CD20","CD8","FoxP3","CK")
TIL$`Sample Name` <- gsub("]_.*", "]", TIL$`Sample Name`)
TIL$`Sample Name` <- gsub('Core...', '[', TIL$`Sample Name`)
TIL$CD4 <- as.numeric(gsub(",", ".", TIL$CD4, fixed = T))
TIL$CD20 <- as.numeric(gsub(",", ".", TIL$CD20, fixed = T))
TIL$CD8 <- as.numeric(gsub(",", ".", TIL$CD8, fixed = T))
TIL$FoxP3 <- as.numeric(gsub(",", ".", TIL$FoxP3, fixed = T))
TIL$CK <- as.numeric(gsub(",", ".", TIL$CK, fixed = T))
TIL$`Cell X Position` <- as.numeric(gsub(",", ".", TIL$`Cell X Position`, fixed = T))
TIL$`Cell Y Position` <- as.numeric(gsub(",", ".", TIL$`Cell Y Position`, fixed = T))
TIL <- TIL[ which(TIL$`Tissue Category` != 'Blanc'), ]                     
Ovca_Lund <- na.omit(TIL)

##
### END Ovca_Lund
##


##
###  UBladder
##
mergedTIL <-data.table::fread("UroBladder_TIL_merged_cell_seg_data.txt")

TIL <- mergedTIL[,c(2,3, 8,9,   72, 77, 82, 41,  97)]  
TIL <- na.omit(TIL)
names(TIL)[c(5:9)]<-c("CD4","CD20","CD8","FoxP3","CK")
TIL$`Sample Name` <- gsub("]_.*", "]", TIL$`Sample Name`)
TIL$`Sample Name` <- gsub('Core...', '[', TIL$`Sample Name`)
TIL$CD4 <- as.numeric(gsub(",", ".", TIL$CD4, fixed = T))
TIL$CD20 <- as.numeric(gsub(",", ".", TIL$CD20, fixed = T))
TIL$CD8 <- as.numeric(gsub(",", ".", TIL$CD8, fixed = T))
TIL$FoxP3 <- as.numeric(gsub(",", ".", TIL$FoxP3, fixed = T))
TIL$CK <- as.numeric(gsub(",", ".", TIL$CK, fixed = T))
TIL$`Cell X Position` <- as.numeric(gsub(",", ".", TIL$`Cell X Position`, fixed = T))
TIL$`Cell Y Position` <- as.numeric(gsub(",", ".", TIL$`Cell Y Position`, fixed = T))
TIL <- TIL[ which(TIL$`Tissue Category` != 'Blanc'), ]                     
UBladder <- na.omit(TIL)

##
### END UBladder
##



#######
######
#######


                            ### CK threshold algorithm
                            
                            CK_calculator <- function(fdata){
                              require(plyr)
                              library(autothresholdr)
                              
                              otsu<- as.integer(fdata$CK*1000)
                              o_CK_fdata <- auto_thresh(otsu, "Otsu")
                              o_CK_fdata <- as.numeric(o_CK_fdata)/1000
                              First_distribution <- fdata[ which(fdata$CK  < o_CK_fdata), ]
                              Second_distribution <- fdata[ which(fdata$CK  > o_CK_fdata), ]
                              peak1 <- density(First_distribution$CK)
                              peak2 <- density(Second_distribution$CK)
                              # get highest Peak
                              peak1 <-peak1$x[which.max(peak1$y)] #Gives you the first highest Peak
                              peak2 <-peak2$x[which.max(peak2$y)] #Gives you the second highest Peak
                              dY <- density(fdata$CK)$y
                              dX <- density(fdata$CK)$x
                              MinYDensity<- min(dY[dX  > peak1 & dX < peak2])
                              minimum_CK <-dX[which(dY == MinYDensity)]
                              
                              return(minimum_CK)
                            }
                            
                            
                                    minimum_CK <- CK_calculator(Endometrium)
                                    Endometrium$CK <- ifelse((Endometrium$CK>minimum_CK), 1, 0)
                                    minimum_CK <- CK_calculator(Esoph)
                                    Esoph$CK <- ifelse((Esoph$CK>minimum_CK), 1, 0)
                                    minimum_CK <- CK_calculator(Lung)
                                    Lung$CK <- ifelse((Lung$CK>minimum_CK), 1, 0)
                                    minimum_CK <- CK_calculator(MEL)
                                    MEL$CK <- ifelse((MEL$CK>minimum_CK), 1, 0)
                                    minimum_CK <- CK_calculator(Ovca_Lund)
                                    Ovca_Lund$CK <- ifelse((Ovca_Lund$CK>minimum_CK), 1, 0)
                                    minimum_CK <- CK_calculator(UBladder)
                                    UBladder$CK <- ifelse((UBladder$CK>minimum_CK), 1, 0)

    
    Lung_match <- as.data.frame(read.xlsx2("matching_reference_Lung_TIL.xlsx", sheetIndex = 1))
    Endometrium_match <- as.data.frame(read.xlsx2("matching_reference_Endometrium_TIL.xlsx", sheetIndex = 1))
    Esoph_match <- as.data.frame(read.xlsx2("matching_reference_Esoph_TIL.xlsx", sheetIndex = 1))
    Ovca_Lund_match <- as.data.frame(read.xlsx2("matching_reference_OvarianLund_TIL.xlsx", sheetIndex = 1))
    UBladder_match <- as.data.frame(read.xlsx2("matching_reference_UroBladder_TIL.xlsx", sheetIndex = 1))
    MEL_match <- as.data.frame(read.xlsx2("matching_reference_Melanoma_TIL.xlsx", sheetIndex = 1))
                                    
                                  Lung$ID <- Lung_match$ID[match(unlist(Lung$`Sample Name`), Lung_match$Sample.Name)]
                                  Lung$sample_region <- Lung_match$sample_region[match(unlist(Lung$`Sample Name`), Lung_match$Sample.Name)]
                                  Lung$`Cell X Position` <- Lung$`Cell X Position`*0.497
                                  Lung$`Cell Y Position` <- Lung$`Cell Y Position`*0.497
                                  Lung <- Lung[ which(Lung$`Tissue Category` != 'Blanc'), ]
                                  Lung <- na.omit(Lung)           

                                              TIL <- NULL
                                              TIL <- Lung
                                              CD4cutoff <-3.096276594
                                              CD8cutoff <-1.7325162674856
                                              CD20cutoff <-1.35090303491679
                                              FoxP3cutoff <-2.13288056121903
                                              
                                              TIL$CD4 <- ifelse((TIL$CD4>CD4cutoff), 1, 0)
                                              TIL$CD8 <- ifelse((TIL$CD8>CD8cutoff), 1, 0)
                                              TIL$CD20 <- ifelse((TIL$CD20>CD20cutoff), 1, 0)
                                              TIL$FoxP3 <- ifelse((TIL$FoxP3>FoxP3cutoff), 1, 0)
                                              Lung <- TIL
                                  
                                  Endometrium$ID <- Endometrium_match$ID[match(unlist(Endometrium$`Sample Name`), Endometrium_match$Sample.Name)]
                                  Endometrium$`Cell X Position` <- Endometrium$`Cell X Position`*0.497
                                  Endometrium$`Cell Y Position` <- Endometrium$`Cell Y Position`*0.497
                                  Endometrium <- Endometrium[ which(Endometrium$`Tissue Category` != 'Blanc'), ]
                                  Endometrium <- na.omit(Endometrium)  
                                  
                                  TIL <- NULL
                                  TIL <- Endometrium
                                  CD4cutoff <-2.173293187044
                                  CD8cutoff <-1.7325162674856
                                  CD20cutoff <-1.35090303491679
                                  FoxP3cutoff <-3.554125281 
                                  
                                  TIL$CD4 <- ifelse((TIL$CD4>CD4cutoff), 1, 0)
                                  TIL$CD8 <- ifelse((TIL$CD8>CD8cutoff), 1, 0)
                                  TIL$CD20 <- ifelse((TIL$CD20>CD20cutoff), 1, 0)
                                  TIL$FoxP3 <- ifelse((TIL$FoxP3>FoxP3cutoff), 1, 0)

                                  Endometrium <- TIL
                                  
                                                  Esoph$ID <- Esoph_match$ID[match(unlist(Esoph$`Sample Name`), Esoph_match$Sample.Name)]
                                                  Esoph$sample_region <- Esoph_match$sample_region[match(unlist(Esoph$`Sample Name`), Esoph_match$Sample.Name)]
                                                  Esoph$`Cell X Position` <- Esoph$`Cell X Position`*0.497
                                                  Esoph$`Cell Y Position` <- Esoph$`Cell Y Position`*0.497
                                                  Esoph <- Esoph[ which(Esoph$`Tissue Category` != 'Blanc'), ]
                                                  Esoph <- na.omit(Esoph)  
                                                  
                                                  TIL <- NULL
                                                  TIL <- Esoph
                                                  CD4cutoff <-4.5591717887172
                                                  CD8cutoff <-1.7325162674856
                                                  CD20cutoff <-1.35090303491679
                                                  FoxP3cutoff <-2.51507028060951
                                                  
                                                  TIL$CD4 <- ifelse((TIL$CD4>CD4cutoff), 1, 0)
                                                  TIL$CD8 <- ifelse((TIL$CD8>CD8cutoff), 1, 0)
                                                  TIL$CD20 <- ifelse((TIL$CD20>CD20cutoff), 1, 0)
                                                  TIL$FoxP3 <- ifelse((TIL$FoxP3>FoxP3cutoff), 1, 0)
                                                  
                                                  Esoph <- TIL

                                  
                                  Ovca_Lund$ID <- Ovca_Lund_match$ID[match(unlist(Ovca_Lund$`Sample Name`), Ovca_Lund_match$Sample.Name)]
                                  Ovca_Lund$sample_region <- Ovca_Lund_match$sample_region[match(unlist(Ovca_Lund$`Sample Name`), Ovca_Lund_match$Sample.Name)]
                                  Ovca_Lund <- Ovca_Lund[ which(Ovca_Lund$`Tissue Category` != 'Blanc'), ]
                                  Ovca_Lund <- na.omit(Ovca_Lund)  
                                  
                                  TIL <- NULL
                                  TIL <- Ovca_Lund
                                  CD4cutoff <-2.173293187044
                                  CD8cutoff <-1.7325162674856
                                  CD20cutoff <-1.35090303491679
                                  FoxP3cutoff <-2.13288056121903
                                  
                                  TIL$CD4 <- ifelse((TIL$CD4>CD4cutoff), 1, 0)
                                  TIL$CD8 <- ifelse((TIL$CD8>CD8cutoff), 1, 0)
                                  TIL$CD20 <- ifelse((TIL$CD20>CD20cutoff), 1, 0)
                                  TIL$FoxP3 <- ifelse((TIL$FoxP3>FoxP3cutoff), 1, 0)
                                  
                                  Ovca_Lund <- TIL
                                  

                                                  UBladder$ID <- UBladder_match$PAD_year[match(unlist(UBladder$`Sample Name`), UBladder_match$Sample.Name)]
                                                  UBladder$`Cell X Position` <- UBladder$`Cell X Position`*0.497
                                                  UBladder$`Cell Y Position` <- UBladder$`Cell Y Position`*0.497
                                                  UBladder <- UBladder[ which(UBladder$`Tissue Category` != 'Blanc'), ]
                                                  UBladder <- na.omit(UBladder) 
                                                  
                                                  TIL <- NULL
                                                  TIL <- UBladder
                                                  CD4cutoff <-2.173293187044
                                                  CD8cutoff <-1.7325162674856
                                                  CD20cutoff <-1.35090303491679
                                                  FoxP3cutoff <-2.13288056121903
                                                  
                                                  TIL$CD4 <- ifelse((TIL$CD4>CD4cutoff), 1, 0)
                                                  TIL$CD8 <- ifelse((TIL$CD8>CD8cutoff), 1, 0)
                                                  TIL$CD20 <- ifelse((TIL$CD20>CD20cutoff), 1, 0)
                                                  TIL$FoxP3 <- ifelse((TIL$FoxP3>FoxP3cutoff), 1, 0)
                                                  
                                                  UBladder <- TIL

                                  MEL$ID <- MEL_match$ID[match(unlist(MEL$`Sample Name`), MEL_match$Sample.Name)]
                                  MEL <- MEL[ which(MEL$`Tissue Category` != 'Blanc'), ]
                                  MEL <- na.omit(MEL) 
                                  
                                  TIL <- NULL
                                  TIL <- MEL
                                  CD4cutoff <-1.704258976
                                  CD8cutoff <-1.7325162674856
                                  CD20cutoff <-1.35090303491679
                                  FoxP3cutoff <-2.13288056121903
                                  
                                  TIL$CD4 <- ifelse((TIL$CD4>CD4cutoff), 1, 0)
                                  TIL$CD8 <- ifelse((TIL$CD8>CD8cutoff), 1, 0)
                                  TIL$CD20 <- ifelse((TIL$CD20>CD20cutoff), 1, 0)
                                  TIL$FoxP3 <- ifelse((TIL$FoxP3>FoxP3cutoff), 1, 0)
                                  
                                  MEL <- TIL

                      Lung[,"cohort"]='Lung'
                      Endometrium[,"sample_region"]='T'
                      Endometrium[,"cohort"]='Endometrium'
                      Esoph[,"cohort"]='Esoph'
                      Ovca_Lund[,"cohort"]='Ovca_Lund'
                      UBladder[,"sample_region"]='T'
                      UBladder[,"cohort"]='UBladder'
                      MEL[,"sample_region"]='T'
                      MEL[,"cohort"]='MEL'
   
VALIDATION_TIL <- rbind(Endometrium, Esoph, Lung, MEL, Ovca_Lund,  UBladder)
VALIDATION_TIL$ID <- paste(VALIDATION_TIL$cohort, "#", VALIDATION_TIL$ID)
####




# processing tissue-based files

Lung_tissue <-data.table::fread("Lung_TIL_merged_tissue_seg_data_summary.txt")
Endometrium_tissue1 <-data.table::fread("Endometrium_TIL_1_merged_tissue_seg_data_summary.txt")
Endometrium_tissue2 <-data.table::fread("Endometrium_TIL_2_merged_tissue_seg_data_summary.txt")
Esoph_tissue <-data.table::fread("Esoph_TIL_merged_tissue_seg_data_summary.txt")
Ovca_Lund_tissue <-data.table::fread("OvcaLund_TIL_merged_tissue_seg_data_summary.txt")
UBladder_tissue <-data.table::fread("UroBladder_TIL_merged_tissue_seg_data_summary.txt")
MEL_tissue <-data.table::fread("Melanoma_TIL_merged_tissue_seg_data_summary.txt")

Lung_tissue          <- Lung_tissue[,c(2, 3,10)]
Endometrium_tissue1  <- Endometrium_tissue1[,c(2, 3,10)]
Endometrium_tissue2  <- Endometrium_tissue2[,c(2, 3,10)]
Esoph_tissue         <- Esoph_tissue[,c(2, 3,10)]
Ovca_Lund_tissue     <- Ovca_Lund_tissue[,c(2, 3,10)]
UBladder_tissue      <- UBladder_tissue[,c(2, 3,10)]
MEL_tissue      <- MEL_tissue[,c(2, 3,10)]


Lung_tissue$`Sample Name` <- gsub("]_.*", "]", Lung_tissue$`Sample Name`)
Lung_tissue$`Sample Name` <- gsub('Core...', '[', Lung_tissue$`Sample Name`)

Endometrium_tissue <-rbind(Endometrium_tissue1, Endometrium_tissue2)
Endometrium_tissue$`Sample Name` <- gsub("]_.*", "]", Endometrium_tissue$`Sample Name`)
Endometrium_tissue$`Sample Name` <- gsub('Core...', '[', Endometrium_tissue$`Sample Name`)

Esoph_tissue$`Sample Name` <- gsub("]_.*", "]", Esoph_tissue$`Sample Name`)
Esoph_tissue$`Sample Name` <- gsub('Core...', '[', Esoph_tissue$`Sample Name`)

Ovca_Lund_tissue$`Sample Name` <- gsub("]_.*", "]", Ovca_Lund_tissue$`Sample Name`)
Ovca_Lund_tissue$`Sample Name` <- gsub('Core...', '[', Ovca_Lund_tissue$`Sample Name`)

UBladder_tissue$`Sample Name` <- gsub("]_.*", "]", UBladder_tissue$`Sample Name`)
UBladder_tissue$`Sample Name` <- gsub('Core...', '[', UBladder_tissue$`Sample Name`)

MEL_tissue$`Sample Name` <- gsub("]_.*", "]", MEL_tissue$`Sample Name`)
MEL_tissue$`Sample Name` <- gsub('Core...', '[', MEL_tissue$`Sample Name`)




Lung <- Lung_tissue
Lung$ID <- Lung_match$ID[match(unlist(Lung$`Sample Name`), Lung_match$Sample.Name)]
Lung$`Region Area (pixels)` <- as.numeric(gsub(",", ".", Lung$`Region Area (pixels)`, fixed = T))
Lung$`Region Area (pixels)` <- Lung$`Region Area (pixels)`*(0.497^2)/1000000
names(Lung)[3] <- "Region.Area..mm2."
Lung <- Lung[ which(Lung$`Tissue Category` != 'Blanc'), ]
Lung <- na.omit(Lung)           

Endometrium <- Endometrium_tissue
Endometrium$ID <- Endometrium_match$ID[match(unlist(Endometrium$`Sample Name`), Endometrium_match$Sample.Name)]
Endometrium$`Region Area (pixels)` <- as.numeric(gsub(",", ".", Endometrium$`Region Area (pixels)`, fixed = T))
Endometrium$`Region Area (pixels)` <- Endometrium$`Region Area (pixels)`*(0.497^2)/1000000
names(Endometrium)[3] <- "Region.Area..mm2."
Endometrium <- Endometrium[ which(Endometrium$`Tissue Category` != 'Blanc'), ]
Endometrium <- na.omit(Endometrium)  

Esoph <- Esoph_tissue
Esoph$ID <- Esoph_match$ID[match(unlist(Esoph$`Sample Name`), Esoph_match$Sample.Name)]
Esoph$sample_region <- Esoph_match$sample_region[match(unlist(Esoph$`Sample Name`), Esoph_match$Sample.Name)]
Esoph$`Region Area (pixels)` <- as.numeric(gsub(",", ".", Esoph$`Region Area (pixels)`, fixed = T))
Esoph$`Region Area (pixels)` <- Esoph$`Region Area (pixels)`*(0.497^2)/1000000
names(Esoph)[3] <- "Region.Area..mm2."
Esoph <- Esoph[ which(Esoph$`Tissue Category` != 'Blanc'), ]
Esoph <- na.omit(Esoph)  

Ovca_Lund <- Ovca_Lund_tissue
Ovca_Lund$ID <- Ovca_Lund_match$ID[match(unlist(Ovca_Lund$`Sample Name`), Ovca_Lund_match$Sample.Name)]
Ovca_Lund$sample_region <- Ovca_Lund_match$sample_region[match(unlist(Ovca_Lund$`Sample Name`), Ovca_Lund_match$Sample.Name)]
Ovca_Lund$`Region Area (square microns)` <- as.numeric(gsub(",", ".", Ovca_Lund$`Region Area (square microns)`, fixed = T))
Ovca_Lund$`Region Area (square microns)` <- Ovca_Lund$`Region Area (square microns)`/1000000
names(Ovca_Lund)[3] <- "Region.Area..mm2."
Ovca_Lund <- Ovca_Lund[ which(Ovca_Lund$`Tissue Category` != 'Blanc'), ]
Ovca_Lund <- na.omit(Ovca_Lund)  

UBladder <- UBladder_tissue
UBladder$ID <- UBladder_match$PAD_year[match(unlist(UBladder$`Sample Name`), UBladder_match$Sample.Name)]
UBladder$`Region Area (pixels)` <- as.numeric(gsub(",", ".", UBladder$`Region Area (pixels)`, fixed = T))
UBladder$`Region Area (pixels)` <- UBladder$`Region Area (pixels)`*(0.497^2)/1000000
names(UBladder)[3] <- "Region.Area..mm2."
UBladder <- UBladder[ which(UBladder$`Tissue Category` != 'Blanc'), ]
UBladder <- na.omit(UBladder) 

MEL <- MEL_tissue
MEL$ID <- MEL_match$ID[match(unlist(MEL$`Sample Name`), MEL_match$Sample.Name)]
MEL$`Region Area (square microns)` <- as.numeric(gsub(",", ".", MEL$`Region Area (square microns)`, fixed = T))
MEL$`Region Area (square microns)` <- MEL$`Region Area (square microns)`/1000000
names(MEL)[3] <- "Region.Area..mm2."
MEL <- MEL[ which(MEL$`Tissue Category` != 'Blanc'), ]
MEL <- na.omit(MEL) 


Lung[,"sample_region"]='T'
Lung[,"cohort"]='Lung'

Endometrium[,"sample_region"]='T'
Endometrium[,"cohort"]='Endometrium'

Esoph[,"cohort"]='Esoph'

Ovca_Lund[,"cohort"]='Ovca_Lund'
head(Ovca_Lund)

UBladder[,"sample_region"]='T'
UBladder[,"cohort"]='UBladder'

MEL[,"sample_region"]='T'
MEL[,"cohort"]='MEL'

Merged_tissues_TIL <- rbind(                Lung,
                                            Endometrium,
                                            Esoph,
                                            Ovca_Lund,
                                            UBladder,
                                            MEL)


Merged_tissues_TIL$ID <- paste(Merged_tissues_TIL$cohort, "#", Merged_tissues_TIL$ID)






##################################################
######## *** computing cell densities *** ########

d <-NULL
d <- VALIDATION_TIL
dataTissue <- Merged_tissues_TIL


# remove all samples except ones from primary tumors and remove small samples

dataTissue$sample_region <- ifelse((dataTissue$Region.Area..mm2.>0.005), dataTissue$sample_region, "SMALL")
d$sample_region <- dataTissue$sample_region[match(unlist(d$`Sample Name`), dataTissue$`Sample Name`)]

dataTissue <- dataTissue[ which(dataTissue$sample_region=='T'), ]
d <- d[ which(d$sample_region=='T'), ]


dmeanssum <-NULL

dmeanssum <-data.frame(d[0,c(2:5)]) 
for(i in unique(d$ID)){
  
  d.i <- d[ which(d$ID==i), ]
  dataTissue.i <- dataTissue[ which(dataTissue$ID==i), ]
  
  ###############  
  totalcount<-sum(dataTissue.i$Region.Area..mm2.)
  
  CD8   <- sum(d.i$CD8)/totalcount
  dmeans <-data.frame(cbind(i,  CD8))
  dmeanssum <- data.frame(rbind(dmeanssum, dmeans))
  
}
names(dmeanssum)[c(1:2)]<-c( "Case","CD8")

VALIDATION_TIL <- dmeanssum
head(VALIDATION_TIL)





                                                ##################################################
                                                ############ ***     NK panel    *** ############
                                                ##################################################




##################################################
######### *** reading original files *** #########

##
### Endometrium
##

mergedNK <-data.table::fread("Endometrium_NK_merged_cell_seg_data.txt")

NK <- mergedNK[,c(2,3, 8,9,    72, 77, 82, 87, 92, 97)]  
NK <- na.omit(NK)
names(NK)[c(5:10)]<-c("CD3","CD56","CD68","NKp46","CD163","CK")
NK$`Sample Name` <- gsub("]_.*", "]", NK$`Sample Name`)
NK$`Sample Name` <- gsub('Core...', '[', NK$`Sample Name`)
NK <- NK[ which(NK$`Tissue Category` != 'Blanc'), ]                     
Endometrium <- na.omit(NK)

##
### END Endometrium
##




##
### Esoph
##

mergedNK <-data.table::fread("Esoph_NK_merged_cell_seg_data.txt")

NK <- mergedNK[,c(2,3, 8,9,    72, 77, 82, 87, 92, 97)]  
NK <- na.omit(NK)
names(NK)[c(5:10)]<-c("CD3","CD56","CD68","NKp46","CD163","CK")
NK$`Sample Name` <- gsub("]_.*", "]", NK$`Sample Name`)
NK$`Sample Name` <- gsub('Core...', '[', NK$`Sample Name`)
NK$CD3   <-  as.numeric(gsub(",", ".", NK$CD3, fixed = T))
NK$CD56  <-  as.numeric(gsub(",", ".", NK$CD56, fixed = T))
NK$CD68   <-  as.numeric(gsub(",", ".", NK$CD68, fixed = T))
NK$NKp46 <-  as.numeric(gsub(",", ".", NK$NKp46, fixed = T))
NK$CD163  <-  as.numeric(gsub(",", ".", NK$CD163, fixed = T))
NK$CK    <-  as.numeric(gsub(",", ".", NK$CK, fixed = T))
NK$`Cell X Position` <- as.numeric(gsub(",", ".", NK$`Cell X Position`, fixed = T))
NK$`Cell Y Position` <- as.numeric(gsub(",", ".", NK$`Cell Y Position`, fixed = T))
NK <- NK[ which(NK$`Tissue Category` != 'Blanc'), ]                     
Esoph <- na.omit(NK)

##
### END Esoph
##





##
###  Lung
##

mergedNK <-data.table::fread("Lung_NK_merged_cell_seg_data.txt")

NK <- mergedNK[,c(2,3, 8,9,    72, 77, 82, 87, 92, 97)]  
NK <- na.omit(NK)
names(NK)[c(5:10)]<-c("CD3","CD56","CD68","NKp46","CD163","CK")
NK$`Sample Name` <- gsub("]_.*", "]", NK$`Sample Name`)
NK$`Sample Name` <- gsub('Core...', '[', NK$`Sample Name`)
NK$CD3   <-  as.numeric(gsub(",", ".", NK$CD3, fixed = T))
NK$CD56  <-  as.numeric(gsub(",", ".", NK$CD56, fixed = T))
NK$CD68   <-  as.numeric(gsub(",", ".", NK$CD68, fixed = T))
NK$NKp46 <-  as.numeric(gsub(",", ".", NK$NKp46, fixed = T))
NK$CD163  <-  as.numeric(gsub(",", ".", NK$CD163, fixed = T))
NK$CK    <-  as.numeric(gsub(",", ".", NK$CK, fixed = T))
NK$`Cell X Position` <- as.numeric(gsub(",", ".", NK$`Cell X Position`, fixed = T))
NK$`Cell Y Position` <- as.numeric(gsub(",", ".", NK$`Cell Y Position`, fixed = T))
NK <- NK[ which(NK$`Tissue Category` != 'Blanc'), ]                     
Lung <- na.omit(NK)


##
### END Lung
##





##
###  MEL
##

mergedNK <-data.table::fread("Melanoma_NK_merged_cell_seg_data.txt")

NK <- mergedNK[,c(2,3, 8,9,    72, 77, 82, 87, 92, 97)]  
NK <- na.omit(NK)
names(NK)[c(5:10)]<-c("CD3","CD56","CD68","NKp46","CD163","CK")
NK$`Sample Name` <- gsub("]_.*", "]", NK$`Sample Name`)
NK$`Sample Name` <- gsub('Core...', '[', NK$`Sample Name`)
NK <- NK[ which(NK$`Tissue Category` != 'Blanc'), ]                     
MEL <- na.omit(NK)

##
### END MEL
##


##
###  Ovca
##

mergedNK <-data.table::fread("Ovarian_NK_merged_cell_seg_data.txt")
NK <- mergedNK[,c(2,3, 8,9,    72, 77, 82, 87, 92, 97)]  
NK <- na.omit(NK)
names(NK)[c(5:10)]<-c("CD3","CD56","CD68","NKp46","CD163","CK")
NK$`Sample Name` <- gsub("]_.*", "]", NK$`Sample Name`)
NK$`Sample Name` <- gsub('Core...', '[', NK$`Sample Name`)
NK <- NK[ which(NK$`Tissue Category` != 'Blanc'), ]                     
NK <- na.omit(NK)
Ovca <- na.omit(NK)

##
### END Ovca
##


##
###  UBladder
##

mergedNK <-data.table::fread("UroNK_merged_cell_seg_data.txt")
NK <- mergedNK[,c(2,3, 8,9,    72, 77, 82, 87, 92, 97)]  
NK <- na.omit(NK)
names(NK)[c(5:10)]<-c("CD3","CD56","CD68","NKp46","CD163","CK")
NK$`Sample Name` <- gsub("]_.*", "]", NK$`Sample Name`)
NK$`Sample Name` <- gsub('Core...', '[', NK$`Sample Name`)
NK$`Sample Name` <- gsub('UroBl', 'UBladder', NK$`Sample Name`)
NK <- NK[ which(NK$`Tissue Category` != 'Blanc'), ]                     
UBladder <- na.omit(NK)

##
### END UBladder
##



#######
######
#######



### CK

CK_calculator <- function(fdata){
  require(plyr)
  library(autothresholdr)
  otsu<- as.integer(fdata$CK*1000)
  o_CK_fdata <- auto_thresh(otsu, "Otsu")
  o_CK_fdata <- as.numeric(o_CK_fdata)/1000
  First_distribution <- fdata[ which(fdata$CK  < o_CK_fdata), ]
  Second_distribution <- fdata[ which(fdata$CK  > o_CK_fdata), ]
  peak1 <- density(First_distribution$CK)
  peak2 <- density(Second_distribution$CK)
  # get highest Peak
  peak1 <-peak1$x[which.max(peak1$y)] #Gives you the first highest Peak
  peak2 <-peak2$x[which.max(peak2$y)] #Gives you the second highest Peak
  dY <- density(fdata$CK)$y
  dX <- density(fdata$CK)$x
  MinYDensity<- min(dY[dX  > peak1 & dX < peak2])
  minimum_CK <-dX[which(dY == MinYDensity)]
  
  return(minimum_CK)
}


minimum_CK <- CK_calculator(Endometrium)
Endometrium$CK <- ifelse((Endometrium$CK>minimum_CK), 1, 0)
minimum_CK <- CK_calculator(Esoph)
Esoph$CK <- ifelse((Esoph$CK>minimum_CK), 1, 0)
minimum_CK <- CK_calculator(Lung)
Lung$CK <- ifelse((Lung$CK>minimum_CK), 1, 0)
minimum_CK <- CK_calculator(MEL)
MEL$CK <- ifelse((MEL$CK>minimum_CK), 1, 0)
minimum_CK <- CK_calculator(Ovca)
Ovca$CK <- ifelse((Ovca$CK>minimum_CK), 1, 0)
minimum_CK <- CK_calculator(UBladder)
UBladder$CK <- ifelse((UBladder$CK>minimum_CK), 1, 0)



#Split Ovarian cohort:

Ovca$cohort <-  Ovca$`Sample Name`
Ovca$cohort <- gsub("Ovarian_NK_1C", "Ovarian_Lund_NK_1C", Ovca$cohort)
Ovca$cohort <- gsub("Ovarian_NK_2C", "Ovarian_Lund_NK_2C", Ovca$cohort)
Ovca$cohort <- gsub("Ovarian_NK_3C", "Ovarian_Lund_NK_3C", Ovca$cohort)
Ovca$cohort <- gsub("Ovarian_NK_4C", "Ovarian_Lund_NK_4C", Ovca$cohort)
Ovca$cohort <- gsub("Ovarian_Lund...*", "Ovca_Lund", Ovca$cohort)
Ovca$cohort <- gsub("Tub_NK...*", "Ovca_Lund", Ovca$cohort)
Ovca$cohort <- gsub("Ovarian_NK_2012049T5...*", "Ovca_Sto", Ovca$cohort)
Ovca$cohort <- gsub("Ovarian_NK_BBK877...*", "Ovca_Sto", Ovca$cohort)
Ovca$cohort <- gsub("Ovarian_NK_F049...*", "Ovca_Sto", Ovca$cohort)
Ovca$cohort <- gsub("Ovarian_NK_TMA12...*", "Ovca_Sto", Ovca$cohort)
Ovca$cohort <- gsub("Ovarian_NK_TMA13...*", "Ovca_Sto", Ovca$cohort)


Lung_match <- as.data.frame(read.xlsx2("matching_reference_Lung_NK.xlsx", sheetIndex = 1))
Endometrium_match <- as.data.frame(read.xlsx2("matching_reference_Endometrium_NK.xlsx", sheetIndex = 1))
Esoph_match <- as.data.frame(read.xlsx2("matching_reference_Esoph_NK.xlsx", sheetIndex = 1))
Ovca_Lund_match <- as.data.frame(read.xlsx2("matching_reference_OvarianLund_NK.xlsx", sheetIndex = 1))
UBladder_match <- as.data.frame(read.xlsx2("matching_reference_UroBladder_NK.xlsx", sheetIndex = 1))
MEL_match <- as.data.frame(read.xlsx2("matching_reference_Melanoma_NK.xlsx", sheetIndex = 1))



Lung$ID <- Lung_match$ID[match(unlist(Lung$`Sample Name`), Lung_match$Sample.Name)]
Lung$sample_region <- Lung_match$sample_region[match(unlist(Lung$`Sample Name`), Lung_match$Sample.Name)]
Lung$`Cell X Position` <- Lung$`Cell X Position`*0.497
Lung$`Cell Y Position` <- Lung$`Cell Y Position`*0.497
Lung <- Lung[ which(Lung$`Tissue Category` != 'Blanc'), ]
Lung <- na.omit(Lung)        
NK <- NULL
NK <- Lung

CD3cutoff <-4.64144580785412
CD56cutoff <-2.71856382706077
CD68cutoff <-2.88207019903426
NKp46cutoff <-2.96957425069096
CD163cutoff <-4.25261617609162

NK$CD3 <- ifelse((NK$CD3>CD3cutoff), 1, 0)
NK$CD56 <- ifelse((NK$CD56>CD56cutoff), 1, 0)
NK$CD68 <- ifelse((NK$CD68>CD68cutoff), 1, 0)
NK$NKp46 <- ifelse((NK$NKp46>NKp46cutoff), 1, 0)
NK$CD163 <- ifelse((NK$CD163>CD163cutoff), 1, 0)
Lung <- NK


Endometrium$ID <- Endometrium_match$ID[match(unlist(Endometrium$`Sample Name`), Endometrium_match$Sample.Name)]
Endometrium$`Cell X Position` <- Endometrium$`Cell X Position`*0.497
Endometrium$`Cell Y Position` <- Endometrium$`Cell Y Position`*0.497
Endometrium <- Endometrium[ which(Endometrium$`Tissue Category` != 'Blanc'), ]
Endometrium <- na.omit(Endometrium)  

NK <- NULL
NK <- Endometrium

CD3cutoff <-2.76141161570824
CD56cutoff <-2.71856382706077
CD68cutoff <-2.96957425069096
NKp46cutoff <-3.43475957957958 
CD163cutoff <-1.38386235218325


NK$CD3 <- ifelse((NK$CD3>CD3cutoff), 1, 0)
NK$CD56 <- ifelse((NK$CD56>CD56cutoff), 1, 0)
NK$CD68 <- ifelse((NK$CD68>CD68cutoff), 1, 0)
NK$NKp46 <- ifelse((NK$NKp46>NKp46cutoff), 1, 0)
NK$CD163 <- ifelse((NK$CD163>CD163cutoff), 1, 0)
Endometrium <- NK



Esoph$ID <- Esoph_match$ID[match(unlist(Esoph$`Sample Name`), Esoph_match$Sample.Name)]
Esoph$sample_region <- Esoph_match$sample_region[match(unlist(Esoph$`Sample Name`), Esoph_match$Sample.Name)]
Esoph$`Cell X Position` <- Esoph$`Cell X Position`*0.497
Esoph$`Cell Y Position` <- Esoph$`Cell Y Position`*0.497
Esoph <- Esoph[ which(Esoph$`Tissue Category` != 'Blanc'), ]
Esoph <- na.omit(Esoph)  

NK <- NULL
NK <- Esoph

CD3cutoff <-2.76141161570824
CD56cutoff <-2.71856382706077
CD68cutoff <-2.88207019903426
NKp46cutoff <-2.96957425069096
CD163cutoff <-1.38386235218325


NK$CD3 <- ifelse((NK$CD3>CD3cutoff), 1, 0)
NK$CD56 <- ifelse((NK$CD56>CD56cutoff), 1, 0)
NK$CD68 <- ifelse((NK$CD68>CD68cutoff), 1, 0)
NK$NKp46 <- ifelse((NK$NKp46>NKp46cutoff), 1, 0)
NK$CD163 <- ifelse((NK$CD163>CD163cutoff), 1, 0)
Esoph <- NK



Ovca_Lund <- Ovca[ which(Ovca$cohort == 'Ovca_Lund'), ]
Ovca_Lund <- Ovca_Lund[,c(-11)]
Ovca_Lund$ID <- Ovca_Lund_match$ID[match(unlist(Ovca_Lund$`Sample Name`), Ovca_Lund_match$Sample.Name)]
Ovca_Lund$sample_region <- Ovca_Lund_match$sample_region[match(unlist(Ovca_Lund$`Sample Name`), Ovca_Lund_match$Sample.Name)]
Ovca_Lund$`Cell X Position` <- Ovca_Lund$`Cell X Position`*0.497
Ovca_Lund$`Cell Y Position` <- Ovca_Lund$`Cell Y Position`*0.497
Ovca_Lund <- Ovca_Lund[ which(Ovca_Lund$`Tissue Category` != 'Blanc'), ]
Ovca_Lund <- na.omit(Ovca_Lund)  

NK <- NULL
NK <- Ovca_Lund

CD3cutoff <-2.76141161570824
CD56cutoff <-2.71856382706077
CD68cutoff <-2.88207019903426
NKp46cutoff <-2.96957425069096
CD163cutoff <-3.78511617609163

NK$CD3 <- ifelse((NK$CD3>CD3cutoff), 1, 0)
NK$CD56 <- ifelse((NK$CD56>CD56cutoff), 1, 0)
NK$CD68 <- ifelse((NK$CD68>CD68cutoff), 1, 0)
NK$NKp46 <- ifelse((NK$NKp46>NKp46cutoff), 1, 0)
NK$CD163 <- ifelse((NK$CD163>CD163cutoff), 1, 0)
Ovca_Lund <- NK




UBladder$ID <- UBladder_match$PAD_year[match(unlist(UBladder$`Sample Name`), UBladder_match$Sample.Name)]
UBladder$`Cell X Position` <- UBladder$`Cell X Position`*0.497
UBladder$`Cell Y Position` <- UBladder$`Cell Y Position`*0.497
UBladder <- UBladder[ which(UBladder$`Tissue Category` != 'Blanc'), ]
UBladder <- na.omit(UBladder) 

NK <- NULL
NK <- UBladder

CD3cutoff <-6.43589080785412
CD56cutoff <-2.71856382706077
CD68cutoff <-2.88207019903426
NKp46cutoff <-2.96957425069096
CD163cutoff <-1.383862352

NK$CD3 <- ifelse((NK$CD3>CD3cutoff), 1, 0)
NK$CD56 <- ifelse((NK$CD56>CD56cutoff), 1, 0)
NK$CD68 <- ifelse((NK$CD68>CD68cutoff), 1, 0)
NK$NKp46 <- ifelse((NK$NKp46>NKp46cutoff), 1, 0)
NK$CD163 <- ifelse((NK$CD163>CD163cutoff), 1, 0)
UBladder <- NK


MEL$ID <- MEL_match$ID[match(unlist(MEL$`Sample Name`), MEL_match$Sample.Name)]
MEL$`Cell X Position` <- MEL$`Cell X Position`*0.497
MEL$`Cell Y Position` <- MEL$`Cell Y Position`*0.497
MEL <- MEL[ which(MEL$`Tissue Category` != 'Blanc'), ]
MEL <- na.omit(MEL) 

NK <- NULL
NK <- MEL

CD3cutoff <-2.76141161570824
CD56cutoff <-2.71856382706077
CD68cutoff <-2.88207019903426
NKp46cutoff <-2.96957425069096
CD163cutoff <-1.38386235218325

NK$CD3 <- ifelse((NK$CD3>CD3cutoff), 1, 0)
NK$CD56 <- ifelse((NK$CD56>CD56cutoff), 1, 0)
NK$CD68 <- ifelse((NK$CD68>CD68cutoff), 1, 0)
NK$NKp46 <- ifelse((NK$NKp46>NKp46cutoff), 1, 0)
NK$CD163 <- ifelse((NK$CD163>CD163cutoff), 1, 0)

MEL <- NK



Lung[,"cohort"]='Lung'

Endometrium[,"sample_region"]='T'
Endometrium[,"cohort"]='Endometrium'

Esoph[,"cohort"]='Esoph'

Ovca_Lund[,"cohort"]='Ovca_Lund'

UBladder[,"sample_region"]='T'
UBladder[,"cohort"]='UBladder'

MEL[,"sample_region"]='T'
MEL[,"cohort"]='MEL'


VALIDATION_NK <- rbind(Endometrium, Esoph,  Lung, MEL, Ovca_Lund, UBladder)


VALIDATION_NK$ID <- paste(VALIDATION_NK$cohort, "#", VALIDATION_NK$ID)

head(VALIDATION_NK)







# processing tissue-based files


Lung_tissue <-data.table::fread("Lung_NK_merged_tissue_seg_data_summary.txt")
Endometrium_tissue <-data.table::fread("Endometrium_NK_merged_tissue_seg_data_summary.txt")
Esoph_tissue <-data.table::fread("Esoph_NK_merged_tissue_seg_data_summary.txt")
Ovca_tissue <-data.table::fread("Ovarian_NK_merged_tissue_seg_data_summary.txt")
UBladder_tissue <-data.table::fread("UroNK_merged_tissue_seg_data_summary.txt")
MEL_tissue <-data.table::fread("Melanoma_NK_merged_tissue_seg_data_summary.txt")

Lung_tissue          <- Lung_tissue[,c(2, 3,10)]
Endometrium_tissue  <- Endometrium_tissue[,c(2, 3,10)]
Esoph_tissue         <- Esoph_tissue[,c(2, 3,10)]
Ovca_tissue     <- Ovca_tissue[,c(2, 3,10)]
UBladder_tissue      <- UBladder_tissue[,c(2, 3,10)]
MEL_tissue      <- MEL_tissue[,c(2, 3,10)]

Lung_tissue$`Sample Name` <- gsub("]_.*", "]", Lung_tissue$`Sample Name`)
Lung_tissue$`Sample Name` <- gsub('Core...', '[', Lung_tissue$`Sample Name`)

Endometrium_tissue$`Sample Name` <- gsub("]_.*", "]", Endometrium_tissue$`Sample Name`)
Endometrium_tissue$`Sample Name` <- gsub('Core...', '[', Endometrium_tissue$`Sample Name`)

Esoph_tissue$`Sample Name` <- gsub("]_.*", "]", Esoph_tissue$`Sample Name`)
Esoph_tissue$`Sample Name` <- gsub('Core...', '[', Esoph_tissue$`Sample Name`)

Ovca_tissue$`Sample Name` <- gsub("]_.*", "]", Ovca_tissue$`Sample Name`)
Ovca_tissue$`Sample Name` <- gsub('Core...', '[', Ovca_tissue$`Sample Name`)

UBladder_tissue$`Sample Name` <- gsub("]_.*", "]", UBladder_tissue$`Sample Name`)
UBladder_tissue$`Sample Name` <- gsub('Core...', '[', UBladder_tissue$`Sample Name`)
UBladder_tissue$`Sample Name` <- gsub('UroBl', 'UBladder', UBladder_tissue$`Sample Name`)

MEL_tissue$`Sample Name` <- gsub("]_.*", "]", MEL_tissue$`Sample Name`)
MEL_tissue$`Sample Name` <- gsub('Core...', '[', MEL_tissue$`Sample Name`)

#Split Ovarian cohort:

Ovca_tissue$cohort <-  Ovca_tissue$`Sample Name`

Ovca_tissue$cohort <- gsub("Ovarian_NK_1C", "Ovarian_Lund_NK_1C", Ovca_tissue$cohort)
Ovca_tissue$cohort <- gsub("Ovarian_NK_2C", "Ovarian_Lund_NK_2C", Ovca_tissue$cohort)
Ovca_tissue$cohort <- gsub("Ovarian_NK_3C", "Ovarian_Lund_NK_3C", Ovca_tissue$cohort)
Ovca_tissue$cohort <- gsub("Ovarian_NK_4C", "Ovarian_Lund_NK_4C", Ovca_tissue$cohort)
Ovca_tissue$cohort <- gsub("Ovarian_Lund...*", "Ovca_Lund", Ovca_tissue$cohort)
Ovca_tissue$cohort <- gsub("Tub_NK...*", "Ovca_Lund", Ovca_tissue$cohort)
Ovca_tissue$cohort <- gsub("Ovarian_NK_2012049T5...*", "Ovca_Sto", Ovca_tissue$cohort)
Ovca_tissue$cohort <- gsub("Ovarian_NK_BBK877...*", "Ovca_Sto", Ovca_tissue$cohort)
Ovca_tissue$cohort <- gsub("Ovarian_NK_F049...*", "Ovca_Sto", Ovca_tissue$cohort)
Ovca_tissue$cohort <- gsub("Ovarian_NK_TMA12...*", "Ovca_Sto", Ovca_tissue$cohort)
Ovca_tissue$cohort <- gsub("Ovarian_NK_TMA13...*", "Ovca_Sto", Ovca_tissue$cohort)



########   ****   transform into mm2   **** ########
Lung <- Lung_tissue
Lung$ID <- Lung_match$ID[match(unlist(Lung$`Sample Name`), Lung_match$Sample.Name)]
Lung$`Region Area (pixels)` <- as.numeric(gsub(",", ".", Lung$`Region Area (pixels)`, fixed = T))
Lung$`Region Area (pixels)` <- Lung$`Region Area (pixels)`*(0.497^2)/1000000
names(Lung)[3] <- "Region.Area..mm2."
Lung <- Lung[ which(Lung$`Tissue Category` != 'Blanc'), ]
Lung <- na.omit(Lung)           

Endometrium <- Endometrium_tissue
Endometrium$ID <- Endometrium_match$ID[match(unlist(Endometrium$`Sample Name`), Endometrium_match$Sample.Name)]
Endometrium$`Region Area (pixels)` <- as.numeric(gsub(",", ".", Endometrium$`Region Area (pixels)`, fixed = T))
Endometrium$`Region Area (pixels)` <- Endometrium$`Region Area (pixels)`*(0.497^2)/1000000
names(Endometrium)[3] <- "Region.Area..mm2."
Endometrium <- Endometrium[ which(Endometrium$`Tissue Category` != 'Blanc'), ]
Endometrium <- na.omit(Endometrium)  

Esoph <- Esoph_tissue
Esoph$ID <- Esoph_match$ID[match(unlist(Esoph$`Sample Name`), Esoph_match$Sample.Name)]
Esoph$sample_region <- Esoph_match$sample_region[match(unlist(Esoph$`Sample Name`), Esoph_match$Sample.Name)]
Esoph$`Region Area (pixels)` <- as.numeric(gsub(",", ".", Esoph$`Region Area (pixels)`, fixed = T))
Esoph$`Region Area (pixels)` <- Esoph$`Region Area (pixels)`*(0.497^2)/1000000
names(Esoph)[3] <- "Region.Area..mm2."
Esoph <- Esoph[ which(Esoph$`Tissue Category` != 'Blanc'), ]
Esoph <- na.omit(Esoph)  

Ovca_Lund_tissue <- Ovca_tissue[ which(Ovca_tissue$cohort == 'Ovca_Lund'), ]
Ovca_Lund <- Ovca_Lund_tissue
Ovca_Lund <- Ovca_Lund[,c(-4)]
Ovca_Lund$ID <- Ovca_Lund_match$ID[match(unlist(Ovca_Lund$`Sample Name`), Ovca_Lund_match$Sample.Name)]
Ovca_Lund$sample_region <- Ovca_Lund_match$sample_region[match(unlist(Ovca_Lund$`Sample Name`), Ovca_Lund_match$Sample.Name)]
Ovca_Lund$`Region Area (pixels)` <- as.numeric(gsub(",", ".", Ovca_Lund$`Region Area (pixels)`, fixed = T))
Ovca_Lund$`Region Area (pixels)` <- Ovca_Lund$`Region Area (pixels)`*(0.497^2)/1000000
names(Ovca_Lund)[3] <- "Region.Area..mm2."
Ovca_Lund <- Ovca_Lund[ which(Ovca_Lund$`Tissue Category` != 'Blanc'), ]
Ovca_Lund <- na.omit(Ovca_Lund)  

UBladder <- UBladder_tissue
UBladder$ID <- UBladder_match$PAD_year[match(unlist(UBladder$`Sample Name`), UBladder_match$Sample.Name)]
UBladder$`Region Area (pixels)` <- as.numeric(gsub(",", ".", UBladder$`Region Area (pixels)`, fixed = T))
UBladder$`Region Area (pixels)` <- UBladder$`Region Area (pixels)`*(0.497^2)/1000000
names(UBladder)[3] <- "Region.Area..mm2."
UBladder <- UBladder[ which(UBladder$`Tissue Category` != 'Blanc'), ]
UBladder <- na.omit(UBladder) 

MEL <- MEL_tissue
MEL$ID <- MEL_match$ID[match(unlist(MEL$`Sample Name`), MEL_match$Sample.Name)]
MEL$`Region Area (pixels)` <- as.numeric(gsub(",", ".", MEL$`Region Area (pixels)`, fixed = T))
MEL$`Region Area (pixels)` <- MEL$`Region Area (pixels)`*(0.497^2)/1000000
names(MEL)[3] <- "Region.Area..mm2."
MEL <- MEL[ which(MEL$`Tissue Category` != 'Blanc'), ]
MEL <- na.omit(MEL) 


Lung[,"sample_region"]='T'
Lung[,"cohort"]='Lung'

Endometrium[,"sample_region"]='T'
Endometrium[,"cohort"]='Endometrium'

Esoph[,"cohort"]='Esoph'

Ovca_Lund[,"cohort"]='Ovca_Lund'

UBladder[,"sample_region"]='T'
UBladder[,"cohort"]='UBladder'

MEL[,"sample_region"]='T'
MEL[,"cohort"]='MEL'

Merged_tissues_NK <- rbind(
  Lung,
  Endometrium,
  Esoph,
  Ovca_Lund,
  UBladder,
  MEL)

Merged_tissues_NK$ID <- paste(Merged_tissues_NK$cohort, "#", Merged_tissues_NK$ID)





##################################################
######## *** computing cell densities *** ########

d <- NULL
d <- VALIDATION_NK

dataTissue <- Merged_tissues_NK


d$СD3  <- ifelse((VALIDATION_NK$CD3==1), 1, 0)
d$M2  <- ifelse((VALIDATION_NK$CD68==1 & VALIDATION_NK$CD163==1), 1, 0)


# filter 

dataTissue$sample_region <- ifelse((dataTissue$Region.Area..mm2.>0.005), dataTissue$sample_region, "SMALL")
d$sample_region <- dataTissue$sample_region[match(unlist(d$`Sample Name`), dataTissue$`Sample Name`)]

dataTissue <- dataTissue[ which(dataTissue$sample_region=='T'), ]
d <- d[ which(d$sample_region=='T'), ]

# filter end

dmeanssum <-NULL
dmeanssum <-data.frame(d[0,c(2:5)]) 
for(i in unique(d$ID)){
  
  d.i <- d[ which(d$ID==i), ]
  dataTissue.i <- dataTissue[ which(dataTissue$ID==i), ]
  
  ###############  
  totalcount<-sum(dataTissue.i$Region.Area..mm2.)
  
  СD3 <- sum(d.i$CD3)/totalcount
  M2  <- sum(d.i$M2)/totalcount
  
  dmeans <-data.frame(cbind(i, СD3, M2 ))
  
  dmeanssum <- data.frame(rbind(dmeanssum, dmeans))
  
}
names(dmeanssum)[c(1:3)]<-c("Case",  "СD3",   "M2")

VALIDATION_NK <- dmeanssum
head(VALIDATION_NK)
###




                  ##################################################
                  ########### *** merge to clinical *** ############
                  ##################################################


########################

library(dplyr)


##### match 

df <- full_join(VALIDATION_TIL, VALIDATION_NK,by="Case")

head(df)



Lung_clin <- as.data.frame(read.xlsx2("Lung_clinical.xlsx", sheetIndex = 1))
Endometrium_clin <- as.data.frame(read.xlsx2("Endometrium_clinical.xlsx", sheetIndex = 1))
Esoph_clin <- as.data.frame(read.xlsx2("Esophagus_clinical.xlsx", sheetIndex = 1))
Ovca_Lund_clin <- as.data.frame(read.xlsx2("Ovca_Lund_clinical.xlsx", sheetIndex = 1))
UBladder_clin <- as.data.frame(read.xlsx2("UBladder_clinical.xlsx", sheetIndex = 1))
MEL_clin <- as.data.frame(read.xlsx2("MEL_clinical.xlsx", sheetIndex = 1))


Lung_clin[,1]<-Lung_clin[,4]
Endometrium_clin[,"TMA_cohort"]='Endometrium'
Esoph_clin[,"TMA_cohort"]='Esoph'
Ovca_Lund_clin[,"TMA_cohort"]='Ovca_Lund'
UBladder_clin[,"TMA_cohort"]='UBladder'
MEL_clin[,"TMA_cohort"]='MEL'

Lung_clin[,"Tumor_type"]='Lung'
Endometrium_clin[,"Tumor_type"]='Endometrium'
Esoph_clin[,"Tumor_type"]='Esoph'
Ovca_Lund_clin[,"Tumor_type"]='Ovca_Lund'
UBladder_clin[,"Tumor_type"]='UBladder'
MEL_clin[,"Tumor_type"]='MEL'


Merged_clinical <- data.frame(rbind(Lung_clin, Endometrium_clin,Esoph_clin,
                                    Ovca_Lund_clin, UBladder_clin, MEL_clin))


names(Merged_clinical)[2]<-paste("Case")

Merged_clinical$Case <- paste(Merged_clinical$Tumor_type, "#", Merged_clinical$Case)


library(tidyverse)

library(dplyr)

Merged <- as.data.frame(left_join(Merged_clinical, df,by="Case"))
Merged <- Merged[ which(Merged$exclude != 'exclude'), ]
Merged <- Merged[ which(Merged$Tumor_type_code != ''), ]
Merged[] <- Map(function(x) replace(x, is.infinite(x), NA), Merged)

setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub")
write.table(Merged, file ="Database_working_Validation_cohort.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)



