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
setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub/datasets/CRC/")



matching_reference  <- read.table("matching_reference_TIL.txt", sep="\t", header=T, fill = T)

CRC_1 <-data.table::fread("merged_cell_seg_data_TIL_1.txt")
CRC_1 <- CRC_1[,c(2,3, 8,9, 72,77,82,41,92,97)]  
CRC_1 <- na.omit(CRC_1)
names(CRC_1)[c(5:10)]<-c("CD4","CD20","CD8","FoxP3","CD45RO","CK")
CRC_1$`Sample Name` <- gsub('UCAN 07', 'UCAN C07', CRC_1$`Sample Name`)
CRC_1$`Sample Name` <- gsub('UCAN ', '', CRC_1$`Sample Name`)
CRC_1$`Sample Name` <- gsub("]_.*", "]", CRC_1$`Sample Name`)
CRC_1$`Sample Name` <- gsub('Core', 'C', CRC_1$`Sample Name`)


CRC_2 <-data.table::fread("merged_cell_seg_data_TIL_2.txt")
CRC_2 <- CRC_2[,c(2,3, 8,9, 72,77,82,41,92,97)] 
CRC_2 <- na.omit(CRC_2)
names(CRC_2)[c(5:10)]<-c("CD4","CD20","CD8","FoxP3","CD45RO","CK")
CRC_2$`Sample Name` <- gsub('Core...', '[', CRC_2$`Sample Name`)
CRC_2$`Sample Name` <- gsub('UCAN ', '', CRC_2$`Sample Name`)
CRC_2$`Sample Name` <- gsub("]_.*", "]", CRC_2$`Sample Name`)



############ *** applying cutoffs


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
        
minimum_CK_1 <- CK_calculator(CRC_1)
minimum_CK_2 <- CK_calculator(CRC_2)

CRC_1$CK <- ifelse((CRC_1$CK>minimum_CK_1), 1, 0)
CRC_2$CK <- ifelse((CRC_2$CK>minimum_CK_2), 1, 0)


# cutofs

CD4cutoff_1 <-2.173293187044
CD4cutoff_2 <-2.173293187044

CD20cutoff_1 <-1.350903035
CD20cutoff_2 <-1.350903035

CD8cutoff_1 <-1.732516267
CD8cutoff_2 <-1.732516267

FoxP3cutoff_1 <-2.132880561
FoxP3cutoff_2 <-2.132880561

CD45ROcutoff_1 <-31
CD45ROcutoff_2 <-42


### CD4
CRC_1$CD4 <- ifelse((CRC_1$CD4>CD4cutoff_1), 1, 0)
CRC_2$CD4 <- ifelse((CRC_2$CD4>CD4cutoff_2), 1, 0)

### CD8
CRC_1$CD8 <- ifelse((CRC_1$CD8>CD8cutoff_1), 1, 0)
CRC_2$CD8 <- ifelse((CRC_2$CD8>CD8cutoff_2), 1, 0)

### CD20
CRC_1$CD20 <- ifelse((CRC_1$CD20>CD20cutoff_1), 1, 0)
CRC_2$CD20 <- ifelse((CRC_2$CD20>CD20cutoff_2), 1, 0)

### CD45RO
CRC_1$CD45RO <- ifelse((CRC_1$CD45RO>CD45ROcutoff_1), 1, 0)
CRC_2$CD45RO <- ifelse((CRC_2$CD45RO>CD45ROcutoff_2), 1, 0)

### FoxP3
CRC_1$FoxP3 <- ifelse((CRC_1$FoxP3>FoxP3cutoff_1), 1, 0)
CRC_2$FoxP3 <- ifelse((CRC_2$FoxP3>FoxP3cutoff_2), 1, 0)


##### merge
CRC<-rbind(CRC_1, CRC_2)
rm(CRC_1)
rm(CRC_2)
####
write.table(CRC, file ="Merged_all_cells_CRC_TIL.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
CRC <- NULL
CRC <- read.table("Merged_all_cells_CRC_TIL.txt", sep="\t", header=T, fill = T)

########   ****   transform coordinates into um   **** ########
CRC$Cell.X.Position <- CRC$Cell.X.Position *0.497
CRC$Cell.Y.Position <- CRC$Cell.Y.Position *0.497
names(CRC)[3] <- "Cell.X.Position.um"
names(CRC)[4] <- "Cell.Y.Position.um"

##### match core names
            CRC$Core.ID <- CRC$Sample.Name
            CRC$Sample.Name <- matching_reference$case_region[match(unlist(CRC$Sample.Name), matching_reference$Sample.Name)]
            CRC$Sample <-CRC$Sample.Name
            CRC$Sample.Region <-CRC$Sample.Name
            CRC$Sample <- gsub("_.*", "", CRC$Sample)
            CRC$Sample.Region <- gsub(".*_", "", CRC$Sample.Region)
            CRC <- CRC[,c(11,1,12,13, 2:10)]
            CRC <- CRC[ which(CRC$Tissue.Category != 'Blanc'), ]
            CRC <- na.omit(CRC)


            
            
            
########################
# processing tissue-based files


temp_1<- data.table::fread("merged_tissue_seg_data_summary_TIL_1.txt", sep="\t", header=T, fill = T)
temp_2<- data.table::fread("merged_tissue_seg_data_summary_TIL_2.txt", sep="\t", header=T, fill = T)
temp_1 <- temp_1 [,c(2, 3,10)]
temp_2 <- temp_2 [,c(2, 3,10)]
tissuemerged <-rbind(temp_1, temp_2)
tissuemerged <- na.omit(tissuemerged)

tissuemerged$`Sample Name` <- gsub('UCAN 07', 'UCAN C07', tissuemerged$`Sample Name`)
tissuemerged$`Sample Name` <- gsub('UCAN ', '', tissuemerged$`Sample Name`)
tissuemerged$`Sample Name` <- gsub("]_.*", "]", tissuemerged$`Sample Name`)
tissuemerged$`Sample Name` <- gsub('Core', 'C', tissuemerged$`Sample Name`)


tissuemerged$`Region Area (pixels)` <- as.numeric(tissuemerged$`Region Area (pixels)`)
###   ****   transform into mm2   **** ########
              tissuemerged$`Region Area (pixels)` <- tissuemerged$`Region Area (pixels)`*(0.497^2)/1000000
              names(tissuemerged)[3] <- "Region.Area..mm2."

##### match core names
tissuemerged$Core.ID <-tissuemerged$`Sample Name`
tissuemerged$`Sample Name` <- matching_reference$case_region[match(unlist(tissuemerged$`Sample Name`), matching_reference$Sample.Name)]
tissuemerged$Sample <-tissuemerged$`Sample Name`
tissuemerged$Sample.Region <-tissuemerged$`Sample Name`
tissuemerged$Sample <- gsub("_.*", "", tissuemerged$Sample)
tissuemerged$Sample.Region <- gsub(".*_", "", tissuemerged$Sample.Region)
tissuemerged <- tissuemerged[ which(tissuemerged$`Tissue Category` != 'Blanc'), ]
tissuemerged <- na.omit(tissuemerged)
tissuemerged[order(tissuemerged$Region.Area..mm2.) , ]
tissuemerged <- tissuemerged[ which(tissuemerged$Region.Area..mm2. > 0), ]





##################################################
######## *** comluting cell densities *** ########


data <- NULL
data <- CRC
d <- NULL
d <- data


d$CD4_Single  <- ifelse((data$CD20==0 & data$CD4==1 & data$CD8==0 & data$FoxP3==0 & data$CD45RO==0), 1, 0)
d$CD4_activated  <- ifelse((data$CD20==0 & data$CD4==1 & data$CD8==0 & data$FoxP3==0 & data$CD45RO==1), 1, 0)
d$CD4_Treg  <- ifelse((data$CD20==0 & data$CD4==1 & data$CD8==0 & data$FoxP3==1), 1, 0)
d$CD8_Single  <- ifelse((data$CD20==0 & data$CD4==0 & data$CD8==1 & data$FoxP3==0 & data$CD45RO==0), 1, 0)
d$CD8_activated  <- ifelse((data$CD20==0 & data$CD4==0 & data$CD8==1 & data$FoxP3==0 & data$CD45RO==1), 1, 0)
d$CD8_Treg  <- ifelse((data$CD20==0 & data$CD4==0 & data$CD8==1 & data$FoxP3==1), 1, 0)
d$B_cells  <- ifelse((data$CD20==1 & data$CD4==0 & data$CD8==0 & data$FoxP3==0), 1, 0)

dmeanssum_ALL <-NULL
dmeanssum_ALLs <-NULL

dataTissue <- tissuemerged
# FOR tot
dmeanssum <-data.frame(d[0,c(2:5)]) 
for(i in unique(d$Sample)){
  
  d.i <- d[ which(d$Sample==i), ]
  #d.i <- d.i[ which(d.i$Tissue.Category=='Stroma' | d.i$Tissue.Category=='Tumor'), ]
  d.i <- d.i[ which(d.i$Sample.Region == 'Center'| d.i$Sample.Region=='Peri'), ]
  
  dataTissue.i <- dataTissue[ which(dataTissue$Sample==i), ]
  dataTissue.i <- dataTissue.i[ which(dataTissue.i$Sample.Region == 'Center'| dataTissue.i$Sample.Region == 'Peri'), ]
  
  ###############  
  totalcount<-sum(dataTissue.i$Region.Area..mm2.)
  
  CD4 <- sum(d.i$CD4)/totalcount
  CD4_Single <- sum(d.i$CD4_Single)/totalcount
  CD4_activated   <- sum(d.i$CD4_activated)/totalcount
  CD4_Treg  <- sum(d.i$CD4_Treg)/totalcount
  CD8 <- sum(d.i$CD8)/totalcount
  CD8_Single<- sum(d.i$CD8_Single)/totalcount
  CD8_activated   <- sum(d.i$CD8_activated)/totalcount
  CD8_Treg  <- sum(d.i$CD8_Treg)/totalcount
  CD45RO   <- sum(d.i$CD45RO)/totalcount
  B_cells   <- sum(d.i$B_cells)/totalcount
  
  dmeans <-data.frame(cbind(i,
                            
                            
                            CD4,
                            CD4_Single,
                            CD4_activated,
                            CD4_Treg,
                            CD8,
                            CD8_Single,
                            CD8_activated,
                            CD8_Treg,
                            CD45RO,
                            B_cells
                            
                            
  ))
  
  colnames(dmeans)
  dmeanssum <- data.frame(rbind(dmeanssum, dmeans))
  
}

names(dmeanssum)[c(1:11)]<-c( "Case",
                              "CD4","CD4_Single","CD4_activated", "CD4_Treg", 
                              "CD8", "CD8_Single", "CD8_activated", "CD8_Treg", "CD45RO","B_cells")
dmeanssum_ALL <- dmeanssum


# FOR Center_T


dmeanssum <-data.frame(d[0,c(2:5)]) 
for(i in unique(d$Sample)){
  d.i <- d[ which(d$Sample==i), ]
  #d.i <- d.i[ which(d.i$Tissue.Category=='Stroma' | d.i$Tissue.Category=='Tumor'), ]
  d.i <- d.i[ which(d.i$Sample.Region == 'Center'), ]
  
  dataTissue.i <- dataTissue[ which(dataTissue$Sample==i), ]
  dataTissue.i <- dataTissue.i[ which(dataTissue.i$Sample.Region == 'Center'), ]
  
  ############### 
  totalcount<-sum(dataTissue.i$Region.Area..mm2.)
  
  CD4 <- sum(d.i$CD4)/totalcount
  CD4_Single <- sum(d.i$CD4_Single)/totalcount
  CD4_activated   <- sum(d.i$CD4_activated)/totalcount
  CD4_Treg  <- sum(d.i$CD4_Treg)/totalcount
  CD8 <- sum(d.i$CD8)/totalcount
  CD8_Single<- sum(d.i$CD8_Single)/totalcount
  CD8_activated   <- sum(d.i$CD8_activated)/totalcount
  CD8_Treg  <- sum(d.i$CD8_Treg)/totalcount
  CD45RO   <- sum(d.i$CD45RO)/totalcount
  B_cells   <- sum(d.i$B_cells)/totalcount
  
  dmeans <-data.frame(cbind(i,
                            
                            CD4,
                            CD4_Single,
                            CD4_activated,
                            CD4_Treg,
                            CD8,
                            CD8_Single,
                            CD8_activated,
                            CD8_Treg,
                            CD45RO,
                            B_cells
                            
  ))
  
  colnames(dmeans)
  dmeanssum <- data.frame(rbind(dmeanssum, dmeans))
}

names(dmeanssum)[c(1:11)]<-c( "Case",
                              "CD4_CT","CD4_Single_CT","CD4_activated_CT", "CD4_Treg_CT", 
                              "CD8_CT", "CD8_Single_CT", "CD8_activated_CT", "CD8_Treg_CT", "CD45RO_CT","B_cells_CT")

dmeanssum_ALL <- data.frame(cbind(dmeanssum_ALL, dmeanssum))


# FOR Center_S

dmeanssum <-data.frame(d[0,c(1:2)])    
for(i in unique(d$Sample)){
  d.i <- d[ which(d$Sample==i), ]
  #d.i <- d.i[ which(d.i$Tissue.Category=='Stroma' | d.i$Tissue.Category=='Tumor'), ]
  d.i <- d.i[ which(d.i$Sample.Region == 'Peri'), ]
  dataTissue.i <- dataTissue[ which(dataTissue$Sample==i), ]
  dataTissue.i <- dataTissue.i[ which(dataTissue.i$Sample.Region == 'Peri'), ]
  
  ###############  
  totalcount<-sum(dataTissue.i$Region.Area..mm2.)
  
  CD4 <- sum(d.i$CD4)/totalcount
  CD4_Single <- sum(d.i$CD4_Single)/totalcount
  CD4_activated   <- sum(d.i$CD4_activated)/totalcount
  CD4_Treg  <- sum(d.i$CD4_Treg)/totalcount
  CD8 <- sum(d.i$CD8)/totalcount
  CD8_Single<- sum(d.i$CD8_Single)/totalcount
  CD8_activated   <- sum(d.i$CD8_activated)/totalcount
  CD8_Treg  <- sum(d.i$CD8_Treg)/totalcount
  CD45RO   <- sum(d.i$CD45RO)/totalcount
  B_cells   <- sum(d.i$B_cells)/totalcount
  
  dmeans <-data.frame(cbind(i,
                            
                            CD4,
                            CD4_Single,
                            CD4_activated,
                            CD4_Treg,
                            CD8,
                            CD8_Single,
                            CD8_activated,
                            CD8_Treg,
                            CD45RO,
                            B_cells
                            
  ))
  
  colnames(dmeans)
  dmeanssum <- data.frame(rbind(dmeanssum, dmeans))
}

names(dmeanssum)[c(1:11)]<-c( "Case",
                              "CD4_IM","CD4_Single_IM","CD4_activated_IM", "CD4_Treg_IM", 
                              "CD8_IM", "CD8_Single_IM", "CD8_activated_IM", "CD8_Treg_IM", "CD45RO_IM","B_cells_IM")

dmeanssum_ALL <- data.frame(cbind(dmeanssum_ALL, dmeanssum))
CRC_TIL <- dmeanssum_ALL[,c(-12, -23)]

setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub/datasets/CRC/")
write.table(CRC_TIL, file ="CRC_TIL.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)






                                                ##################################################
                                                ############ ***     NK panel    *** ############
                                                ##################################################




##################################################
######### *** reading original files *** #########
matching_reference  <- read.table("matching_reference_NK.txt", sep="\t", header=T, fill = T)


CRC_NK <-data.table::fread("merged_cell_seg_data_NK.txt")
CRC_NK <- CRC_NK[,c(2,3, 8,9, 72,77,82,87,92,97)]  
CRC_NK <- na.omit(CRC_NK)
names(CRC_NK)[c(5:10)]<-c("CD3","CD56","CD68","NKp46","CD163","CK")
CRC_NK$`Sample Name` <- gsub('UCAN CRC 07', 'UCAN CRC C07', CRC_NK$`Sample Name`)
CRC_NK$`Sample Name` <- gsub('UCAN CRC c', 'UCAN CRC C', CRC_NK$`Sample Name`)
CRC_NK$`Sample Name` <- gsub('UCAN CRC ', '', CRC_NK$`Sample Name`)
CRC_NK$`Sample Name` <- gsub("]_.*", "]", CRC_NK$`Sample Name`)
CRC_NK$`Sample Name` <- gsub('.....Core', 'Core', CRC_NK$`Sample Name`)
CRC_NK$`Sample Name` <- gsub(' Core...', '_[', CRC_NK$`Sample Name`)
CRC_NK$CD3 <- as.numeric(gsub(",", ".", CRC_NK$CD3, fixed = T))
CRC_NK$CD56 <- as.numeric(gsub(",", ".", CRC_NK$CD56, fixed = T))
CRC_NK$CD68 <- as.numeric(gsub(",", ".", CRC_NK$CD68, fixed = T))
CRC_NK$NKp46 <- as.numeric(gsub(",", ".", CRC_NK$NKp46, fixed = T))
CRC_NK$CD163 <- as.numeric(gsub(",", ".", CRC_NK$CD163, fixed = T))
CRC_NK$CK <- as.numeric(gsub(",", ".", CRC_NK$CK, fixed = T))
CRC_NK$`Cell X Position` <- as.numeric(gsub(",", ".", CRC_NK$`Cell X Position`, fixed = T))
CRC_NK$`Cell Y Position` <- as.numeric(gsub(",", ".", CRC_NK$`Cell Y Position`, fixed = T))

####
CRC <- na.omit(CRC_NK)
rm(CRC_NK)

############ *** applying cutoffs

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

minimum_CK <- CK_calculator(CRC)
CRC$CK <- ifelse((CRC$CK>minimum_CK), 1, 0)

# cutofs

CD3cutoff <-2.761411616
CD56cutoff<-2.718563827
CD68cutoff<-2.882070199
NKp46cutoff<-2.969574251
CD163cutoff<-1.383862352


### CD3
CRC$CD3 <- ifelse((CRC$CD3>CD3cutoff), 1, 0)

### CD56
CRC$CD56 <- ifelse((CRC$CD56>CD56cutoff), 1, 0)

### CD68
CRC$CD68 <- ifelse((CRC$CD68>CD68cutoff), 1, 0)

### CDNKp46
CRC$NKp46 <- ifelse((CRC$NKp46>NKp46cutoff), 1, 0)

### CD163
CRC$CD163 <- ifelse((CRC$CD163>CD163cutoff), 1, 0)


############## **** ############

mean(CRC$CD3)
mean(CRC$CD56)
mean(CRC$NKp46)
mean(CRC$CD68)
mean(CRC$CD163)
mean(CRC$CK)

####
write.table(CRC, file ="Merged_all_cells_CRC_NK.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
CRC <- NULL
CRC <- read.table("Merged_all_cells_CRC_NK.txt", sep="\t", header=T, fill = T)

########   ****   transform coordinates into um   **** ########
CRC$Cell.X.Position <- CRC$Cell.X.Position *0.496
CRC$Cell.Y.Position <- CRC$Cell.Y.Position *0.496
names(CRC)[3] <- "Cell.X.Position.um"
names(CRC)[4] <- "Cell.Y.Position.um"


##### match core names
CRC$Core.ID <- CRC$Sample.Name
CRC$Sample.Name <- matching_reference$case_region[match(unlist(CRC$Sample.Name), matching_reference$Sample.Name)]
CRC$Sample <-CRC$Sample.Name
CRC$Sample.Region <-CRC$Sample.Name
CRC$Sample <- gsub("_.*", "", CRC$Sample)
CRC$Sample.Region <- gsub(".*_", "", CRC$Sample.Region)
CRC <- CRC[,c(11,1,12,13, 2:10)]
CRC <- CRC[ which(CRC$Tissue.Category != 'Blanc'), ]
CRC <- na.omit(CRC)




# processing tissue-based files


########################

temp<- data.table::fread("merged_tissue_seg_data_summary_NK.txt", sep="\t", header=T, fill = T)
tissuemerged <- temp [,c(2, 3,10)]
tissuemerged <- na.omit(tissuemerged)
tissuemerged$`Sample Name` <- gsub('UCAN CRC 07', 'UCAN CRC C07', tissuemerged$`Sample Name`)
tissuemerged$`Sample Name` <- gsub('UCAN CRC c', 'UCAN CRC C', tissuemerged$`Sample Name`)
tissuemerged$`Sample Name` <- gsub('UCAN CRC ', '', tissuemerged$`Sample Name`)
tissuemerged$`Sample Name` <- gsub("]_.*", "]", tissuemerged$`Sample Name`)
tissuemerged$`Sample Name` <- gsub('.....Core', 'Core', tissuemerged$`Sample Name`)
tissuemerged$`Sample Name` <- gsub(' Core...', '_[', tissuemerged$`Sample Name`)
tissuemerged$`Region Area (pixels)` <- as.numeric(tissuemerged$`Region Area (pixels)`)

########   ****   transform into mm2   **** ########
tissuemerged$`Region Area (pixels)` <- tissuemerged$`Region Area (pixels)`*(0.496^2)/1000000
names(tissuemerged)[3] <- "Region.Area..mm2."

##### match core names
tissuemerged$Core.ID <-tissuemerged$`Sample Name`
tissuemerged$`Sample Name` <- matching_reference$case_region[match(unlist(tissuemerged$`Sample Name`), matching_reference$Sample.Name)]
tissuemerged$Sample <-tissuemerged$`Sample Name`
tissuemerged$Sample.Region <-tissuemerged$`Sample Name`
tissuemerged$Sample <- gsub("_.*", "", tissuemerged$Sample)
tissuemerged$Sample.Region <- gsub(".*_", "", tissuemerged$Sample.Region)
tissuemerged <- tissuemerged[ which(tissuemerged$`Tissue Category` != 'Blanc'), ]
tissuemerged <- na.omit(tissuemerged)
tissuemerged[order(tissuemerged$Region.Area..mm2.) , ]
tissuemerged <- tissuemerged[ which(tissuemerged$Region.Area..mm2. > 0), ]
#########




##################################################
######## *** comluting cell densities *** ########

data <- NULL
data <- CRC
d <- NULL
d <- data


d$CD3  <- ifelse((data$CD3==1), 1, 0)
d$NK  <- ifelse((data$CD3==0 & data$CD56==1 & data$NKp46==1), 1, 0)
d$NKT  <- ifelse((data$CD3==1 & data$CD56==1 & data$NKp46==1), 1, 0)
d$M1  <- ifelse((data$CD68==1 & data$CD163==0), 1, 0)
d$Myeloid  <- ifelse((data$CD68==0 & data$CD163==1), 1, 0)
d$M2  <- ifelse(( data$CD68==1 & data$CD163==1), 1, 0)
d$CK_single  <- ifelse((data$CD3==0 & data$CD56==0 & data$NKp46==0 & data$CD68==0 & data$CD163==0 & data$CK==1), 1, 0)


dmeanssum_ALL <-NULL
dmeanssum_ALLs <-NULL

dataTissue <- tissuemerged

# FOR _tot
dmeanssum <-data.frame(d[0,c(2:5)]) 
for(i in unique(d$Sample)){
  d.i <- d[ which(d$Sample==i), ]
  #d.i <- d.i[ which(d.i$Tissue.Category=='Stroma' | d.i$Tissue.Category=='Tumor'), ]
  d.i <- d.i[ which(d.i$Sample.Region == 'Center'| d.i$Sample.Region=='Peri'), ]
  dataTissue.i <- dataTissue[ which(dataTissue$Sample==i), ]
  dataTissue.i <- dataTissue.i[ which(dataTissue.i$Sample.Region == 'Center'| dataTissue.i$Sample.Region == 'Peri'), ]
  
  ###############  
  totalcount<-sum(dataTissue.i$Region.Area..mm2.)
  
  CD3 <- sum(d.i$CD3)/totalcount
  NK <- sum(d.i$NK)/totalcount
  NKT <- sum(d.i$NKT)/totalcount
  CD68  <- sum(d.i$CD68)/totalcount
  CD163  <- sum(d.i$CD163)/totalcount
  M1<- sum(d.i$M1)/totalcount
  Myeloid  <- sum(d.i$Myeloid)/totalcount
  M2  <- sum(d.i$M2)/totalcount
  
  dmeans <-data.frame(cbind(i,
                            CD3,
                            NK,
                            NKT,
                            CD68,
                            CD163,
                            M1,
                            Myeloid,
                            M2
                            
  ))
  dmeanssum <- data.frame(rbind(dmeanssum, dmeans))
  
}

names(dmeanssum)[c(1:9)]<-c( "Case",
                             "CD3",
                             "NK",
                             "NKT",
                             "CD68",
                             "CD163",
                             "M1",
                             "Myeloid",
                             "M2"   )

dmeanssum_ALL <- dmeanssum

# FOR Center
dmeanssum <-data.frame(d[0,c(2:5)]) 
for(i in unique(d$Sample)){
  
  d.i <- d[ which(d$Sample==i), ]
  #d.i <- d.i[ which(d.i$Tissue.Category=='Stroma' | d.i$Tissue.Category=='Tumor'), ]
  d.i <- d.i[ which(d.i$Sample.Region == 'Center'), ]
  dataTissue.i <- dataTissue[ which(dataTissue$Sample==i), ]
  dataTissue.i <- dataTissue.i[ which(dataTissue.i$Sample.Region == 'Center'), ]
  
  ###############  
  totalcount<-sum(dataTissue.i$Region.Area..mm2.)
  
  CD3 <- sum(d.i$CD3)/totalcount
  NK <- sum(d.i$NK)/totalcount
  NKT <- sum(d.i$NKT)/totalcount
  CD68  <- sum(d.i$CD68)/totalcount
  CD163  <- sum(d.i$CD163)/totalcount
  M1<- sum(d.i$M1)/totalcount
  Myeloid  <- sum(d.i$Myeloid)/totalcount
  M2  <- sum(d.i$M2)/totalcount
  
  dmeans <-data.frame(cbind(i,
                            
                            CD3,
                            NK,
                            NKT,
                            CD68,
                            CD163,
                            M1,
                            Myeloid,
                            M2
                            
  ))
  
  dmeanssum <- data.frame(rbind(dmeanssum, dmeans))
}

names(dmeanssum)[c(1:9)]<-c( "Case",
                             "CD3_CT",
                             "NK_CT",
                             "NKT_CT",
                             "CD68_CT",
                             "CD163_CT",
                             "M1_CT",
                             "Myeloid_CT",
                             "M2_CT"   )


dmeanssum_ALL <- data.frame(cbind(dmeanssum_ALL, dmeanssum))

# FOR Peri
dmeanssum <-data.frame(d[0,c(2:5)]) 
for(i in unique(d$Sample)){
  
  d.i <- d[ which(d$Sample==i), ]
  #d.i <- d.i[ which(d.i$Tissue.Category=='Stroma' | d.i$Tissue.Category=='Tumor'), ]
  d.i <- d.i[ which(d.i$Sample.Region == 'Peri'), ]
  dataTissue.i <- dataTissue[ which(dataTissue$Sample==i), ]
  dataTissue.i <- dataTissue.i[ which(dataTissue.i$Sample.Region == 'Peri'), ]
  
  ###############  
  totalcount<-sum(dataTissue.i$Region.Area..mm2.)
  
  CD3 <- sum(d.i$CD3)/totalcount
  NK <- sum(d.i$NK)/totalcount
  NKT <- sum(d.i$NKT)/totalcount
  CD68  <- sum(d.i$CD68)/totalcount
  CD163  <- sum(d.i$CD163)/totalcount
  M1<- sum(d.i$M1)/totalcount
  Myeloid  <- sum(d.i$Myeloid)/totalcount
  M2  <- sum(d.i$M2)/totalcount
  
  dmeans <-data.frame(cbind(i,
                            
                            CD3,
                            NK,
                            NKT,
                            CD68,
                            CD163,
                            M1,
                            Myeloid,
                            M2
                            
  ))
  
  dmeanssum <- data.frame(rbind(dmeanssum, dmeans))
  
  
}

names(dmeanssum)[c(1:9)]<-c( "Case_IM",
                             "CD3_IM",
                             "NK_IM",
                             "NKT_IM",
                             "CD68_IM",
                             "CD163_IM",
                             "M1_IM",
                             "Myeloid_IM",
                             "M2_IM"   )

dmeanssum_ALL <- data.frame(cbind(dmeanssum_ALL, dmeanssum))

CRC_NK <- dmeanssum_ALL[,c(-10, -19)]
setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub/datasets/CRC/")
write.table(CRC_NK, file ="CRC_NK.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)




                                    ##################################################
                                    ########### *** merge to clinical *** ############
                                    ##################################################
library(dplyr)

setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub/datasets/CRC/")
dat <- NULL
dat <- read.xlsx2("CRC clinical data.xlsx", sheetIndex = 1)
CRC_NK <- read.table("CRC_NK.txt", sep="\t", header=T, fill = T)
CRC_TIL <- read.table("CRC_TIL.txt", sep="\t", header=T, fill = T)
##### match 

df <- full_join(CRC_TIL, CRC_NK,by="Case")
df$Case <- as.factor(df$Case)

names(dat)[16]<-"Case"
data <- left_join(dat, df,by="Case", all = F)

######  make immunoscore-like metric
dat2 <- data
colnames(dat2)
for(i in names(dat2[,c(215,225,239,247)] )){
  x  <-  as.data.frame(dat2[!is.na(dat2[i]),])
  x[i] <- (rank(x[[i]]) - 1) / (length(x[[i]])-1)
  f <- "UCAN_TMA_Patnr"
  dat2 <- merge( dat2,  x[,c(f,i)], by = f, all = T)
}

dat2$ISLike <- rowMeans(dat2[,c("CD8_CT.y","CD8_IM.y","CD3_CT.y","CD3_IM.y")], na.rm=TRUE)
head(dat2)
dat2$ISLike_temp <-  dat2$ISLike
dat2$ISLike_temp <-  ifelse((dat2$ISLike>0.25), 2, 1)
dat2$ISLike <-  ifelse((dat2$ISLike>0.70), 3, dat2$ISLike_temp)
data <- dat2[,c(-255,-256,-257,-258,-260)]

head(data)
#### END of make immunoscore

setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub")
write.table(data, file ="Database_working_CRC.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)




