########################
########################
setwd("/Volumes/Samsung_X5/ENCIMA")

require(gdata)
require(xlsx)
library(scales) 
library(ggplot2)
library(reshape)
library(plyr)
library(robustbase)
library(data.table)
library(WriteXLS)
library(tidyverse)
library(readxl)
library(DescTools)
require("survival")
library("survminer")
library(survival)
library(schoRsch)
library("pracma")






setwd("/Volumes/Samsung_X5/ENCIMA/cohort output/")

cohort <- NULL

cohort <- data.table::fread("Merged_cohort.txt")

#cohort$CD8  <- as.numeric(nthroot(cohort$CD8,4))
#cohort$M2  <- as.numeric(nthroot(cohort$M2,4))

cohort$SIA <- (cohort$CD8)/(cohort$CD8+cohort$M2) 
#cohort$SIA  <- as.numeric(nthroot(cohort$SIA,4))
#cohort$SIA  <- cohort$SIA^8
#cohort$SIA <- (cohort$CD8-cohort$M2)/(cohort$CD8+cohort$M2)


unique(cohort$Tumor_type_code)
cohort <- cohort[ which(cohort$Tumor_type_code != 'DCIS'), ]
#cohort <- cohort[ which(cohort$exclude != 'exclude'), ]
#cohort <- cohort[ which(cohort$PreOp_treatment_yesno != 'Yes'), ]
#cohort <- cohort[ which(cohort$pT_stge != '0'), ]
#cohort <- cohort[ which(cohort$clinicl_stge != '0'), ]
#cohort <-cohort[ which(!is.na(cohort$Tumor_type_code)), ]

cohort <-cohort[ which(cohort$Tumor_type_code == 'MEL'), ]


library(tidyr) # for gather
library(ggplot2)
colnames(cohort)



########
########            OS
########
datS<-as.data.frame(cohort)
str(datS)
datS <- datS[ which(datS$Time_Diagnosis_Last_followup..OS. != "missing"), ]
datS$Time_Diagnosis_Last_followup..OS. <- as.numeric(datS$Time_Diagnosis_Last_followup..OS.)
datS <- datS[ which(!is.na(datS$Time_Diagnosis_Last_followup..OS.)), ]

datS$status <-datS$Event_last_followup
datS$SurvObj <- with(datS, Surv(datS$Time_Diagnosis_Last_followup..OS., status == "Dead"))
datS <- datS[ which(datS$Time_Diagnosis_Last_followup..OS. > 0), ]
head(datS)
unique(datS$TMA_cohort)

########
########            TTR
########

#datS<-as.data.frame(cohort)
#str(datS)

   # datS$TTR <- as.numeric(datS$TTR)
  #  datS$statusTTR <-datS$TTR.event
  #  datS$SurvObj <- with(datS, Surv(datS$TTR, statusTTR == "Yes"))
   # datS <- datS[ which(datS$TTR > 0), ]


  #  datS$RFS <- as.numeric(datS$RFS)
  #  datS$statusRFS <-datS$RFS.event
  #  datS$SurvObj <- with(datS, Surv(datS$RFS, statusRFS == "Yes"))
  #  datS <- datS[ which(datS$RFS > 0), ]
    

    
###########
##########
###########
        
    library(tidyr) # for gather
    library(ggplot2)  
           
dat <- datS
dat <- dat[!is.na(dat$CD8),]
dat <- dat[!is.na(dat$M2),]

colnames(dat)

for(i in names(dat[,c( 40,45,49,54, 66, 71)] )){
  
  
  x  <-  as.data.frame(dat[!is.na(dat[i]),])
  
  #x[i] <- ntiles(x, dv = i, bins=3)
  str(x[i])
  mean<- mean(x[[i]])
  med <- median(x[[i]])
  print(i)
  #print(med)
  print(mean)
  unique(x[i])
  x[i] <-  ifelse(( x[i]>med), 1, 0)
  
  #  x[i] <-  ifelse(( x[i]>mean), 2, 1)
  
  
  f <- "Case"
  
  
  
  dat <- left_join( dat,  x[,c(f,i)], by = f, all = T)
  
  
}






colnames(dat)


d <- dat[,c(2,75:ncol(dat))]
colnames(d) 
str(d)
library(reshape)
d$Case <- as.character(d$Case)
d_melt <-reshape2::melt(d, id="Case")
head(d_melt)


markers<-unique(d_melt$variable)
cases<-unique(d_melt$Case)




myplots <- list()

merged <- NULL
for(j in unique(cases))
{
  
  #j<-87
  J.list <-NULL
  J.list <- as.data.frame(d_melt[ which(d_melt$Case==j), ])
  J.list
  
  scores <- NULL
  for(i in unique(markers))
  {
    
    #i<-"CD4_Single" 
    x <-NULL
    x <- J.list[ which(J.list$variable==i), ]
    x
    
    score <- sum(x$value, na.rm = T)
    #score <- max(x$value)
    #score <- mean(x$value, na.rm = T)
    devider <- length(which(!is.na(x$value)))
    score <- score/devider*2
    
    
    score <-  round(score,digits=5)
    
    
    marker_ID  <-i
    Case  <-j
    score <- data.frame(Case,marker_ID, score)
    
    scores <- data.frame(rbind(scores, score))  
  }
  
  
  merged <- data.frame(rbind(merged, scores))
  
}



final <-cast(merged, Case~marker_ID)



dat <- left_join( dat,  final, by = "Case", all = T)






        #dev.off()
        par(mfrow=c(3,6))

        datF <- dat


        
        datF$Tumor_type_code
        
        fit <- survfit(SurvObj ~ SIA_tri, data = datF, conf.type = "log-log")

        c<-coxph(SurvObj ~ SIA_tri, data = datF)
        c<-summary(c)$coefficients[5]
        c
        c<-as.character(round(c, digits = 15))
        c<-paste("p=", c, sep = "")
        legend(x="bottomleft",legend=c, bty = "n",cex=1.5)
        
        ##
        
        ggsurvplot(fit, risk.table = T, lwd=1,pval = TRUE,palette = "jco")
        

        #dev.off()   
        par(mfrow=c(4,4))
