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

#cohort$SIA <- (cohort$CD8_T)/(cohort$CD8_T+cohort$M2_T) + (cohort$CD8_S)/(cohort$CD8_S+cohort$M2_S)
 
#cohort$SIA  <- as.numeric(nthroot(cohort$SIA,4))
#cohort$SIA  <- cohort$SIA^8
#cohort$SIA <- (cohort$CD8-cohort$M2)/(cohort$CD8+cohort$M2)

colnames(cohort)
write.table(cohort, file ="Multi_tumor_database.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
unique(cohort$TMA_cohort)
cohort <- cohort[ which(cohort$Tumor_type_code != 'DCIS'), ]
cohort <- cohort[ which(cohort$exclude != 'exclude'), ]
#cohort <- cohort[ which(cohort$PreOp_treatment_yesno != 'Yes'), ]
#cohort <- cohort[ which(cohort$pT_stge != '4'), ]
#cohort <- cohort[ which(cohort$clinicl_stge != '4'), ]
#cohort <-cohort[ which(!is.na(cohort$Tumor_type_code)), ]
#cohort <-cohort[ which((cohort$TMA_cohort=="Esoph" )), ]
#cohort <-cohort[ which((cohort$Tumor_type_code=="STAD" )), ]
#cohort <-cohort[ which((cohort$TMA_cohort=="UBladder" | cohort$TMA_cohort=="Lung" | cohort$TMA_cohort=="MEL" | cohort$TMA_cohort=="Esoph" )), ]
#cohort <-cohort[ which((cohort$TMA_cohort=="UBladder" | cohort$TMA_cohort=="Lung" | cohort$TMA_cohort=="MEL"| cohort$TMA_cohort=="Kidney")), ]
cohort <-cohort[ which((cohort$TMA_cohor=="Ovca_Lund" | cohort$TMA_cohort=="UBladder" | cohort$TMA_cohort=="Lung" | cohort$TMA_cohort=="MEL"| cohort$TMA_cohort=="Kidney"| cohort$TMA_cohort=="Prostate")), ]

unique (cohort$TMA_cohort)
library(dplyr)
df_mod <-cohort %>%
  group_by(Tumor_type_code) %>%
  summarise(Count = n())



#Draw pie
pie(df_mod$Count, labels = c(df_mod$Tumor_type_code,df_mod$Count))


library(plotrix)

lbls <- paste(df_mod$Tumor_type_code, "\n", df_mod$Count, sep="")
pie3D(df_mod$Count,labels=lbls,explode=0.1,
      main="Pie Chart of Countries ")

pie(df_mod$Count, labels = lbls, 
    main="Pie Chart of Species\n (with sample sizes)")


#k <- ggplot(cohort, aes( x = Fraction, color = Tumor_type_code)) + stat_ecdf()

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
datS$status   <- ifelse((datS$status =="Dead"), 1, 0)
datS$SurvObj <- with(datS, Surv(datS$Time_Diagnosis_Last_followup..OS., status == "1"))
datS <- datS[ which(datS$Time_Diagnosis_Last_followup..OS. > 0), ]
head(datS)
unique(datS$TMA_cohort)



########
########            RFS
########

#datS<-as.data.frame(cohort)
#str(datS)



    datS$RFS <- as.numeric(datS$RFS)
    datS$status <-datS$RFS.event
    datS$status   <- ifelse((datS$status =="Yes"), 1, 0)
    datS$SurvObj <- with(datS, Surv(datS$RFS, status == "1"))
    datS <- datS[ which(datS$RFS > 0), ]
    

    
###########
##########
###########
        


datS$SIA <- datS$CD8/(datS$CD8+datS$M2)
   # datS$SIA <- (datS$CD8)/(datS$CD8+datS$M2+datS$iDC)

#datS$SIA <- datS$iDC
#datS$SIA <- (datS$CD8+datS$iDC)/(datS$M2)


    library(tidyr) # for gather
    library(ggplot2)  
           

        #dev.off()
        par(mfrow=c(3,6))

        datF <- datS
        datF <-datF[ which(!is.na(datF$SIA)), ]
        datF$SIA_tri <- ntiles(datF, dv = "SIA", bins=3)

        
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

        library("devtools")

        datF <- datS  
        
                                    for(i in unique(datF$TMA_cohort)){                            
                                        
                                        d<-NULL
                                        d <- datF[ which(datF$TMA_cohort == i), ]
                                        head(d)
                                        colnames(d)
                                        d<-d[,c(81:84)]
                                        d <- na.omit(d)
                                        par(mar=c(2,2,2,2))
                                        
                                        d$SIA_tri <- ntiles(d, dv = "SIA", bins=3)
                                        
                                        quant_cut <- function(x, n) {
                                          qs <- quantile(x, 1:(n-1)/n)
                                          brks <- c(-Inf, qs, Inf)
                                          cut(x, breaks=brks, labels=FALSE)
                                        }
                                        
                                       # qs <- quantile(d$SIA, 1:(3-1)/3)
                                        
                                        rank(d$SIA)
                    
                                        #d$SIA_tri <-quant_cut(d$SIA, 3)
                                        
                                        
                                        #d$SIA_tri <-  ifelse(( d$SIA_tri>2), 2, 1)
                                       #d$SIA_tri <-  ifelse(( d$SIA>mean(d$SIA)), 2, 1)
                                        
                                      #    med <- mean(d$M2)
                                      #    d$M2bin <-  ifelse(( d$M2>med), 2, 1)
                                      #  med <- mean(d$CD8)
                                      #  d$CD8bin <-  ifelse(( d$CD8>med), 2, 1)
                                      #  d$SIA_tri<- d$CD8bin/(d$CD8bin+d$M2bin)
                                        
                                        km <- survfit(SurvObj ~ SIA_tri, data = d, conf.type = "log-log")
                                        plot(km, mark.time = T, lwd=2, col=1:5,
                                             main=i)
                                        #fit <- survfit(SurvObj ~ SIA_tri, data = d, conf.type = "log-log")
                                       #ggsurvplot(km, risk.table = F, pval = TRUE)
                                        
                                        c<-coxph(SurvObj ~ SIA_tri, data = d)
                                        c<-summary(c)$coefficients[5]
                                        c
                                        c<-as.character(round(c, digits = 5))
                                        c<-paste("p=", c, sep = "")
                                        legend(x="bottomleft",legend=c, bty = "n",cex=1.5)
                                        
                                        #d$SIA_tri <- ifelse((d$SIA_tri>2), 1, 0)
                                       # d$status <- ifelse((d$status=="No"), 1, 0)
                                       # d$status <- ifelse((d$status=="Alive"), 1, 0)
                                       # PRROC_obj <- roc.curve(scores.class0 = d$SIA, weights.class0=d$status,
                                                 #              curve=TRUE)
                                       # plot(PRROC_obj)      

                                      }
                                      
                                      legend(x="topright", col=1:5, lwd=2, legend=(1:3))
                                      
                         
                                      
                                      
                                      
                                      
                                      
                                      
                                      ############################################
                                      ############################################            Univariate
                                      ###########################################
                                      
                                      par(mfrow=c(1,1))
                                      

                                      datF2 <- datS
                                      
                                      
                                
                                      
                                      
                                      for(i in unique(datF2$TMA_cohort)){                            
                                        
                                        d<-NULL
                                        d <- datF2[ which(datF$TMA_cohort == i), ]
                                        head(d)
                                        colnames(d)
                                        d<-d[,c(81:83)]
                                        d <- na.omit(d)
                                        par(mar=c(2,2,2,2))
                                        
                                        fit0 <- coxph(SurvObj ~ SIA,data = d, na.action=na.omit )
                                 
                                        print(fit0)
                                      }
                                      
                                      
                                      
                                      
    #######
                                      ####
                                      #####
                                      #####
                                      
                                      
                                      
        
        ############# dat[ which(dat$UCAN_TMA_Patnr <520),  ]
                                      
                                      datF2 <- datF[ which(datF$TMA_cohort == "Lung"),  ]
                                      datF2 <- datF[ which(datF$TMA_cohort == "MEL"),  ]
                                      datF2 <- datF[ which(datF$TMA_cohort == "UBladder"),  ]

                                      
                                      dat <-  datF2
                                      colnames(dat)     
                                      head(dat)    
                                      x  <-  dat[,c(2,83)]
                                      x  <-  na.omit(x)
                                      head(x)
                                      x$SIA_tri <- ntiles(x, dv = "SIA", bins=3)
                                      x<- x[-2]
                                      dat <- right_join(dat, x,by="Case")

                                      datF2 <- dat
  
                                      

        
        # https://stats.stackexchange.com/questions/86535/bootstrap-to-evaluate-variance-of-auc-roc
        library(pROC)
        library(boot)
        library(MASS)
                                      library(risksetROC)
        i <- NULL
       
        
        
        test <- datF2$Time_Diagnosis_Last_followup..OS.
        #test <- datF2$RFS
        #score <- datF2$SIA_tri
        score <- as.factor(as.character(datF2$SIA_tri))
        #score <- as.numeric(datF2$SIA)
        #score <- as.factor(as.character(datF2$pT_stge))

        dx <- datF2$status
        df <- data.frame(cbind(test, dx,score)) 
        str(df)
        r <- 1000
        boot.f2 <- function(d, i){
          data <- d[i,]
          #data <- df[1,]
          print(i)
          survival.time <- data$test
          score <- data$score
          survival.status <- data$dx
          surv.prob <- unique(survfit(Surv(survival.time,survival.status)~1)$surv)
          fit0 <- coxph( Surv(survival.time,survival.status)
                         ~ score , na.action=na.omit )
          eta <- fit0$linear.predictor
          model.score <- eta
          
          utimes <- unique( survival.time[ survival.status == 1 ] )
          utimes <- utimes[ order(utimes) ]
          
          ## find AUC at unique failure times
          AUC <- rep( NA, length(utimes) )
          for( j in 1:length(utimes) )
          {
            out <- CoxWeights( eta, survival.time, survival.status,utimes[j])
            AUC[j] <- out$AUC
          }
          ## integrated AUC to get concordance measure
          print(IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" ))
          IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" )
          
          
          
        } 
        roc.boot<-NULL
        roc.boot<-boot(df, boot.f2, r) 
        roc.boot 
       median(roc.boot$t)
        
                               
       
       
       
       
       
       
       
       ##########
       ##########
       ##########
       ##########                          iAUC for all factors separatelly
       #                                                                           https://rdrr.io/cran/risksetROC/man/IntegrateAUC.html
       
  
       
       
       
                                         ##########             M  E  L  A  N  O  M  A   
                                         ##########
                                         ##########
                                         ##########                          iAUC for all factors separatelly
                                         #                                                                           https://rdrr.io/cran/risksetROC/man/IntegrateAUC.html
                                                                                  
                                         library(pROC)
                                         library(boot)
                                         library(MASS)
                                         library(risksetROC)
                                         
                                         datF <- datS 
                                         
                                         datF2 <- datF[ which(datF$TMA_cohort == "MEL"),  ]

                                         
                                         
                                         dat <-  datF2
                                         colnames(dat)     
                                         head(dat)    
                                         x  <-  dat[,c(2,83)]
                                         x  <-  na.omit(x)
                                         head(x)
                                         x$SIA_tri <- ntiles(x, dv = "SIA", bins=3)
                                         x<- x[-2]
                                         dat <- right_join(dat, x,by="Case")
                                         
                                         datF2 <- dat
                                         
                                         
                                         
                                         ### age dichotomisation
                                         dat <-  datF2
                                         colnames(dat)     
                                         head(dat)    
                                         x  <-  dat[,c(2,14)]
                                         x  <-  na.omit(x)
                                         head(x)
                                         med <- NULL
                                         med <- median(as.numeric(x$Age))
                                         print(i)
                                         print(med)
                                   
                                       
                                         x$Age <-  ifelse(( x$Age >med), 2, 1)
                                         
                                         dat <- right_join(dat, x,by="Case")
                                         
                                         datF2 <- dat
                                         
                                         
                                         
                                         datF2$SIA_tri <- as.factor(datF2$SIA_tri)
                                         datF2$SIA_tri <- relevel(datF2$SIA_tri, ref = "1")
                                         
                                         f <- cph(SurvObj ~ SIA_tri + 
                                                    pT_stge+Gender+Age.y,
                                                  x=TRUE, y=TRUE, data=datF2)
                                         plot(anova(f), what='proportion chisq')
                                         
                                         
                                         
                                         factors<-colnames(datF2[,c("Gender","pT_stge","Age.y" , "SIA_tri")])
                                         roc.boot_ALL <- NULL
                                         for(j in unique(factors)){
                                           
                                           datF3<-NULL
                                           datF3 <- datF2[j]
                                           str(datF2)
                                           datF3$test <- datF2$Time_Diagnosis_Last_followup..OS.
                                          # datF2$test <- datF2$RFS
                                           
                                           
                                           datF3$dx <-datF2$status
                                           
                                           
                                           datF3 <- na.omit(datF3)
                                           
                                           
                                           # test <- datF2$OS_surgery_to_last_followup
                                           
                                           head(datF3)
                                           # variables
                                           names(datF3)[1] <- "score"
                                           datF3$test <- as.numeric(datF3$test)
                                           datF3$dx <- as.numeric(datF3$dx)
                                           datF3$score <- as.factor(datF3$score)
                                           colnames(datF3)
                                           
                                           
                                           
                                           str(datF3)
                                           r <- 1000
                                           boot.f2 <- function(d, i){
                                             data <- d[i,]
                                             #data <- datF2[i,]
                                             print(i)
                                             survival.time <- data$test
                                             score <- data$score
                                             survival.status <- data$dx
                                             surv.prob <- unique(survfit(Surv(survival.time,survival.status)~1)$surv)
                                             fit0 <- coxph( Surv(survival.time,survival.status)
                                                            ~ score, na.action=na.omit )
                                             eta <- fit0$linear.predictor
                                             model.score <- eta
                                             
                                             utimes <- unique( survival.time[ survival.status == 1 ] )
                                             utimes <- utimes[ order(utimes) ]
                                             
                                             ## find AUC at unique failure times
                                             AUC <- rep( NA, length(utimes) )
                                             for( j in 1:length(utimes) )
                                             {
                                               out <- CoxWeights( eta, survival.time, survival.status,utimes[j])
                                               AUC[j] <- out$AUC
                                             }
                                             ## integrated AUC to get concordance measure
                                             print(IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" ))
                                             IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" )
                                             
                                             
                                             
                                           } 
                                           roc.boot<-NULL
                                           roc.boot<-boot(datF3, boot.f2, r) 
                                           roc.boot 
                                           
                                           roc.boot_ALL[[j]] <- roc.boot  
                                           
                                         }
                                         roc.boot_ALL
                                         Sex<- roc.boot_ALL$Gender$t
                                         clinicl_stge<- roc.boot_ALL$clinicl_stge$t
                                         T_stage<- roc.boot_ALL$pT_stge$t
                                         N_stage<- roc.boot_ALL$pN_stage$t
                                         Age<- roc.boot_ALL$Age.y$t
                                         Tumor_type<- roc.boot_ALL$Tumor_type$t
                                         #Adjuvant_treatment<- roc.boot_ALL$T_Adjuvant$t
                                         SIA<- roc.boot_ALL$SIA_tri$t     
                                         
                                         
                                         median(SIA)
                                         
                                         Factprs_AUC <- cbind(Sex,Age,T_stage,SIA)
                                         Factprs_AUC <- as.data.frame(Factprs_AUC)
                                         names(Factprs_AUC)[c(1:4)]<-c("Sex","Age","T_stage","SIA")
                                         Factprs_AUC<-tibble::rownames_to_column(Factprs_AUC)
                                         head(Factprs_AUC)                                      
                                         
                                         All_clinical_together_and_SIA
                                         df_plot <- reshape2::melt(Factprs_AUC, id.vars="rowname")
                                         head(df_plot)
                                         
                                         df_plot$variable <- fct_rev(df_plot$variable)
                                         
                                         library(scales)
                                         library(RColorBrewer)
                                         
                                         
                                         k <- ggplot(df_plot, aes(x=reorder(variable, desc(variable)), y=value, fill=variable))+
                                           scale_fill_brewer(palette="RdBu")
                                         
                                         
                                         
                                         k+ 
                                           theme_bw()+ theme(panel.border = element_blank())+ theme(axis.line = element_line(colour = "black"))+
                                           # geom_point(aes(colour=variable,group = variable,y=value), size=2.3, alpha=1, position=position_dodge(width=0))+
                                           geom_boxplot(outlier.shape = NA, stat = "boxplot", width=0.7, alpha=1, position=position_dodge(0.9))+
                                           theme(legend.position = "none",
                                                 text = element_text(size=10),
                                                 axis.text.x=element_text(size=rel(2), angle=45))
                                         
       
       
       
       
       
                                         
                                         



##########             L   U   N   G   
##########
##########
##########                          iAUC for all factors separatelly
#                                                                           https://rdrr.io/cran/risksetROC/man/IntegrateAUC.html


library(pROC)
library(boot)
library(MASS)
library(risksetROC)

datF <- datS 

datF2 <- datF[ which(datF$TMA_cohort == "Lung"),  ]



dat <-  datF2
colnames(dat)     
head(dat)    
x  <-  dat[,c(2,83)]
x  <-  na.omit(x)
head(x)
x$SIA_tri <- ntiles(x, dv = "SIA", bins=3)
x<- x[-2]
dat <- right_join(dat, x,by="Case")

datF2 <- dat

f <- cph(SurvObj ~ SIA_tri + 
           clinicl_stge+Gender,
         x=TRUE, y=TRUE, data=datF2)
plot(anova(f), what='proportion chisq')


### age dichotomisation
dat <-  datF2
colnames(dat)     
head(dat)    
x  <-  dat[,c(2,14)]
x  <-  na.omit(x)
head(x)
x$Age <- as.numeric(x$Age)
med <- median(x$Age)
print(i)
print(med)
print(mean)

x$Age <-  ifelse(( x$Age >med), 2, 1)

dat <- right_join(dat, x,by="Case")

datF2 <- dat






factors<-colnames(datF2[,c("Gender", "clinicl_stge","Smoking","WHO_PS","pT_stge","pN_stage","Age.y" ,"SIA_tri")])
roc.boot_ALL <- NULL
for(j in unique(factors)){

datF3<-NULL
datF3 <- datF2[j]
str(datF2)
datF3$test <- datF2$Time_Diagnosis_Last_followup..OS.
#datF2$test <- datF2$RFS


datF3$dx <-datF2$status


datF3 <- na.omit(datF3)


# test <- datF2$OS_surgery_to_last_followup

head(datF3)
# variables
names(datF3)[1] <- "score"
datF3$test <- as.numeric(datF3$test)
datF3$dx <- as.numeric(datF3$dx)
datF3$score <- as.factor(datF3$score)
colnames(datF3)



str(datF3)
r <- 1000
boot.f2 <- function(d, i){
 data <- d[i,]
 #data <- datF2[i,]
 print(i)
 survival.time <- data$test
 score <- data$score
 survival.status <- data$dx
 surv.prob <- unique(survfit(Surv(survival.time,survival.status)~1)$surv)
 fit0 <- coxph( Surv(survival.time,survival.status)
                ~ score, na.action=na.omit )
 eta <- fit0$linear.predictor
 model.score <- eta
 
 utimes <- unique( survival.time[ survival.status == 1 ] )
 utimes <- utimes[ order(utimes) ]
 
 ## find AUC at unique failure times
 AUC <- rep( NA, length(utimes) )
 for( j in 1:length(utimes) )
 {
   out <- CoxWeights( eta, survival.time, survival.status,utimes[j])
   AUC[j] <- out$AUC
 }
 ## integrated AUC to get concordance measure
 print(IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" ))
 IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" )
 
 
 
} 
roc.boot<-NULL
roc.boot<-boot(datF3, boot.f2, r) 
roc.boot 

roc.boot_ALL[[j]] <- roc.boot  

}
roc.boot_ALL
Sex<- roc.boot_ALL$Gender$t
clinicl_stge<- roc.boot_ALL$clinicl_stge$t
T_stage<- roc.boot_ALL$pT_stge$t
N_stage<- roc.boot_ALL$pN_stage$t
Age<- roc.boot_ALL$Age.y$t
Tumor_type<- roc.boot_ALL$Tumor_type$t
#Adjuvant_treatment<- roc.boot_ALL$T_Adjuvant$t
SIA<- roc.boot_ALL$SIA_tri$t     
Smoking<- roc.boot_ALL$Smoking$t
WHO_PS<- roc.boot_ALL$WHO_PS$t
Sex
Tumor_type
Age
Smoking
WHO_PS
T_stage
N_stage
SIA
median(SIA)

Factprs_AUC <- cbind(Sex,Age,Smoking,WHO_PS,T_stage,N_stage,SIA)
Factprs_AUC <- as.data.frame(Factprs_AUC)
names(Factprs_AUC)[c(1:7)]<-c("Sex","Age","Smoking","WHO_PS","T_stage","N_stage","SIA")
Factprs_AUC<-tibble::rownames_to_column(Factprs_AUC)
head(Factprs_AUC)                                      

All_clinical_together_and_SIA
df_plot <- reshape2::melt(Factprs_AUC, id.vars="rowname")
head(df_plot)

df_plot$variable <- fct_rev(df_plot$variable)

library(scales)
library(RColorBrewer)


k <- ggplot(df_plot, aes(x=reorder(variable, desc(variable)), y=value, fill=variable))+
scale_fill_brewer(palette="RdBu")



k+ 
theme_bw()+ theme(panel.border = element_blank())+ theme(axis.line = element_line(colour = "black"))+
# geom_point(aes(colour=variable,group = variable,y=value), size=2.3, alpha=1, position=position_dodge(width=0))+
geom_boxplot(outlier.shape = NA, stat = "boxplot", width=0.7, alpha=1, position=position_dodge(0.9))+
theme(legend.position = "none",
     text = element_text(size=10),
     axis.text.x=element_text(size=rel(2), angle=45))









                                              
                                              ##########             B L A D D E R    
                                              ##########
                                              ##########
                                              ##########                          iAUC for all factors separatelly
                                              #                                                                           https://rdrr.io/cran/risksetROC/man/IntegrateAUC.html
                                              
                                              
                                              library(pROC)
                                              library(boot)
                                              library(MASS)
                                              library(risksetROC)
                                              
                                              datF <- datS 
                                              datF2 <- datF[ which(datF$TMA_cohort == "UBladder"),  ]

                                              
                                              
                                              
                                              dat <-  datF2
                                              colnames(dat)     
                                              head(dat)    
                                              x  <-  dat[,c(2,83)]
                                              x  <-  na.omit(x)
                                              head(x)
                                              x$SIA_tri <- ntiles(x, dv = "SIA", bins=3)
                                              x<- x[-2]
                                              dat <- right_join(dat, x,by="Case")
                                              
                                              datF2 <- dat
                                              
                                              f <- cph(SurvObj ~ SIA_tri + 
                                                         pT_stge+pN_stage+Gender,
                                                       x=TRUE, y=TRUE, data=datF2)
                                              plot(anova(f), what='proportion chisq')
                                              
                                              ### age dichotomisation
                                              dat <-  datF2
                                              colnames(dat)     
                                              head(dat)    
                                              x  <-  dat[,c(2,14)]
                                              x  <-  na.omit(x)
                                              head(x)
                                              med <- NULL
                                              med <- median(as.numeric(x$Age))
                                              print(i)
                                              print(med)
                                              
                                        
                                              x$Age <-  ifelse(( x$Age >med), 2, 1)
                                              
                                              dat <- right_join(dat, x,by="Case")
                                              
                                              datF2 <- dat
                                              
                                              
                                              
                                              
                                              
                                              
                                              factors<-colnames(datF2[,c("Gender","Diff_grade","pT_stge","Age.y" , "SIA_tri")])
                                              roc.boot_ALL <- NULL
                                              for(j in unique(factors)){
                                                
                                                datF3<-NULL
                                                datF3 <- datF2[j]
                                                str(datF2)
                                                datF3$test <- datF2$Time_Diagnosis_Last_followup..OS.
                                                #datF2$test <- datF2$RFS
                                                
                                                
                                                datF3$dx <-datF2$status
                                                
                                                
                                                datF3 <- na.omit(datF3)
                                                
                                                
                                                # test <- datF2$OS_surgery_to_last_followup
                                                
                                                head(datF3)
                                                # variables
                                                names(datF3)[1] <- "score"
                                                datF3$test <- as.numeric(datF3$test)
                                                datF3$dx <- as.numeric(datF3$dx)
                                                datF3$score <- as.factor(datF3$score)
                                                colnames(datF3)
                                                
                                                
                                                
                                                str(datF3)
                                                r <- 1000
                                                boot.f2 <- function(d, i){
                                                  data <- d[i,]
                                                  #data <- datF2[i,]
                                                  print(i)
                                                  survival.time <- data$test
                                                  score <- data$score
                                                  survival.status <- data$dx
                                                  surv.prob <- unique(survfit(Surv(survival.time,survival.status)~1)$surv)
                                                  fit0 <- coxph( Surv(survival.time,survival.status)
                                                                 ~ score, na.action=na.omit )
                                                  eta <- fit0$linear.predictor
                                                  model.score <- eta
                                                  
                                                  utimes <- unique( survival.time[ survival.status == 1 ] )
                                                  utimes <- utimes[ order(utimes) ]
                                                  
                                                  ## find AUC at unique failure times
                                                  AUC <- rep( NA, length(utimes) )
                                                  for( j in 1:length(utimes) )
                                                  {
                                                    out <- CoxWeights( eta, survival.time, survival.status,utimes[j])
                                                    AUC[j] <- out$AUC
                                                  }
                                                  ## integrated AUC to get concordance measure
                                                  print(IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" ))
                                                  IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" )
                                                  
                                                  
                                                  
                                                } 
                                                roc.boot<-NULL
                                                roc.boot<-boot(datF3, boot.f2, r) 
                                                roc.boot 
                                                
                                                roc.boot_ALL[[j]] <- roc.boot  
                                                
                                              }
                                              roc.boot_ALL
                                              Sex<- roc.boot_ALL$Gender$t
                                              clinicl_stge<- roc.boot_ALL$clinicl_stge$t
                                              T_stage<- roc.boot_ALL$pT_stge$t
                                              N_stage<- roc.boot_ALL$pN_stage$t
                                              Age<- roc.boot_ALL$Age.y$t
                                              Tumor_type<- roc.boot_ALL$Tumor_type$t
                                              #Adjuvant_treatment<- roc.boot_ALL$T_Adjuvant$t
                                              SIA<- roc.boot_ALL$SIA_tri$t     
                                              Diff_grade<- roc.boot_ALL$Diff_grade$t
                                              median(SIA)
                                              
                                              
                                              Factprs_AUC <- cbind(Sex,Age,Diff_grade,T_stage,SIA)
                                              Factprs_AUC <- as.data.frame(Factprs_AUC)
                                              names(Factprs_AUC)[c(1:5)]<-c("Sex","Age","Diff_grade","T_stage","SIA")
                                              Factprs_AUC<-tibble::rownames_to_column(Factprs_AUC)
                                            head(Factprs_AUC)                                      
                                              
                                              All_clinical_together_and_SIA
                                              df_plot <- reshape2::melt(Factprs_AUC, id.vars="rowname")
                                              head(df_plot)
                                              
                                              df_plot$variable <- fct_rev(df_plot$variable)
                                              
                                              library(scales)
                                              library(RColorBrewer)
                                              
                                              
                                              k <- ggplot(df_plot, aes(x=reorder(variable, desc(variable)), y=value, fill=variable))+
                                                scale_fill_brewer(palette="RdBu")
                                              
                                              
                                              
                                              k+ 
                                                theme_bw()+ theme(panel.border = element_blank())+ theme(axis.line = element_line(colour = "black"))+
                                                # geom_point(aes(colour=variable,group = variable,y=value), size=2.3, alpha=1, position=position_dodge(width=0))+
                                                geom_boxplot(outlier.shape = NA, stat = "boxplot", width=0.7, alpha=1, position=position_dodge(0.9))+
                                                theme(legend.position = "none",
                                                      text = element_text(size=10),
                                                      axis.text.x=element_text(size=rel(2), angle=45))
                                              
                                              
                                              
                                              


                                         
                        