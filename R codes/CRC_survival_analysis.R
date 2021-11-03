########################
#
require(data.table)
library(metafor)
library(tidyr)
library(ggplot2) 
library(coin)
library(scales)
library(RColorBrewer)
library(data.table)
library(survival)
library(rms)
require(xlsx)
library(tidyverse)
library("survminer")

setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub")
dat <- NULL
dat <- read.table("Database_working_CRC.txt", sep="\t", header=T, fill = T)
dat$UCAN_TMA_Patnr <- as.character(dat$UCAN_TMA_Patnr)


            ###########################################################################
            ########### *** immune cells survival colon cancer st I-III *** ###########
            ###########################################################################



            data <- dat
            data <- data[ which(data$ColonRectum== 'Colon'), ]
            data <- data[ which( data$STAGE_AM == '1' | data$STAGE_AM == '2' | data$STAGE_AM == '3'), ]
            data <- data[ which(data$T_neoadjuvant != 'Yes'), ]
            dat <- data

            for(i in names(dat[,c( 201:211,231:238)] )){
              x  <-  as.data.frame(dat[!is.na(dat[i]),])
              
              quant_cut <- function(x, n) {
                qs <- quantile(x, 1:(n-1)/n)
                brks <- c(-Inf, qs, Inf)
                cut(x, breaks=brks, labels=FALSE)
              }
              
              x[i] <-quant_cut(x[[i]], 3)
              f <- "UCAN_TMA_Patnr"
              dat <- merge( dat,  x[,c(f,i)], by = f, all = T)
            }
            
            data <- dat
            
            data$OS_surgery_to_last_followup <- as.numeric(data$OS_surgery_to_last_followup)
            data$status <-(data$DEAD_ALIVE.1)
            data$status   <- ifelse((data$status =="Dead"), 1, 0)
            data$OS_surgery_to_last_followup <-  ifelse(( data$OS_surgery_to_last_followup >0), data$OS_surgery_to_last_followup, NA) 
            data <- data[!is.na(data$OS_surgery_to_last_followup), ] 
            data$SurvObj <- with(data, Surv(data$OS_surgery_to_last_followup, status == "1"))

            names(data) <-  gsub(".y", "", names(data), fixed = TRUE)

            
            colnames(data)
            CD4 <- coxph(SurvObj ~ CD4, data = data)
            CD4_activated <- coxph(SurvObj ~ CD4_activated, data = data)
            CD4_Treg <- coxph(SurvObj ~ CD4_Treg, data = data)
            CD8 <- coxph(SurvObj ~ CD8, data = data)
            CD8_activated <- coxph(SurvObj ~ CD8_activated, data = data)
            CD8_Treg <- coxph(SurvObj ~ CD8_Treg, data = data)
            B_cells <- coxph(SurvObj ~ B_cells, data = data)
            CD45RO <- coxph(SurvObj ~ CD45RO, data = data)
            NK <- coxph(SurvObj ~ NK, data = data)
            NKT <- coxph(SurvObj ~ NKT, data = data)
            CD68 <- coxph(SurvObj ~ CD68, data = data)
            CD163 <- coxph(SurvObj ~ CD163, data = data)
            M1 <- coxph(SurvObj ~ M1, data = data)
            Myeloid <- coxph(SurvObj ~ Myeloid, data = data)
            M2 <- coxph(SurvObj ~ M2, data = data)

            
            list<- list(CD4,
                        CD4_activated,
                        CD4_Treg,
                        CD8,
                        CD8_activated,
                        CD8_Treg,
                        B_cells,
                        CD45RO,
                        NK,
                        NKT,
                        CD68,
                        CD163,
                        M1,
                        Myeloid,
                        M2)
            
            celllist <-  NULL
            res<- NULL
            COX_results <-NULL
            for( i in unique(list) )
            {
              res<- NULL
              xxx<-as.data.frame(i$coefficients)
              name<-rownames(xxx)
              x <- summary(i)
              n<-x$n
              coeff<- as.data.table(x$coefficients)
              coeff<-coeff[1,]
              coeff$`Pr(>|z|)`
              p.value<-signif(coeff$`Pr(>|z|)`, digits=3)
              HR   <- coeff$`exp(coef)`
              conf.int<- as.data.table(x$conf.int)
              conf.int<-conf.int[1,]
              lower.95  <- conf.int$`lower .95`
              upper.95  <- conf.int$`upper .95`
              coef<-signif(x$coef[1], digits=3);#coeficient beta
              se.coef<-signif(x$coef[3], digits=3);#coeficient se(coef)
              res<-c(name,HR, lower.95, upper.95,p.value, n, coef, se.coef)
              names(res)<-c("cell","HR", "lower.95", "upper.95", "p.value", "n", "coef", "se.coef")

              COX_results <- data.frame(rbind(COX_results, res))
            }
            
            write.csv(x=COX_results, file="COX_results.csv")
            Cox.df <-read_csv("COX_results.csv")
            names <- as.vector(COX_results$cell)
            labs <- names
            yi   <- as.numeric(COX_results$coef)
            sei  <- as.numeric(COX_results$se.coef)
            p.vals  <- as.numeric(COX_results$p.value)
            
            # Combine data into summary estimate
            res  <- rma(yi=yi, sei=sei, method="FE")
            summary(res)
            
            # Format pvalues so only those bellow 0.01 are scientifically notated
            p.vals<- as.numeric(p.vals)
            p.vals <- ifelse(p.vals < 0.001, 
                             format(p.vals,digits = 3,scientific = TRUE,trim = TRUE),
                             format(round(p.vals, 3), nsmall=2, trim=TRUE))
            
            p.vals <- gsub('e(.*)', ' x 10^\\1', p.vals)
            
            ### ferestplot:

            forest(res, transf=exp, refline=1, xlab="HR (95%CI)", 
                   slab=labs, ilab = p.vals, ilab.xpos = 2, mlab="Summary Estimate", alim=c(0,2), xlim=c(-1,4),steps=5, cex=1)
            
            
            ### another ferestplot:
            
            library(forestplot)         

                          data_PLOT <- structure(list(HR  = c(COX_results$HR), 
                                                      low = c(COX_results$lower.95),
                                                      high = c(COX_results$upper.95)),
                                                 .Names = c("HR", "low", "high"), 
                                                 row.names = c(COX_results$cell), 
                                                 class = "data.frame")
                          
                          print("....Creating the plot....")                 
                          jpeg(filename="Hazard_ratio_plot.jpg",units="cm",width=20,height=17, res=800)
                          
                          dev.off()
                          forestplot(data_PLOT,
                                     
                                     mean=c(data_PLOT$HR), 
                                     lower=as.numeric(c(data_PLOT$low)), upper=as.numeric(c(data_PLOT$high)),
                                     graph.pos = 2,
                                     zero      = 1,
                                     #xlog      = TRUE,
                                     boxsize=0.2,
                                     col = fpColors(box="red",line="darkblue"),
                                     ci.vertices = F,
                                     
                                     title="Hazard Ratio",
                                     lineheight = "auto",
                                     clip = c(-Inf,500),
                                     
                                     txt_gp=fpTxtGp(label=gpar(cex=0.5),
                                                    ticks=gpar(cex=1.1),
                                                    xlab=gpar(cex = 0.2),
                                                    title=gpar(cex = 0.2)))
            
            

                          
                          
                          

###########################################################################
################ *** SIA and IS colon cancer st I-III *** #################
################ ***               OS                 *** #################
###########################################################################



setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub")
dat <- NULL
dat <- read.table("Database_working_CRC.txt", sep="\t", header=T, fill = T)
dat$UCAN_TMA_Patnr <- as.character(dat$UCAN_TMA_Patnr)
dat$SIA <- dat$CD8/(dat$CD8+dat$M2)

data <- dat
data <- data[ which(data$ColonRectum== 'Colon'), ]
data <- data[ which( data$STAGE_AM == '1' | data$STAGE_AM == '2' | data$STAGE_AM == '3'), ]
data <- data[ which(data$T_neoadjuvant != 'Yes'), ]
dat <- data
colnames(dat)     

x  <-  dat[,c(1,256)]
x  <-  na.omit(x)
head(x)

quant_cut <- function(x, n) {
  qs <- quantile(x, 1:(n-1)/n)
  brks <- c(-Inf, qs, Inf)
  cut(x, breaks=brks, labels=FALSE)
}

x$SIA_tri <-quant_cut(x$SIA, 3)
x<- x[-2]
dat <- merge(dat, x,by="UCAN_TMA_Patnr")
data <- dat

data$OS_surgery_to_last_followup <- as.numeric(data$OS_surgery_to_last_followup)
data$status <-(data$DEAD_ALIVE.1)
data$status   <- ifelse((data$status =="Dead"), 1, 0)
data$OS_surgery_to_last_followup <-  ifelse(( data$OS_surgery_to_last_followup >0), data$OS_surgery_to_last_followup, NA) 
data <- data[!is.na(data$OS_surgery_to_last_followup), ] 
data$TIME <-data$OS_surgery_to_last_followup
data$SurvObj <- with(data, Surv(data$OS_surgery_to_last_followup, status == "1"))


# Kaplan-Meier for SIA OS

datF <- data
datF <-datF[ which(!is.na(datF$SIA_tri)), ]
datF$SIA_tri <- as.factor(datF$SIA_tri)
datF$SIA_tri <- relevel(datF$SIA_tri, ref = "1")

fit <- survfit(SurvObj ~ SIA_tri, data = datF, conf.type = "log-log")

c<-coxph(SurvObj ~ SIA_tri, data = datF)
c
ggsurvplot(fit, risk.table = T, lwd=1,pval = TRUE)

#####

# Kaplan-Meier for IS OS

par(mfrow=c(3,6))

datF <- data
datF <-datF[ which(!is.na(datF$ISLike)), ]
datF$ISLike <- as.factor(datF$ISLike)
datF$ISLike <- relevel(datF$ISLike, ref = "1")

fit <- survfit(SurvObj ~ ISLike, data = datF, conf.type = "log-log")
c<-coxph(SurvObj ~ ISLike, data = datF)
c
ggsurvplot(fit, risk.table = T, lwd=1,pval = TRUE)


#############         multivariable Cox
datF <- data
datF <-datF[ which(!is.na(datF$SIA_tri)), ]
data_multi <- datF


data_multi$T_stage <- as.character(data_multi$T_stage)
data_multi$T_stage <- as.factor(data_multi$T_stage)
data_multi$T_stage <- relevel(data_multi$T_stage, ref = "1")

data_multi <- data_multi[ which(data_multi$N_stage != 'missing'), ]
data_multi$N_stage <- factor(data_multi$N_stage, levels = c("0", "1"))
data_multi$N_stage <- relevel(data_multi$N_stage, ref = "0")

data_multi$DIFFGRAD.1 = factor(data_multi$DIFFGRAD.1, levels = c("low grade", "missing", "high grade"))
data_multi$DIFFGRAD.1 = relevel(data_multi$DIFFGRAD.1, ref = "low grade")    

data_multi$T_Adjuvant = factor(data_multi$T_Adjuvant, levels = c("No","missing","Yes"))
data_multi$T_Adjuvant = relevel(data_multi$T_Adjuvant, ref = "No")   

data_multi$AgeDiagnosis <- ifelse(data_multi$AgeDiagnosis>75, 1,0)

data_multi$Left_Right <- factor(data_multi$Left_Right, levels = c('left', 'right'))

data_multi$MSS_T = factor(data_multi$MSS_T, levels = c("MSS","missing","MSI"))
data_multi$MSS_T = relevel(data_multi$MSS_T, ref = "MSS")  



data_multi$ISLike <- as.factor(data_multi$ISLike)
data_multi$ISLike <- relevel(data_multi$ISLike, ref = "1")

data_multi$SIA_tri <- as.factor(data_multi$SIA_tri)
data_multi$SIA_tri <- relevel(data_multi$SIA_tri, ref = "1")


multi.cox<-coxph(SurvObj ~ SIA_tri + ISLike+
                   T_stage +N_stage+
                   AgeDiagnosis +Kön+MSS_T, data = data_multi)


multi.cox
coefficients <-summary(multi.cox)$coefficient

Cox2 <-as.data.frame(coefficients)

write.csv(x=Cox2, file="Multi_Cox.csv")





###########################################
################    AUCs   ################
###########################################

library(pROC)
library(boot)
library(MASS)
library(risksetROC)

###                AUC             ###
#### all covariables separatelly  ####
###                                ###


i <- NULL
datF2 <- data_multi

factors<-colnames(datF2[,c("Kön","MSS_T", "Left_Right" ,"DIFFGRAD.1", "T_stage", "N_stage","ISLike", "SIA_tri")])
roc.boot_ALL <- NULL
for(j in unique(factors)){
  
  datF3<-NULL
  datF3 <- datF2[j]
  datF3$test <- datF2$TIME
  datF3$dx <-datF2$status
  datF3 <- na.omit(datF3)
  # variables
  names(datF3)[1] <- "score"
  datF3$test <- as.numeric(datF3$test)
  datF3$dx <- as.numeric(datF3$dx)
  datF3$score <- as.factor(datF3$score)
  
  r <- 1000
  boot.f2 <- function(d, i){
    data <- d[i,]
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
    for( jj in 1:length(utimes) )
    {
      out <- CoxWeights( eta, survival.time, survival.status,utimes[jj])
      AUC[jj] <- out$AUC
    }
    print(IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" ))
    IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" )
    
  } 
  roc.boot<-NULL
  roc.boot<-boot(datF3, boot.f2, r) 
  roc.boot_ALL[[j]] <- roc.boot  
  
}

Sex<- roc.boot_ALL$Kön$t
MSI<- roc.boot_ALL$MSS_T$t
sideness<- roc.boot_ALL$Left_Right$t
Differentiation<- roc.boot_ALL$DIFFGRAD.1$t
T_stage<- roc.boot_ALL$T_stage$t
N_stage<- roc.boot_ALL$N_stage$t
ISLike<- roc.boot_ALL$ISLike$t 
SIA<- roc.boot_ALL$SIA_tri$t     


###                AUC             ###
####     all clinical together    ####
###                                ###

i <- NULL
datF2 <- data_multi

datF3<-NULL
datF3 <- datF2
datF3$test <- datF2$TIME
datF3$dx <-datF2$status
datF3$test <- as.numeric(datF3$test)
datF3$dx <- as.numeric(datF3$dx)
datF3$Sex<- datF3$Kön
datF3$MSI<- datF3$MSS_T
datF3$sideness<- datF3$Left_Right
datF3$Differentiation<- datF3$DIFFGRAD.1
datF3$T_stage<- datF3$T_stage
datF3$N_stage<- datF3$N_stage
datF3$ISLike<- datF3$ISLike
datF3$SIA<- datF3$SIA_tri  

r <- 1000
boot.f2 <- function(d, i){
  data <- d[i,]
  print(i)
  survival.time <- data$test
  score <- data$score
  survival.status <- data$dx
  surv.prob <- unique(survfit(Surv(survival.time,survival.status)~1)$surv)
  fit0 <- coxph( Surv(survival.time,survival.status)
                 ~data$Sex+
                   data$MSI+
                   data$sideness+
                   data$Differentiation+
                   data$T_stage+
                   data$N_stage, na.action=na.omit )
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

ALL_clinical <- roc.boot$t



###                AUC             ###
#### all clinical AND IS together ####
###                                ###

i <- NULL

datF2 <- data_multi
datF3<-NULL
datF3 <- datF2
datF3$test <- datF2$TIME
datF3$dx <-datF2$status
datF3$test <- as.numeric(datF3$test)
datF3$dx <- as.numeric(datF3$dx)
datF3$Sex<- datF3$Kön
datF3$MSI<- datF3$MSS_T
datF3$sideness<- datF3$Left_Right
datF3$Differentiation<- datF3$DIFFGRAD.1
datF3$T_stage<- datF3$T_stage
datF3$N_stage<- datF3$N_stage
datF3$ISLike<- datF3$ISLike
datF3$SIA<- datF3$SIA_tri  

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
                 ~data$Sex+
                   data$MSI+
                   data$sideness+
                   data$Differentiation+
                   data$T_stage+
                   data$N_stage+
                   data$ISLike, na.action=na.omit )
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

ALL_clinical_and_IS <- roc.boot$t


###                AUC             ###
#### all clinical AND SIA together ####
###                                ###

i <- NULL

datF2 <- data_multi
datF3<-NULL
datF3 <- datF2
datF3$test <- datF2$TIME
datF3$dx <-datF2$status
datF3$test <- as.numeric(datF3$test)
datF3$dx <- as.numeric(datF3$dx)
datF3$Sex<- datF3$Kön
datF3$MSI<- datF3$MSS_T
datF3$sideness<- datF3$Left_Right
datF3$Differentiation<- datF3$DIFFGRAD.1
datF3$T_stage<- datF3$T_stage
datF3$N_stage<- datF3$N_stage
datF3$ISLike<- datF3$ISLike
datF3$SIA<- datF3$SIA_tri  

r <- 1000
boot.f2 <- function(d, i){
  data <- d[i,]
  print(i)
  survival.time <- data$test
  score <- data$score
  survival.status <- data$dx
  surv.prob <- unique(survfit(Surv(survival.time,survival.status)~1)$surv)
  fit0 <- coxph( Surv(survival.time,survival.status)
                 ~data$Sex+
                   data$MSI+
                   data$sideness+
                   data$Differentiation+
                   data$T_stage+
                   data$N_stage+
                   data$SIA, na.action=na.omit )
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

ALL_clinical_and_SIA <- roc.boot$t


###                AUC                   ###
#### all clinical, IS and SIA together  ####
###                                      ###

i <- NULL

datF2 <- data_multi
datF3<-NULL
datF3 <- datF2
datF3$test <- datF2$TIME
datF3$dx <-datF2$status
datF3$test <- as.numeric(datF3$test)
datF3$dx <- as.numeric(datF3$dx)
datF3$Sex<- datF3$Kön
datF3$MSI<- datF3$MSS_T
datF3$sideness<- datF3$Left_Right
datF3$Differentiation<- datF3$DIFFGRAD.1
datF3$T_stage<- datF3$T_stage
datF3$N_stage<- datF3$N_stage
datF3$ISLike<- datF3$ISLike
datF3$SIA<- datF3$SIA_tri  

r <- 1000
boot.f2 <- function(d, i){
  data <- d[i,]
  print(i)
  survival.time <- data$test
  score <- data$score
  survival.status <- data$dx
  surv.prob <- unique(survfit(Surv(survival.time,survival.status)~1)$surv)
  fit0 <- coxph( Surv(survival.time,survival.status)
                 ~data$Sex+
                   data$MSI+
                   data$sideness+
                   data$Differentiation+
                   data$T_stage+
                   data$N_stage+
                   data$ISLike+
                   data$SIA, na.action=na.omit )
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

ALL_clinical_IS_and_SIA <- roc.boot$t


#####
#####       PLOT AUCs with all combinations
#####


Factprs_AUC <- cbind(Sex,MSI,sideness, Differentiation,T_stage, N_stage,ISLike, SIA, ALL_clinical,ALL_clinical_and_IS, ALL_clinical_and_SIA, ALL_clinical_IS_and_SIA)
Factprs_AUC <- as.data.frame(Factprs_AUC)
names(Factprs_AUC)[c(1:12)]<-c("Sex","MSI","Sideness", "Differentiation","T stage", "N stage","ISLike", "SIA","All clinical", "ALL_clinical_and_IS", "All clinical + SIA", "ALL_clinical_IS_and_SIA")
Factprs_AUC<-tibble::rownames_to_column(Factprs_AUC)

summary<- summary(Factprs_AUC) 
write.table(summary, file ="AUCs summary.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)

df_plot <- reshape2::melt(Factprs_AUC, id.vars="rowname")
df_plot$variable <- fct_rev(df_plot$variable)
unique(df_plot$variable)

k <- ggplot(df_plot, aes(x=reorder(variable, desc(variable)), y=value, fill=variable))+
  scale_fill_brewer(palette="RdBu")

k+ 
  theme_bw()+ theme(panel.border = element_blank())+ theme(axis.line = element_line(colour = "black"))+
  # geom_point(aes(colour=variable,group = variable,y=value), size=2.3, alpha=1, position=position_dodge(width=0))+
  geom_boxplot(outlier.shape = NA, stat = "boxplot", width=0.7, alpha=1, position=position_dodge(0.9))+
  theme(legend.position = "none",
        text = element_text(size=10),
        axis.text.x=element_text(size=rel(2), angle=45))+ 
  ylim(0.49, 0.77)


############# statistics

ALLclinical <- as.data.frame(ALL_clinical)
ALLclinical_and_SIA <- as.data.frame(ALL_clinical_and_SIA)
ALLclinical$group <- "ALLclinical"
ALLclinical_and_SIA$group <- "ALLclinical_and_SIA"
All_and_SIA <- rbind(ALLclinical, ALLclinical_and_SIA)

wilcox.test(V1 ~ group, All_and_SIA,
            exact = FALSE, paired=F)


ALL_clinical_IS_and_SIA <- as.data.frame(ALL_clinical_IS_and_SIA)
ALL_clinical_IS_and_SIA$group <- "ALL_clinical_IS_and_SIA"
All_and_SIA <- rbind(ALLclinical_and_SIA, ALL_clinical_IS_and_SIA)

wilcox.test(V1 ~ group, All_and_SIA,
            exact = FALSE, paired=F)



SIA <- as.data.frame(SIA)
SIA$group <- "SIA"                  

N_stage<- as.data.frame(T_stage)
N_stage$group <- "N_stage"

Nstage_and_SIA <- rbind(N_stage, SIA)
wilcox.test(V1 ~ group, Nstage_and_SIA,
            exact = FALSE, paired=F)


#########
#########         relative contribution      
#########

data_multi <- data
data_multi <-data_multi[ which(!is.na(data_multi$SIA_tri)), ]

data_multi$T_stage <- as.character(data_multi$T_stage)

data_multi$T_stage <- as.factor(data_multi$T_stage)
data_multi$T_stage <- relevel(data_multi$T_stage, ref = "1")

data_multi <- data_multi[ which(data_multi$N_stage != 'missing'), ]
data_multi$N_stage <- factor(data_multi$N_stage, levels = c("0", "1"))
data_multi$N_stage <- relevel(data_multi$N_stage, ref = "0")

data_multi$DIFFGRAD.1 = factor(data_multi$DIFFGRAD.1, levels = c("low grade", "missing", "high grade"))
data_multi$DIFFGRAD.1 = relevel(data_multi$DIFFGRAD.1, ref = "low grade")    

data_multi$T_Adjuvant = factor(data_multi$T_Adjuvant, levels = c("No","missing","Yes"))
data_multi$T_Adjuvant = relevel(data_multi$T_Adjuvant, ref = "No")   

data_multi$AgeDiagnosis <- ifelse(data_multi$AgeDiagnosis>75, 1,0)

data_multi$Left_Right <- factor(data_multi$Left_Right, levels = c('left', 'right'))


data_multi$MSS_T = factor(data_multi$MSS_T, levels = c("MSS","missing","MSI"))
data_multi$MSS_T = relevel(data_multi$MSS_T, ref = "MSS")  


colnames(data)
#dev.off()
par(mfrow=c(2,2))

data_multi$ISLike <- as.factor(data_multi$ISLike)
data_multi$ISLike <- relevel(data_multi$ISLike, ref = "1")

data_multi$SIA_tri <- as.factor(data_multi$SIA_tri)
data_multi$SIA_tri <- relevel(data_multi$SIA_tri, ref = "1")



####  OS withought IS
f <- cph(SurvObj ~ SIA_tri + 
           T_stage +N_stage+DIFFGRAD.1+Left_Right+Kön,
         x=TRUE, y=TRUE, data=data_multi)

CHI <- anova(f)
n<-dim(CHI )[1]
CHI <- as.data.frame(CHI)
setDT(CHI, keep.rownames = TRUE)[]
CHI <-CHI [1:(n-1),]

write.table(CHI, file ="CHI OS without IS .txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
#Draw pie
pie(CHI$`Chi-Square`, labels = CHI$rn)



####  OS withIS
f <- cph(SurvObj ~ SIA_tri + ISLike+
           T_stage +N_stage+DIFFGRAD.1+Left_Right+Kön,
         x=TRUE, y=TRUE, data=data_multi)

CHI <- anova(f)
n<-dim(CHI )[1]
CHI <- as.data.frame(CHI)
setDT(CHI, keep.rownames = TRUE)[]
CHI <-CHI [1:(n-1),]

write.table(CHI, file ="CHI OS with IS .txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
#Draw pie
pie(CHI$`Chi-Square`, labels = CHI$rn)
























###########################################################################
################ *** SIA and IS colon cancer st I-III *** #################
################ ***               RFS                *** #################
###########################################################################



setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub")
dat <- NULL
dat <- read.table("Database_working_CRC.txt", sep="\t", header=T, fill = T)
dat$UCAN_TMA_Patnr <- as.character(dat$UCAN_TMA_Patnr)
dat$SIA <- dat$CD8/(dat$CD8+dat$M2)

data <- dat
data <- data[ which(data$ColonRectum== 'Colon'), ]
data <- data[ which( data$STAGE_AM == '1' | data$STAGE_AM == '2' | data$STAGE_AM == '3'), ]
data <- data[ which(data$T_neoadjuvant != 'Yes'), ]
dat <- data
colnames(dat)     

x  <-  dat[,c(1,256)]
x  <-  na.omit(x)
head(x)

quant_cut <- function(x, n) {
  qs <- quantile(x, 1:(n-1)/n)
  brks <- c(-Inf, qs, Inf)
  cut(x, breaks=brks, labels=FALSE)
}

x$SIA_tri <-quant_cut(x$SIA, 3)
x<- x[-2]
dat <- merge(dat, x,by="UCAN_TMA_Patnr")
data <- dat

data$DFS_surgery_to_last_followup <- as.numeric(data$DFS_surgery_to_last_followup)
data$DFS_surgery_to_last_followup<-  ifelse(( data$DFS_surgery_to_last_followup >0), data$DFS_surgery_to_last_followup, NA) 
data <- data[!is.na(data$DFS_surgery_to_last_followup), ] 
data$status <-(data$DFS_event)
data$status   <- ifelse((data$status =="Yes"), 1, 0)
data$TIME <-data$DFS_surgery_to_last_followup
data$SurvObj <- with(data, Surv(data$DFS_surgery_to_last_followup, status == "1"))


# Kaplan-Meier for SIA OS

datF <- data
datF <-datF[ which(!is.na(datF$SIA_tri)), ]
datF$SIA_tri <- as.factor(datF$SIA_tri)
datF$SIA_tri <- relevel(datF$SIA_tri, ref = "1")

fit <- survfit(SurvObj ~ SIA_tri, data = datF, conf.type = "log-log")

c<-coxph(SurvObj ~ SIA_tri, data = datF)
c
ggsurvplot(fit, risk.table = T, lwd=1,pval = TRUE)

#####

# Kaplan-Meier for IS OS

par(mfrow=c(3,6))

datF <- data
datF <-datF[ which(!is.na(datF$ISLike)), ]
datF$ISLike <- as.factor(datF$ISLike)
datF$ISLike <- relevel(datF$ISLike, ref = "1")

fit <- survfit(SurvObj ~ ISLike, data = datF, conf.type = "log-log")
c<-coxph(SurvObj ~ ISLike, data = datF)
c
ggsurvplot(fit, risk.table = T, lwd=1,pval = TRUE)


#############         multivariable Cox
datF <- data
datF <-datF[ which(!is.na(datF$SIA_tri)), ]
data_multi <- datF


data_multi$T_stage <- as.character(data_multi$T_stage)
data_multi$T_stage <- as.factor(data_multi$T_stage)
data_multi$T_stage <- relevel(data_multi$T_stage, ref = "1")

data_multi <- data_multi[ which(data_multi$N_stage != 'missing'), ]
data_multi$N_stage <- factor(data_multi$N_stage, levels = c("0", "1"))
data_multi$N_stage <- relevel(data_multi$N_stage, ref = "0")

data_multi$DIFFGRAD.1 = factor(data_multi$DIFFGRAD.1, levels = c("low grade", "missing", "high grade"))
data_multi$DIFFGRAD.1 = relevel(data_multi$DIFFGRAD.1, ref = "low grade")    

data_multi$T_Adjuvant = factor(data_multi$T_Adjuvant, levels = c("No","missing","Yes"))
data_multi$T_Adjuvant = relevel(data_multi$T_Adjuvant, ref = "No")   

data_multi$AgeDiagnosis <- ifelse(data_multi$AgeDiagnosis>75, 1,0)

data_multi$Left_Right <- factor(data_multi$Left_Right, levels = c('left', 'right'))

data_multi$MSS_T = factor(data_multi$MSS_T, levels = c("MSS","missing","MSI"))
data_multi$MSS_T = relevel(data_multi$MSS_T, ref = "MSS")  



data_multi$ISLike <- as.factor(data_multi$ISLike)
data_multi$ISLike <- relevel(data_multi$ISLike, ref = "1")

data_multi$SIA_tri <- as.factor(data_multi$SIA_tri)
data_multi$SIA_tri <- relevel(data_multi$SIA_tri, ref = "1")


multi.cox<-coxph(SurvObj ~ SIA_tri + ISLike+
                   T_stage +N_stage+
                   AgeDiagnosis +Kön+MSS_T, data = data_multi)


multi.cox
coefficients <-summary(multi.cox)$coefficient

Cox2 <-as.data.frame(coefficients)

write.csv(x=Cox2, file="Multi_Cox.csv")





###########################################
################    AUCs   ################
###########################################

library(pROC)
library(boot)
library(MASS)
library(risksetROC)

###                AUC             ###
#### all covariables separatelly  ####
###                                ###


i <- NULL
datF2 <- data_multi

factors<-colnames(datF2[,c("Kön","MSS_T", "Left_Right" ,"DIFFGRAD.1", "T_stage", "N_stage","ISLike", "SIA_tri")])
roc.boot_ALL <- NULL
for(j in unique(factors)){
  
  datF3<-NULL
  datF3 <- datF2[j]
  datF3$test <- datF2$TIME
  datF3$dx <-datF2$status
  datF3 <- na.omit(datF3)
  # variables
  names(datF3)[1] <- "score"
  datF3$test <- as.numeric(datF3$test)
  datF3$dx <- as.numeric(datF3$dx)
  datF3$score <- as.factor(datF3$score)
  
  r <- 1000
  boot.f2 <- function(d, i){
    data <- d[i,]
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
    for( jj in 1:length(utimes) )
    {
      out <- CoxWeights( eta, survival.time, survival.status,utimes[jj])
      AUC[jj] <- out$AUC
    }
    print(IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" ))
    IntegrateAUC( AUC, utimes, surv.prob, tmax=1000, weight="rescale" )
    
  } 
  roc.boot<-NULL
  roc.boot<-boot(datF3, boot.f2, r) 
  roc.boot_ALL[[j]] <- roc.boot  
  
}

Sex<- roc.boot_ALL$Kön$t
MSI<- roc.boot_ALL$MSS_T$t
sideness<- roc.boot_ALL$Left_Right$t
Differentiation<- roc.boot_ALL$DIFFGRAD.1$t
T_stage<- roc.boot_ALL$T_stage$t
N_stage<- roc.boot_ALL$N_stage$t
ISLike<- roc.boot_ALL$ISLike$t 
SIA<- roc.boot_ALL$SIA_tri$t     


###                AUC             ###
####     all clinical together    ####
###                                ###

i <- NULL
datF2 <- data_multi

datF3<-NULL
datF3 <- datF2
datF3$test <- datF2$TIME
datF3$dx <-datF2$status
datF3$test <- as.numeric(datF3$test)
datF3$dx <- as.numeric(datF3$dx)
datF3$Sex<- datF3$Kön
datF3$MSI<- datF3$MSS_T
datF3$sideness<- datF3$Left_Right
datF3$Differentiation<- datF3$DIFFGRAD.1
datF3$T_stage<- datF3$T_stage
datF3$N_stage<- datF3$N_stage
datF3$ISLike<- datF3$ISLike
datF3$SIA<- datF3$SIA_tri  

r <- 1000
boot.f2 <- function(d, i){
  data <- d[i,]
  print(i)
  survival.time <- data$test
  score <- data$score
  survival.status <- data$dx
  surv.prob <- unique(survfit(Surv(survival.time,survival.status)~1)$surv)
  fit0 <- coxph( Surv(survival.time,survival.status)
                 ~data$Sex+
                   data$MSI+
                   data$sideness+
                   data$Differentiation+
                   data$T_stage+
                   data$N_stage, na.action=na.omit )
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

ALL_clinical <- roc.boot$t



###                AUC             ###
#### all clinical AND IS together ####
###                                ###

i <- NULL

datF2 <- data_multi
datF3<-NULL
datF3 <- datF2
datF3$test <- datF2$TIME
datF3$dx <-datF2$status
datF3$test <- as.numeric(datF3$test)
datF3$dx <- as.numeric(datF3$dx)
datF3$Sex<- datF3$Kön
datF3$MSI<- datF3$MSS_T
datF3$sideness<- datF3$Left_Right
datF3$Differentiation<- datF3$DIFFGRAD.1
datF3$T_stage<- datF3$T_stage
datF3$N_stage<- datF3$N_stage
datF3$ISLike<- datF3$ISLike
datF3$SIA<- datF3$SIA_tri  

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
                 ~data$Sex+
                   data$MSI+
                   data$sideness+
                   data$Differentiation+
                   data$T_stage+
                   data$N_stage+
                   data$ISLike, na.action=na.omit )
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

ALL_clinical_and_IS <- roc.boot$t


###                AUC             ###
#### all clinical AND SIA together ####
###                                ###

i <- NULL

datF2 <- data_multi
datF3<-NULL
datF3 <- datF2
datF3$test <- datF2$TIME
datF3$dx <-datF2$status
datF3$test <- as.numeric(datF3$test)
datF3$dx <- as.numeric(datF3$dx)
datF3$Sex<- datF3$Kön
datF3$MSI<- datF3$MSS_T
datF3$sideness<- datF3$Left_Right
datF3$Differentiation<- datF3$DIFFGRAD.1
datF3$T_stage<- datF3$T_stage
datF3$N_stage<- datF3$N_stage
datF3$ISLike<- datF3$ISLike
datF3$SIA<- datF3$SIA_tri  

r <- 1000
boot.f2 <- function(d, i){
  data <- d[i,]
  print(i)
  survival.time <- data$test
  score <- data$score
  survival.status <- data$dx
  surv.prob <- unique(survfit(Surv(survival.time,survival.status)~1)$surv)
  fit0 <- coxph( Surv(survival.time,survival.status)
                 ~data$Sex+
                   data$MSI+
                   data$sideness+
                   data$Differentiation+
                   data$T_stage+
                   data$N_stage+
                   data$SIA, na.action=na.omit )
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

ALL_clinical_and_SIA <- roc.boot$t


###                AUC                   ###
#### all clinical, IS and SIA together  ####
###                                      ###

i <- NULL

datF2 <- data_multi
datF3<-NULL
datF3 <- datF2
datF3$test <- datF2$TIME
datF3$dx <-datF2$status
datF3$test <- as.numeric(datF3$test)
datF3$dx <- as.numeric(datF3$dx)
datF3$Sex<- datF3$Kön
datF3$MSI<- datF3$MSS_T
datF3$sideness<- datF3$Left_Right
datF3$Differentiation<- datF3$DIFFGRAD.1
datF3$T_stage<- datF3$T_stage
datF3$N_stage<- datF3$N_stage
datF3$ISLike<- datF3$ISLike
datF3$SIA<- datF3$SIA_tri  

r <- 1000
boot.f2 <- function(d, i){
  data <- d[i,]
  print(i)
  survival.time <- data$test
  score <- data$score
  survival.status <- data$dx
  surv.prob <- unique(survfit(Surv(survival.time,survival.status)~1)$surv)
  fit0 <- coxph( Surv(survival.time,survival.status)
                 ~data$Sex+
                   data$MSI+
                   data$sideness+
                   data$Differentiation+
                   data$T_stage+
                   data$N_stage+
                   data$ISLike+
                   data$SIA, na.action=na.omit )
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

ALL_clinical_IS_and_SIA <- roc.boot$t


#####
#####       PLOT AUCs with all combinations
#####


Factprs_AUC <- cbind(Sex,MSI,sideness, Differentiation,T_stage, N_stage,ISLike, SIA, ALL_clinical,ALL_clinical_and_IS, ALL_clinical_and_SIA, ALL_clinical_IS_and_SIA)
Factprs_AUC <- as.data.frame(Factprs_AUC)
names(Factprs_AUC)[c(1:12)]<-c("Sex","MSI","Sideness", "Differentiation","T stage", "N stage","ISLike", "SIA","All clinical", "ALL_clinical_and_IS", "All clinical + SIA", "ALL_clinical_IS_and_SIA")
Factprs_AUC<-tibble::rownames_to_column(Factprs_AUC)

summary<- summary(Factprs_AUC) 
write.table(summary, file ="AUCs summary.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)

df_plot <- reshape2::melt(Factprs_AUC, id.vars="rowname")
df_plot$variable <- fct_rev(df_plot$variable)
unique(df_plot$variable)

k <- ggplot(df_plot, aes(x=reorder(variable, desc(variable)), y=value, fill=variable))+
  scale_fill_brewer(palette="RdBu")

k+ 
  theme_bw()+ theme(panel.border = element_blank())+ theme(axis.line = element_line(colour = "black"))+
  # geom_point(aes(colour=variable,group = variable,y=value), size=2.3, alpha=1, position=position_dodge(width=0))+
  geom_boxplot(outlier.shape = NA, stat = "boxplot", width=0.7, alpha=1, position=position_dodge(0.9))+
  theme(legend.position = "none",
        text = element_text(size=10),
        axis.text.x=element_text(size=rel(2), angle=45))+ 
  ylim(0.49, 0.77)


############# statistics

ALLclinical <- as.data.frame(ALL_clinical)
ALLclinical_and_SIA <- as.data.frame(ALL_clinical_and_SIA)
ALLclinical$group <- "ALLclinical"
ALLclinical_and_SIA$group <- "ALLclinical_and_SIA"
All_and_SIA <- rbind(ALLclinical, ALLclinical_and_SIA)

wilcox.test(V1 ~ group, All_and_SIA,
            exact = FALSE, paired=F)


ALL_clinical_IS_and_SIA <- as.data.frame(ALL_clinical_IS_and_SIA)
ALL_clinical_IS_and_SIA$group <- "ALL_clinical_IS_and_SIA"
All_and_SIA <- rbind(ALLclinical_and_SIA, ALL_clinical_IS_and_SIA)

wilcox.test(V1 ~ group, All_and_SIA,
            exact = FALSE, paired=F)



SIA <- as.data.frame(SIA)
SIA$group <- "SIA"                  

N_stage<- as.data.frame(T_stage)
N_stage$group <- "N_stage"

Nstage_and_SIA <- rbind(N_stage, SIA)
wilcox.test(V1 ~ group, Nstage_and_SIA,
            exact = FALSE, paired=F)


#########
#########         relative contribution      
#########

data_multi <- data
data_multi <-data_multi[ which(!is.na(data_multi$SIA_tri)), ]

data_multi$T_stage <- as.character(data_multi$T_stage)

data_multi$T_stage <- as.factor(data_multi$T_stage)
data_multi$T_stage <- relevel(data_multi$T_stage, ref = "1")

data_multi <- data_multi[ which(data_multi$N_stage != 'missing'), ]
data_multi$N_stage <- factor(data_multi$N_stage, levels = c("0", "1"))
data_multi$N_stage <- relevel(data_multi$N_stage, ref = "0")

data_multi$DIFFGRAD.1 = factor(data_multi$DIFFGRAD.1, levels = c("low grade", "missing", "high grade"))
data_multi$DIFFGRAD.1 = relevel(data_multi$DIFFGRAD.1, ref = "low grade")    

data_multi$T_Adjuvant = factor(data_multi$T_Adjuvant, levels = c("No","missing","Yes"))
data_multi$T_Adjuvant = relevel(data_multi$T_Adjuvant, ref = "No")   

data_multi$AgeDiagnosis <- ifelse(data_multi$AgeDiagnosis>75, 1,0)

data_multi$Left_Right <- factor(data_multi$Left_Right, levels = c('left', 'right'))


data_multi$MSS_T = factor(data_multi$MSS_T, levels = c("MSS","missing","MSI"))
data_multi$MSS_T = relevel(data_multi$MSS_T, ref = "MSS")  


colnames(data)
#dev.off()
par(mfrow=c(2,2))

data_multi$ISLike <- as.factor(data_multi$ISLike)
data_multi$ISLike <- relevel(data_multi$ISLike, ref = "1")

data_multi$SIA_tri <- as.factor(data_multi$SIA_tri)
data_multi$SIA_tri <- relevel(data_multi$SIA_tri, ref = "1")



####  OS withought IS
f <- cph(SurvObj ~ SIA_tri +
           T_stage +N_stage+DIFFGRAD.1+Left_Right+Kön,
         x=TRUE, y=TRUE, data=data_multi)

CHI <- anova(f)
n<-dim(CHI )[1]
CHI <- as.data.frame(CHI)
setDT(CHI, keep.rownames = TRUE)[]
CHI <-CHI [1:(n-1),]

write.table(CHI, file ="CHI without IS .txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
#Draw pie
pie(CHI$`Chi-Square`, labels = CHI$rn)



####  OS withIS
f <- cph(SurvObj ~ SIA_tri + ISLike+
           T_stage +N_stage+DIFFGRAD.1+Left_Right+Kön,
         x=TRUE, y=TRUE, data=data_multi)

CHI <- anova(f)
n<-dim(CHI )[1]
CHI <- as.data.frame(CHI)
setDT(CHI, keep.rownames = TRUE)[]
CHI <-CHI [1:(n-1),]

write.table(CHI, file ="CHI with IS .txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)
#Draw pie
pie(CHI$`Chi-Square`, labels = CHI$rn)














###########################################################################
################ ***    SIA in colon cancer st II     *** #################
################ ***                                  *** #################
###########################################################################


setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub")
dat <- NULL
dat <- read.table("Database_working_CRC.txt", sep="\t", header=T, fill = T)
dat$UCAN_TMA_Patnr <- as.character(dat$UCAN_TMA_Patnr)
dat$SIA <- dat$CD8/(dat$CD8+dat$M2)

data <- dat
data <- data[ which(data$ColonRectum== 'Colon'), ]
data <- data[ which( data$STAGE_AM == '1' | data$STAGE_AM == '2' | data$STAGE_AM == '3'), ]
data <- data[ which(data$T_neoadjuvant != 'Yes'), ]
dat <- data
colnames(dat)     

x  <-  dat[,c(1,256)]
x  <-  na.omit(x)
head(x)

quant_cut <- function(x, n) {
  qs <- quantile(x, 1:(n-1)/n)
  brks <- c(-Inf, qs, Inf)
  cut(x, breaks=brks, labels=FALSE)
}

x$SIA_tri <-quant_cut(x$SIA, 3)
x<- x[-2]
dat <- merge(dat, x,by="UCAN_TMA_Patnr")
data <- dat

data$OS_surgery_to_last_followup <- as.numeric(data$OS_surgery_to_last_followup)
data$status <-(data$DEAD_ALIVE.1)
data$status   <- ifelse((data$status =="Dead"), 1, 0)
data$OS_surgery_to_last_followup <-  ifelse(( data$OS_surgery_to_last_followup >0), data$OS_surgery_to_last_followup, NA) 
data <- data[!is.na(data$OS_surgery_to_last_followup), ] 
data$TIME <-data$OS_surgery_to_last_followup
data$SurvObj <- with(data, Surv(data$OS_surgery_to_last_followup, status == "1"))


# select stage II only
data <- data[ which( data$STAGE_AM == '2' ), ]
# Kaplan-Meier for SIA OS

datF <- data
datF <-datF[ which(!is.na(datF$SIA_tri)), ]
datF$SIA_tri <- as.factor(datF$SIA_tri)
datF$SIA_tri <- relevel(datF$SIA_tri, ref = "1")

fit <- survfit(SurvObj ~ SIA_tri, data = datF, conf.type = "log-log")

c<-coxph(SurvObj ~ SIA_tri, data = datF)
c
ggsurvplot(fit, risk.table = T, lwd=1,pval = TRUE)

#####











###########################################################################
################ ***         SIA in CRC st IV         *** #################
################ ***                                  *** #################
###########################################################################


setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub")
dat <- NULL
dat <- read.table("Database_working_CRC.txt", sep="\t", header=T, fill = T)
dat$UCAN_TMA_Patnr <- as.character(dat$UCAN_TMA_Patnr)
dat$SIA <- dat$CD8/(dat$CD8+dat$M2)

data <- dat

data <- data[ which( data$STAGE_AM == '4' ), ]
data <- data[ which(data$T_neoadjuvant != 'Yes'), ]
dat <- data
colnames(dat)     

x  <-  dat[,c(1,256)]
x  <-  na.omit(x)
head(x)

quant_cut <- function(x, n) {
  qs <- quantile(x, 1:(n-1)/n)
  brks <- c(-Inf, qs, Inf)
  cut(x, breaks=brks, labels=FALSE)
}

x$SIA_tri <-quant_cut(x$SIA, 3)
x<- x[-2]
dat <- merge(dat, x,by="UCAN_TMA_Patnr")
data <- dat

data$OS_surgery_to_last_followup <- as.numeric(data$OS_surgery_to_last_followup)
data$status <-(data$DEAD_ALIVE.1)
data$status   <- ifelse((data$status =="Dead"), 1, 0)
data$OS_surgery_to_last_followup <-  ifelse(( data$OS_surgery_to_last_followup >0), data$OS_surgery_to_last_followup, NA) 
data <- data[!is.na(data$OS_surgery_to_last_followup), ] 
data$TIME <-data$OS_surgery_to_last_followup
data$SurvObj <- with(data, Surv(data$OS_surgery_to_last_followup, status == "1"))


# Kaplan-Meier for SIA OS

datF <- data
datF <-datF[ which(!is.na(datF$SIA_tri)), ]
datF$SIA_tri <- as.factor(datF$SIA_tri)
datF$SIA_tri <- relevel(datF$SIA_tri, ref = "1")

fit <- survfit(SurvObj ~ SIA_tri, data = datF, conf.type = "log-log")

c<-coxph(SurvObj ~ SIA_tri, data = datF)
c
ggsurvplot(fit, risk.table = T, lwd=1,pval = TRUE)

#####



