########################
########################
require(gdata)
library(reshape)
library(plyr)
library(scales)
library(survival)
library(boot)
library(risksetROC)
library(ggplot2)
library(dplyr)
library(tidyverse)

setwd("/Volumes/Samsung_X5/Cohorts/CRC/Paper current SIA/V 2.0/codes for GitHub")

datS <- NULL
datS <- data.table::fread("Database_working_Validation_cohort.txt")


datS<-as.data.frame(datS)







########
########            OS

datS <- datS[ which(datS$Time_Diagnosis_Last_followup..OS. != "missing"), ]
datS$Time_Diagnosis_Last_followup..OS. <- as.numeric(datS$Time_Diagnosis_Last_followup..OS.)
datS$status <-datS$Event_last_followup
datS$status   <- ifelse((datS$status =="Dead"), 1, 0)
datS$SurvObj <- with(datS, Surv(datS$Time_Diagnosis_Last_followup..OS., status == "1"))
datS <- datS[ which(!is.na(datS$Time_Diagnosis_Last_followup..OS.)), ]
datS <- datS[ which(datS$Time_Diagnosis_Last_followup..OS. > 0), ]

datS$TIME <-datS$Time_Diagnosis_Last_followup..OS.

########       compute SIA

datS$SIA <- datS$CD8/(datS$CD8+datS$M2)

head(datS)
colnames(datS)






########   Survival for melanoma (dichotomised)
datF <- datS  

  d<-NULL
  d <- datF[ which(datF$TMA_cohort == "MEL"), ]

  d<-d[,c(39:37)]
  d <- na.omit(d)
  
  hist(d$SIA, breaks=20, col="red")
  sum(d$SIA==1)
  

  d$SIA_tri <- ifelse((d$SIA>=median(d$SIA)), 2, 1)
  km <- survfit(SurvObj ~ SIA_tri, data = d, conf.type = "log-log")
  plot(km, mark.time = T, lwd=2, col=1:2)
  
  c<-coxph(SurvObj ~ SIA_tri, data = d)
  c<-summary(c)$coefficients[5]
  summary(c)
  c
  c<-as.character(round(c, digits = 5))
  c<-paste("p=", c, sep = "")
  legend(x="bottomleft",legend=c, bty = "n",cex=1.5)
  legend(x="topright", col=1:2, lwd=2, legend=(1:2))



#######
#######    Survival for SIA-sensitive cancers stratified by SIA-groups
######

datF <- datS
  datF <-datF[ which((datF$TMA_cohort=="UBladder" |  datF$TMA_cohort=="Esoph" |  datF$TMA_cohort=="Lung adenocarcinoma")), ]
  
par(mfrow=c(3,4))
for(i in unique(datF$TMA_cohort)){                            
  d<-NULL
  d <- datF[ which(datF$TMA_cohort == i), ]
  d<-d[,c(39:37)]
  d <- na.omit(d)
  
  quant_cut <- function(x, n) {
    qs <- quantile(x, 1:(n-1)/n)
    brks <- c(-Inf, qs, Inf)
    cut(x, breaks=brks, labels=FALSE)
  }
  
  d$SIA_tri <-quant_cut(d$SIA, 3)
  
  d$SIA_tri <- as.factor(d$SIA_tri)
  d$SIA_tri <- relevel(d$SIA_tri, ref = "1")
  
  km <- survfit(SurvObj ~ SIA_tri, data = d, conf.type = "log-log")
  plot(km, mark.time = T, lwd=2, col=1:5,
       main=i)
  
  c<-coxph(SurvObj ~ SIA_tri, data = d)
  p_Med_to_Low<-summary(c)$coefficients[9]
  p_Med_to_Low<-as.character(round(p_Med_to_Low, digits = 3))
  p_High_to_Low<-summary(c)$coefficients[10]
  p_High_to_Low<-as.character(round(p_High_to_Low, digits = 3))
  
  
  c<-paste("p=", p_Med_to_Low , "\np=", p_High_to_Low ,  sep = "")
  legend(x="bottomleft",legend=c, bty = "n",cex=1.5)   
  
}

legend(x="topright", col=1:5, lwd=2, legend=(1:3))




#######
#######    Survival for SIA-unsensitive cancers not stratified by SIA-groups
######





datF <- datS
datF <-datF[ which((datF$TMA_cohort=="Endometrium" |  datF$TMA_cohort=="Ovca_Lund" |  datF$TMA_cohort=="Lung squamous cell carcinoma")), ]


par(mfrow=c(3,4))
for(i in unique(datF$Tumor_type)){                            
  d<-NULL
  d <- datF[ which(datF$Tumor_type == i), ]
  d<-d[,c(39:37)]
  d <- na.omit(d)
  par(mar=c(2,2,2,2))
  
  quant_cut <- function(x, n) {
    qs <- quantile(x, 1:(n-1)/n)
    brks <- c(-Inf, qs, Inf)
    cut(x, breaks=brks, labels=FALSE)
  }
  
  d$SIA_tri <-quant_cut(d$SIA, 3)
  
  
  km <- survfit(SurvObj ~ SIA_tri, data = d, conf.type = "log-log")
  plot(km, mark.time = T, lwd=2, col=1:3,
       main=i)
  


  c<-coxph(SurvObj ~ SIA_tri, data = d)
  c<-summary(c)$coefficients[5]
  c<-as.character(round(c, digits = 5))
  c<-paste("p=", c, sep = "")
  legend(x="bottomleft",legend=c, bty = "n",cex=1.5)   
  
}
legend(x="topright", col=1:3, lwd=2, legend=(1:3))








###########################################
###########################################           compute iAUCs for SIA
###########################################

###            ###
#### for MEL  ####
###            ###
          datF <- datS
          datF2 <- datF[ which(datF$TMA_cohort == "MEL"),  ]
          
                          head(datF2)
                          dat <-  datF2
                          colnames(dat)     
                          head(dat)    
                          x  <-  dat[,c(2,39)]
                          x  <-  na.omit(x)
                          dat <- right_join(dat, x,by="Case")
                          dat$SIA_tri <- ifelse((dat$SIA.x>=median(dat$SIA.x)), 2, 1)
                          dat$SIA_tri <- as.factor(dat$SIA_tri)
                          dat$SIA_tri <- relevel(dat$SIA_tri, ref = "1")
                          datF2 <- dat
          
          i <- NULL
          test <- datF2$TIME
          score <- datF2$SIA_tri
          
          dx <- datF2$status
          df <- data.frame(cbind(test, dx,score)) 
          str(df)
          r <- 1000
          boot.f2 <- function(d, i){
            data <- d[i,]
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
          SIA <-roc.boot$t
          
          SIA_MEL <- SIA 
          
  ###         ###
  #### ** ** ####
  ###         ###         


                ###              ###
                #### for Esohp  ####
                ###              ###
                datF <- datS
                
                datF2 <- datF[ which(datF$TMA_cohort == "Esoph"),  ]

                head(datF2)
                dat <-  datF2
                colnames(dat)     
                head(dat)    
                x  <-  dat[,c(2,39)]
                x  <-  na.omit(x)
                head(x)
                quant_cut <- function(x, n) {
                  qs <- quantile(x, 1:(n-1)/n)
                  brks <- c(-Inf, qs, Inf)
                  cut(x, breaks=brks, labels=FALSE)
                }
                x$SIA_tri <-quant_cut(x$SIA, 3)
                dat <- right_join(dat, x,by="Case")
                dat$SIA_tri <- as.factor(dat$SIA_tri)
                dat$SIA_tri <- relevel(dat$SIA_tri, ref = "1")
                datF2 <- dat
                
                
                i <- NULL
                test <- datF2$TIME
                score <- datF2$SIA_tri
                
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
                SIA <-roc.boot$t
      
                SIA_Esoph <- SIA  

                ###         ###
                #### ** ** ####
                ###         ###         
                
                
                
                    ###                   ###
                    #### for Lung adeno  ####
                    ###                   ###
                    datF <- datS

                    datF2 <- datF[ which(datF$TMA_cohort == "Lung adenocarcinoma"),  ]
                    

                    head(datF2)
                    dat <-  datF2
                    colnames(dat)     
                    head(dat)    
                    x  <-  dat[,c(2,39)]
                    x  <-  na.omit(x)
                    head(x)
                    quant_cut <- function(x, n) {
                      qs <- quantile(x, 1:(n-1)/n)
                      brks <- c(-Inf, qs, Inf)
                      cut(x, breaks=brks, labels=FALSE)
                    }
                    x$SIA_tri <-quant_cut(x$SIA, 3)
                    dat <- right_join(dat, x,by="Case")
                    dat$SIA_tri <- as.factor(dat$SIA_tri)
                    dat$SIA_tri <- relevel(dat$SIA_tri, ref = "1")
                    datF2 <- dat
                    
                    
                    i <- NULL
                    test <- datF2$TIME
                    score <- datF2$SIA_tri
                    
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
                    SIA <-roc.boot$t
                    
                    SIA_Lung<- SIA

                    
                    ###         ###
                    #### ** ** ####
                    ###         ###         
                    
          
                              ###                ###
                              #### for UBladder ####
                              ###                ###
                              datF <- datS

                              datF2 <- datF[ which(datF$TMA_cohort == "UBladder"),  ]

                              
                              head(datF2)
                              dat <-  datF2
                              colnames(dat)     
                              head(dat)    
                              x  <-  dat[,c(2,39)]
                              x  <-  na.omit(x)
                              head(x)
                              quant_cut <- function(x, n) {
                                qs <- quantile(x, 1:(n-1)/n)
                                brks <- c(-Inf, qs, Inf)
                                cut(x, breaks=brks, labels=FALSE)
                              }
                              x$SIA_tri <-quant_cut(x$SIA, 3)
                              dat <- right_join(dat, x,by="Case")
                              dat$SIA_tri <- as.factor(dat$SIA_tri)
                              dat$SIA_tri <- relevel(dat$SIA_tri, ref = "1")
                              datF2 <- dat
                              
                              
                              i <- NULL
                              test <- datF2$TIME
                              score <- datF2$SIA_tri
                              
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
                              SIA <-roc.boot$t

                              SIA_Bladder <- SIA 
                              
                              ###         ###
                              #### ** ** ####
                              ###         ###         
                    

###########################################
###########################################          AUCs   for ISlike
###########################################


### compute Immunoscore-like

dat2 <- datS

for(i in names(dat2[,c(33,34)] )){
  
  x  <-  as.data.frame(dat2[!is.na(dat2[i]),])
  x[i] <- (rank(x[[i]]) - 1) / (length(x[[i]])-1)
  
  f <- "Case"
  dat2 <- merge( dat2,  x[,c(f,i)], by = f, all = T)
}

dat2$ISLike <- rowMeans(dat2[,c(40,41)], na.rm=TRUE)
dat2$ISLike_temp <-  dat2$ISLike
dat2$ISLike_temp <-  ifelse((dat2$ISLike>0.25), 2, 1)
dat2$ISLike <-  ifelse((dat2$ISLike>0.70), 3, dat2$ISLike_temp)

#####   done Immunoscore   ####


        ###                ###
        #### for Esoph    ####
        ###                ###
        
        datF <- dat2
        
        datF2 <- datF[ which(datF$TMA_cohort == "Esoph"),  ]

        dat <-  datF2
        x  <-  dat[,c(1,42)]
        x  <-  na.omit(x)
        dat <- right_join(dat, x,by="Case")
        
        dat$ISLike_tri <- as.factor(dat$ISLike.y)
        dat$ISLike_tri <- relevel(dat$ISLike_tri, ref = "1")
        datF2 <- dat
        
        i <- NULL
        
        
        test <- datF2$TIME
        score <- datF2$ISLike_tri
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
        ISLike <-roc.boot$t
        
        ISLike_Esoph <- ISLike 

        ###         ###
        #### ** ** ####
        ###         ###  


              ###                   ###
              #### for Lung adeno ####
              ###                    ###
              
              datF <- dat2

              datF2 <- datF[ which(datF$TMA_cohort == "Lung adenocarcinoma"),  ]
              
              dat <-  datF2
              x  <-  dat[,c(1,42)]
              x  <-  na.omit(x)
              dat <- right_join(dat, x,by="Case")
              
              dat$ISLike_tri <- as.factor(dat$ISLike.y)
              dat$ISLike_tri <- relevel(dat$ISLike_tri, ref = "1")
              datF2 <- dat
              
              i <- NULL
              
              
              test <- datF2$TIME
              score <- datF2$ISLike_tri
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
              ISLike <-roc.boot$t
              
              ISLike_Lung<- ISLike

              
              ###         ###
              #### ** ** ####
              ###         ###  
              
                  ###                ###
                  #### for MEL      ####
                  ###                ###
                  
                  datF <- dat2
                  
                  datF2 <- datF[ which(datF$TMA_cohort == "MEL"),  ]
                  
                  dat <-  datF2
                  x  <-  dat[,c(1,42)]
                  x  <-  na.omit(x)
                  dat <- right_join(dat, x,by="Case")
                  
                  dat$ISLike_tri <- as.factor(dat$ISLike.y)
                  dat$ISLike_tri <- relevel(dat$ISLike_tri, ref = "1")
                  datF2 <- dat
                  
                  i <- NULL
                  
                  
                  test <- datF2$TIME
                  score <- datF2$ISLike_tri
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
                  ISLike <-roc.boot$t

                 ISLike_MEL <- ISLike 
                  

                  ###         ###
                  #### ** ** ####
                  ###         ###  
                  
                  
                        ###                ###
                        #### for UBladder ####
                        ###                ###
                        
                        datF <- dat2

                        datF2 <- datF[ which(datF$TMA_cohort == "UBladder"),  ]
                        
                        dat <-  datF2
                        x  <-  dat[,c(1,42)]
                        x  <-  na.omit(x)
                        dat <- right_join(dat, x,by="Case")
                        
                        dat$ISLike_tri <- as.factor(dat$ISLike.y)
                        dat$ISLike_tri <- relevel(dat$ISLike_tri, ref = "1")
                        datF2 <- dat
                        
                        i <- NULL
                        
                        
                        test <- datF2$TIME
                        score <- datF2$ISLike_tri
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
                        ISLike <-roc.boot$t
    
                        ISLike_Bladder <- ISLike
                        
                        
                        
                        ###         ###
                        #### ** ** ####
                        ###         ###  


                        ########################
                        ########################

                
###
###
###         SIA and ISLike
###

median(SIA_Bladder)
median(SIA_MEL)
median(ISLike_Bladder)
median(SIA_MEL)
median(SIA_Lung)

Factprs_AUC <- cbind(ISLike_Bladder, SIA_Bladder, ISLike_Esoph, SIA_Esoph,
                     ISLike_Lung,SIA_Lung, ISLike_MEL,SIA_MEL)
Factprs_AUC <- as.data.frame(Factprs_AUC)
names(Factprs_AUC)[c(1:8)]<-c("Bladder_ISLike", "Bladder_SIA",
                              "Gastro-esoph_ISLike","Gastro-esoph_SIA", 
                              "Lung_ISLike","Lung_SIA","Melanoma_ISLike", "Melanoma_SIA")
Factprs_AUC<-tibble::rownames_to_column(Factprs_AUC)
head(Factprs_AUC)                                      

summary<- summary(Factprs_AUC) 
write.table(summary, file ="four AUC summary.txt", append = FALSE, sep="\t", row.names = F, col.names = TRUE)

df_plot <- reshape2::melt(Factprs_AUC, id.vars="rowname")
df_plot$variable <- fct_rev(df_plot$variable)

k <- ggplot(df_plot, aes(x=reorder(variable, desc(variable)), y=value, fill=variable))+
  scale_fill_brewer(palette="RdBu")


k+ 
  theme_bw()+ theme(panel.border = element_blank())+ theme(axis.line = element_line(colour = "black"))+
  # geom_point(aes(colour=variable,group = variable,y=value), size=2.3, alpha=1, position=position_dodge(width=0))+
  geom_boxplot(outlier.shape = NA, stat = "boxplot", width=0.7, alpha=1, position=position_dodge(0.9))+
  theme(legend.position = "none",
        text = element_text(size=12),
        axis.text.x=element_text(size=rel(1), angle=45))+ scale_y_continuous(limits = c(0.5,0.7), breaks = seq(0.45,0.78, by = 0.05))


# for statistics:

SIA_Lung <- as.data.frame(SIA_Esoph)
SIA_Lung$group <- "SIA"          

ISLike_Lung <- as.data.frame(ISLike_Esoph)
ISLike_Lung$group <- "IS" 

IS_and_SIA <- rbind(ISLike_Lung, SIA_Lung)
wilcox.test(V1 ~ group,IS_and_SIA,
            exact = FALSE, paired=F)




                        