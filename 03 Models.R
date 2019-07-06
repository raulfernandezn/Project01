#Models
#Define as directory the route where objects from Script 01 Data Wrangling are saved
#Libraries
library(lme4)
library(Amelia)
library(merTools)
library(plyr)
library(lmerTest)
library(RLRsim)
library(reshape2)
#Functions
#Function p val
compute.pval <- function(dfe,dre,opt)
{
  val <- dim(dfe)[1]+3
  dfe$t.stat <- dfe$mean/dfe$sd
  dre$t.stat <- dre$mean/dre$sd
  dfe$pval <- ifelse(dfe$t.stat>=0,2*(1-pt(dfe$t.stat,df = 379-val)),2*pt(dfe$t.stat,df = 379-val))
  dre$pval <- ifelse(dre$t.stat>=0,2*(1-pt(dre$t.stat,df = 379-val)),2*pt(dre$t.stat,df = 379-val))
  if(opt==1)
  {
    return(dfe)
  } else
  {
    return(dre)
  }
}
#Load Data
load("CleanDF.RData")
#################################(Complete Case PCA)################################################
#Preprocess
names(FullCETA)[2:5] <- paste0("PC",1:4)
names(FullCETA)[6:9] <- paste0("PC",1:4,".OtherS")
#Model Education vs Aetiology
mod.cc.CETA0 <- lmerTest::lmer(CET.A.Raw.Score~1+PC1+PC2+PC3+PC4
                               +PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                 PC4.OtherS+Range.Age+Gender
                               +(1+Range.Education|Aetiology),data = FullCETA,REML = T)
#Model Gender vs Aetiology
mod.cc.CETA1 <- lmerTest::lmer(CET.A.Raw.Score~1+PC1+PC2+PC3+PC4
                               +PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                 PC4.OtherS+Range.Age+Range.Education
                               +(1+Gender|Aetiology),data = FullCETA,REML = T)
#Model Age vs Aetiology
mod.cc.CETA2 <- lmerTest::lmer(CET.A.Raw.Score~1+PC1+PC2+PC3+PC4
                               +PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                 PC4.OtherS+Range.Age+Range.Education
                               +(1+Range.Age|Aetiology),data = FullCETA,REML = T)
##################################(Mult. Imp. PCA CET A)####################################
#Full imputed CET A with PCA
#Pre process
names(TotalCETA)[3:6] <- paste0("PC",1:4)
names(TotalCETA)[11:14] <- paste0("PC",1:4,".OtherS")
#Impute
set.seed(1)
imCETA <- amelia(TotalCETA,idvars = c("Psychology.number","Gender","Aetiology","Range.Age","Range.Education"),m = 10)
#models
#Model Education vs Aetiology
mods0 <- lapply(imCETA$imputations,
                function(d) lmer(CET.A.Raw.Score~1+PC1+PC2+PC3+PC4
                                 +PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                   PC4.OtherS+Range.Age+Gender
                                 +(1+Range.Education|Aetiology),data=d,REML = T))
#Model Gender vs Aetiology
mods1 <- lapply(imCETA$imputations,
                function(d) lmer(CET.A.Raw.Score~1+PC1+PC2+PC3+PC4
                                 +PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                   PC4.OtherS+Range.Age+Range.Education
                                 +(1+Gender|Aetiology),data=d,REML = T))
#Model Age vs Aetiology
mods2 <- lapply(imCETA$imputations,
                function(d) lmer(CET.A.Raw.Score~1+PC1+PC2+PC3+PC4
                                 +PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                   PC4.OtherS+Range.Age+Range.Education
                                 +(1+Range.Age|Aetiology),data=d,REML = T))
####################################Global Effects##############################################
#mods0
imputeFEs0 <- ldply(mods0, FEsim,seed=1)
imputeREs0 <- ldply(mods0, REsim,seed=1)
imputeFEs0 <- ddply(imputeFEs0, .(term), summarize, mean = mean(mean), 
                    median = mean(median), sd = mean(sd))
imputeREs0 <- ddply(imputeREs0, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                    median = mean(median), sd = mean(sd))
#Model Bootstrapped CET A Education vs Aetiology
imputeFEs0 <- compute.pval(imputeFEs0,imputeREs0,1)
imputeREs0 <- compute.pval(imputeFEs0,imputeREs0,0)
#mods1
imputeFEs1 <- ldply(mods1, FEsim)
imputeREs1 <- ldply(mods1, REsim)
imputeFEs1 <- ddply(imputeFEs1, .(term), summarize, mean = mean(mean), 
                    median = mean(median), sd = mean(sd))
imputeREs1 <- ddply(imputeREs1, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                    median = mean(median), sd = mean(sd))
#Model Bootstrapped CET A Gender vs Aetiology
imputeFEs1 <- compute.pval(imputeFEs1,imputeREs1,1)
imputeREs1 <- compute.pval(imputeFEs1,imputeREs1,0)
#mods2
imputeFEs2 <- ldply(mods2, FEsim)
imputeREs2 <- ldply(mods2, REsim)
imputeFEs2 <- ddply(imputeFEs2, .(term), summarize, mean = mean(mean), 
                    median = mean(median), sd = mean(sd))
imputeREs2 <- ddply(imputeREs2, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                    median = mean(median), sd = mean(sd))
#Model Bootstrapped CET A Age vs Aetiology
imputeFEs2 <- compute.pval(imputeFEs2,imputeREs2,1)
imputeREs2 <- compute.pval(imputeFEs2,imputeREs2,0)
##########Comparison Model Education vs Aetiology
#CCPCA CET A
#Fixed effect
feccpca.CETA <- as.data.frame(summary(mod.cc.CETA0)$coefficients)
#Random effect
reccpca.CETA <- as.data.frame(ranef(mod.cc.CETA0))
reccpca.CETA$grp <- as.numeric(reccpca.CETA$grp)
s1 <- dcast(grpvar+grp~term,data = reccpca.CETA[,-5])
#BPCA CET A
imputeFEs0
imputeREs0
##########Comparison Model Gender vs Aetiology
#Fixed effect
feccpca.CETA.ga <- as.data.frame(summary(mod.cc.CETA1)$coefficients)
#Random effect
reccpca.CETA.ga <- as.data.frame(ranef(mod.cc.CETA1))
reccpca.CETA.ga$grp <- as.numeric(reccpca.CETA.ga$grp)
s1ga <- dcast(grpvar+grp~term,data = reccpca.CETA.ga[,-5])
#BPCA CET A
imputeFEs1
imputeREs1
##########Comparison Models Age vs Aetiology
#Fixed effect
feccpca.CETA.aa <- as.data.frame(summary(mod.cc.CETA2)$coefficients)
#Random effect
reccpca.CETA.aa <- as.data.frame(ranef(mod.cc.CETA2))
reccpca.CETA.aa$grp <- as.numeric(reccpca.CETA.aa$grp)
s1aa <- dcast(grpvar+grp~term,data = reccpca.CETA.aa[,-5])
#BPCA CET A
imputeFEs2
imputeREs2

###############################(Complete case MCA CET A)#########################################
#Full CET A Factors
names(RFullCETA) <- gsub(" ",".",names(RFullCETA))
names(RFullCETA)[2:5] <- paste0("PC",1:4,".OtherS")
names(RFullCETA)[11:13] <- paste0("F",1:3)
#Model
mod.cc.CETAF0 <- lmerTest::lmer(CET.A.Raw.Score~1+PC1.OtherS+PC2.OtherS+
                                  PC3.OtherS+PC4.OtherS+Gender+
                                  Range.Age+F1+F2+F3
                                +(1+Range.Education|Aetiology),data = RFullCETA,REML = T)
mod.cc.CETAF1 <- lmerTest::lmer(CET.A.Raw.Score~1+PC1.OtherS+PC2.OtherS+
                                  PC3.OtherS+PC4.OtherS+Range.Age+
                                  Range.Education+F1+F2+F3
                                +(1+Gender|Aetiology),data = RFullCETA,REML = T)
mod.cc.CETAF2 <- lmerTest::lmer(CET.A.Raw.Score~1+PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                  PC4.OtherS+Gender+Range.Education+F1+F2+F3
                                +(1+Range.Age|Aetiology),data = RFullCETA,REML = T)
###############################(Mult. Imp. MCA CET A)#########################################
#Full imputed CET A with MCA
#format
names(TotalCETA.Factors)<-gsub(" ",".",names(TotalCETA.Factors))
names(TotalCETA.Factors)[7:10] <- paste0("PC",1:4,".OtherS")
names(TotalCETA.Factors)[11:13] <- paste0("F",1:3)
#impute
set.seed(1)
imCETAF <- amelia(TotalCETA.Factors,idvars = c("Psychology.number","Gender","Aetiology","Range.Age","Range.Education"),m = 10)
#models
modsf0 <- lapply(imCETAF$imputations,
                 function(d) lmer(CET.A.Raw.Score~1+Gender+Range.Age+
                                    PC1.OtherS+PC2.OtherS+PC3.OtherS+PC4.OtherS+
                                    F1+F2+F3
                                  +(1+Range.Education|Aetiology),data=d,REML = T))
modsf1 <- lapply(imCETAF$imputations,
                 function(d) lmer(CET.A.Raw.Score~1+Range.Age+Range.Education+PC1.OtherS+
                                    PC2.OtherS+PC3.OtherS+PC4.OtherS+
                                    F1+F2+F3
                                  +(1+Gender|Aetiology),data=d,REML = T))
modsf2 <- lapply(imCETAF$imputations,
                 function(d) lmer(CET.A.Raw.Score~1+Gender+Range.Education+
                                    PC1.OtherS+PC2.OtherS+PC3.OtherS+PC4.OtherS+
                                    F1+F2+F3
                                  +(1+Range.Age|Aetiology),data=d,REML = T))
#################################Global Effects##################################################
#mods0
imputeFEsf0 <- ldply(modsf0, FEsim,seed=1)
imputeREsf0 <- ldply(modsf0, REsim,seed=1)
imputeFEsf0 <- ddply(imputeFEsf0, .(term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeREsf0 <- ddply(imputeREsf0, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeFEsf0 <- compute.pval(imputeFEsf0,imputeREsf0,1)
imputeREsf0 <- compute.pval(imputeFEsf0,imputeREsf0,0)

#mods1
imputeFEsf1 <- ldply(modsf1, FEsim,seed=1)
imputeREsf1 <- ldply(modsf1, REsim,seed=1)
imputeFEsf1 <- ddply(imputeFEsf1, .(term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeREsf1 <- ddply(imputeREsf1, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeFEsf1 <- compute.pval(imputeFEsf1,imputeREsf1,1)
imputeREsf1 <- compute.pval(imputeFEsf1,imputeREsf1,0)
#mods2
imputeFEsf2 <- ldply(modsf2, FEsim,seed=1)
imputeREsf2 <- ldply(modsf2, REsim,seed=1)
imputeFEsf2 <- ddply(imputeFEsf2, .(term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeREsf2 <- ddply(imputeREsf2, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeFEsf2 <- compute.pval(imputeFEsf2,imputeREsf2,1)
imputeREsf2 <- compute.pval(imputeFEsf2,imputeREsf2,0)
##########Comparison Model Education vs Aetiology
#CCFA CET A
#Fixed effect
feccpca.CETAF0 <- as.data.frame(summary(mod.cc.CETAF0)$coefficients)
#Random effect
reccpca.CETAF0 <- as.data.frame(ranef(mod.cc.CETAF0))
reccpca.CETAF0$grp <- as.numeric(reccpca.CETAF0$grp)
t1 <- dcast(grpvar+grp~term,data = reccpca.CETAF0[,-5])
#BPCA CET A
imputeFEsf0
imputeREsf0
##########Comparison Model Gender vs Aetiology
#Fixed effect
feccpca.CETAF1.ga <- as.data.frame(summary(mod.cc.CETAF1)$coefficients)
#Random effect
reccpca.CETAF1.ga <- as.data.frame(ranef(mod.cc.CETAF1))
reccpca.CETAF1.ga$grp <- as.numeric(reccpca.CETAF1.ga$grp)
t1ga <- dcast(grpvar+grp~term,data = reccpca.CETAF1.ga[,-5])
#BPCA CET A
imputeFEsf1
imputeREsf1
##########Comparison Models Age vs Aetiology
#Fixed effect
feccpca.CETAF2.aa <- as.data.frame(summary(mod.cc.CETAF2)$coefficients)
#Random effect
reccpca.CETAF2.aa <- as.data.frame(ranef(mod.cc.CETAF2))
reccpca.CETAF2.aa$grp <- as.numeric(reccpca.CETAF2.aa$grp)
t1aa <- dcast(grpvar+grp~term,data = reccpca.CETAF2.aa[,-5])
#BPCA CET A
imputeFEsf2
imputeREsf2
#################################(Complete case PCA CET B)#######################################
#CC CET B PCA
#FULLCETB
#Preprocess
names(FullCETB)[2:5] <- paste0("PC",1:4)
names(FullCETB)[6:9] <- paste0("PC",1:4,".OtherS")

#Model
mod.cc.CETB0 <- lmerTest::lmer(CET.B.Raw.Score~1+PC1+PC2+PC3+PC4+
                                 PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                 PC4.OtherS+Gender+Range.Age
                               +(1+Range.Education|Aetiology),data = FullCETB,REML = T)
mod.cc.CETB1 <- lmerTest::lmer(CET.B.Raw.Score~1+PC1+PC2+PC3+PC4+
                                 PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                 PC4.OtherS+
                                 Range.Age+Range.Education
                               +(1+Gender|Aetiology),data = FullCETB,REML = T)
mod.cc.CETB2 <- lmerTest::lmer(CET.B.Raw.Score~1+PC1+PC2+PC3+PC4+
                                 PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                 PC4.OtherS+
                                 Gender+Range.Education
                               +(1+Range.Age|Aetiology),data = FullCETB,REML = T)
################################(Mult. Imp. PCA CET B)#############################################
#Full imputed CET B with PCA
#Preprocess
names(TotalCETB)[3:6] <- paste0("PC",1:4)
names(TotalCETB)[11:14] <- paste0("PC",1:4,".OtherS")
#Impute
set.seed(1)
imCETB <- amelia(TotalCETB,idvars = c("Psychology.number","Gender","Aetiology","Range.Age","Range.Education"),m = 10)
#models
modsb0 <- lapply(imCETB$imputations,
                 function(d) lmer(CET.B.Raw.Score~1+PC1+PC2+PC3+PC4+
                                    PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                    PC4.OtherS+Gender+Range.Age
                                  +(1+Range.Education|Aetiology),data=d,REML = T))
modsb1 <- lapply(imCETB$imputations,
                 function(d) lmer(CET.B.Raw.Score~1+PC1+PC2+PC3+PC4+
                                    PC1.OtherS+PC2.OtherS+PC3.OtherS+PC4.OtherS+
                                    Range.Age+Range.Education
                                  +(1+Gender|Aetiology),data=d,REML = T))
modsb2 <- lapply(imCETB$imputations,
                 function(d) lmer(CET.B.Raw.Score~1+PC1+PC2+PC3+PC4+
                                    PC1.OtherS+PC2.OtherS+PC3.OtherS+PC4.OtherS+
                                    Gender+Range.Education
                                  +(1+Range.Age|Aetiology),data=d,REML = T))
#####################################Global Effects############################################
#mods0
imputeFEsb0 <- ldply(modsb0, FEsim)
imputeREsb0 <- ldply(modsb0, REsim)
imputeFEsb0 <- ddply(imputeFEsb0, .(term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeREsb0 <- ddply(imputeREsb0, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeFEsb0 <- compute.pval(imputeFEsb0,imputeREsb0,1)
imputeREsb0 <- compute.pval(imputeFEsb0,imputeREsb0,0)
#mods1
imputeFEsb1 <- ldply(modsb1, FEsim)
imputeREsb1 <- ldply(modsb1, REsim)
imputeFEsb1 <- ddply(imputeFEsb1, .(term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeREsb1 <- ddply(imputeREsb1, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeFEsb1 <- compute.pval(imputeFEsb1,imputeREsb1,1)
imputeREsb1 <- compute.pval(imputeFEsb1,imputeREsb1,0)
#mods2
imputeFEsb2 <- ldply(modsb2, FEsim)
imputeREsb2 <- ldply(modsb2, REsim)
imputeFEsb2 <- ddply(imputeFEsb2, .(term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeREsb2 <- ddply(imputeREsb2, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                     median = mean(median), sd = mean(sd))
imputeFEsb2 <- compute.pval(imputeFEsb2,imputeREsb2,1)
imputeREsb2 <- compute.pval(imputeFEsb2,imputeREsb2,0)
##########Comparison Model Education vs Aetiology
#CCFA CET B
#Fixed effect
feccpca.CETB0 <- as.data.frame(summary(mod.cc.CETB0)$coefficients)
#Random effect
reccpca.CETB0 <- as.data.frame(ranef(mod.cc.CETB0))
reccpca.CETB0$grp <- as.numeric(reccpca.CETB0$grp)
u1 <- dcast(grpvar+grp~term,data = reccpca.CETB0[,-5])
#BPCA CET B
imputeFEsb0
imputeREsb0
##########Comparison Model Gender vs Aetiology
#Fixed effect
feccpca.CETB1.ga <- as.data.frame(summary(mod.cc.CETB1)$coefficients)
#Random effect
reccpca.CETB1.ga <- as.data.frame(ranef(mod.cc.CETB1))
reccpca.CETB1.ga$grp <- as.numeric(reccpca.CETB1.ga$grp)
u1ga <- dcast(grpvar+grp~term,data = reccpca.CETB1.ga[,-5])
#BPCA CET B
imputeFEsb1
imputeREsb1
##########Comparison Models Age vs Aetiology
#Fixed effect
feccpca.CETB2.aa <- as.data.frame(summary(mod.cc.CETB2)$coefficients)
#Random effect
reccpca.CETB2.aa <- as.data.frame(ranef(mod.cc.CETB2))
reccpca.CETB2.aa$grp <- as.numeric(reccpca.CETB2.aa$grp)
u1aa <- dcast(grpvar+grp~term,data = reccpca.CETB2.aa[,-5])
#BPCA CET B
imputeFEsb2
imputeREsb2
##################################(Complete Case MCA CET B)############################################
#Full CET B Factors
##CC CET B FA
names(RFullCETB) <- gsub(" ",".",names(RFullCETB))
names(RFullCETB)[2:5] <- paste0("PC",1:4,".OtherS")
names(RFullCETB)[11:13] <- paste0("F",1:3)
#Model
mod.cc.CETBF0 <- lmerTest::lmer(CET.B.Raw.Score~1+PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                  PC4.OtherS+Gender+Range.Age+F1+F2+F3
                                +(1+Range.Education|Aetiology),data = RFullCETB,REML = T)
mod.cc.CETBF1 <- lmerTest::lmer(CET.B.Raw.Score~1+PC1.OtherS+PC2.OtherS+PC3.OtherS+
                                  PC4.OtherS+Range.Age+Range.Education+F1+F2+F3
                                +(1+Gender|Aetiology),data = RFullCETB,REML = T)
mod.cc.CETBF2 <- lmerTest::lmer(CET.B.Raw.Score~1+PC1.OtherS+PC2.OtherS+
                                  PC3.OtherS+PC4.OtherS+Gender+
                                  Range.Education+F1+F2+F3
                                +(1+Range.Age|Aetiology),data = RFullCETB,REML = T)
##################################(Mult. Imp MCA CET B)##########################################
#format
names(TotalCETB.Factors)<-gsub(" ",".",names(TotalCETB.Factors))
names(TotalCETB.Factors)[7:10] <- paste0("PC",1:4,".OtherS")
names(TotalCETB.Factors)[11:13] <- paste0("F",1:3)
#impute
set.seed(1)
imCETBF <- amelia(TotalCETB.Factors,idvars = c("Psychology.number","Gender","Aetiology","Range.Age","Range.Education"),m = 10)
#models
modsfb0 <- lapply(imCETBF$imputations,
                  function(d) lmer(CET.B.Raw.Score~1+Gender+Range.Age+PC1.OtherS+
                                     PC2.OtherS+PC3.OtherS+PC4.OtherS+
                                     F1+F2+F3
                                   +(1+Range.Education|Aetiology),data=d,REML = T))
modsfb1 <- lapply(imCETBF$imputations,
                  function(d) lmer(CET.B.Raw.Score~1+Range.Age+Range.Education+
                                     PC1.OtherS+PC2.OtherS+PC3.OtherS+PC4.OtherS+
                                     F1+F2+F3
                                   +(1+Gender|Aetiology),data=d,REML = T))
modsfb2 <- lapply(imCETBF$imputations,
                  function(d) lmer(CET.B.Raw.Score~1+Gender+Range.Education+PC1.OtherS+
                                     PC2.OtherS+PC3.OtherS+PC4.OtherS+
                                     F1+F2+F3
                                   +(1+Range.Age|Aetiology),data=d,REML = T))
########################################Global Effects###########################################
#mods0
imputeFEsfb0 <- ldply(modsfb0, FEsim)
imputeREsfb0 <- ldply(modsfb0, REsim)
imputeFEsfb0 <- ddply(imputeFEsfb0, .(term), summarize, mean = mean(mean), 
                      median = mean(median), sd = mean(sd))
imputeREsfb0 <- ddply(imputeREsfb0, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                      median = mean(median), sd = mean(sd))
imputeFEsfb0 <- compute.pval(imputeFEsfb0,imputeREsfb0,1)
imputeREsfb0 <- compute.pval(imputeFEsfb0,imputeREsfb0,0)

#mods1
imputeFEsfb1 <- ldply(modsfb1, FEsim)
imputeREsfb1 <- ldply(modsfb1, REsim)
imputeFEsfb1 <- ddply(imputeFEsfb1, .(term), summarize, mean = mean(mean), 
                      median = mean(median), sd = mean(sd))
imputeREsfb1 <- ddply(imputeREsfb1, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                      median = mean(median), sd = mean(sd))
imputeFEsfb1 <- compute.pval(imputeFEsfb1,imputeREsfb1,1)
imputeREsfb1 <- compute.pval(imputeFEsfb1,imputeREsfb1,0)

#mods2
imputeFEsfb2 <- ldply(modsfb2, FEsim)
imputeREsfb2 <- ldply(modsfb2, REsim)
imputeFEsfb2 <- ddply(imputeFEsfb2, .(term), summarize, mean = mean(mean), 
                      median = mean(median), sd = mean(sd))
imputeREsfb2 <- ddply(imputeREsfb2, .(groupFctr, groupID, term), summarize, mean = mean(mean), 
                      median = mean(median), sd = mean(sd))
imputeFEsfb2 <- compute.pval(imputeFEsfb2,imputeREsfb2,1)
imputeREsfb2 <- compute.pval(imputeFEsfb2,imputeREsfb2,0)
##########Comparison Model Education vs Aetiology
#BFA CET B
#Fixed effect
feccpca.CETBF0 <- as.data.frame(summary(mod.cc.CETBF0)$coefficients)
#Random effect
reccpca.CETBF0 <- as.data.frame(ranef(mod.cc.CETBF0))
reccpca.CETBF0$grp <- as.numeric(reccpca.CETBF0$grp)
v1 <- dcast(grpvar+grp~term,data = reccpca.CETBF0[,-5])
#BPCA CET B
imputeFEsfb0
imputeREsfb0
##########Comparison Model Gender vs Aetiology
#Fixed effect
feccpca.CETBF1.ga <- as.data.frame(summary(mod.cc.CETBF1)$coefficients)
#Random effect
reccpca.CETBF1.ga <- as.data.frame(ranef(mod.cc.CETBF1))
reccpca.CETBF1.ga$grp <- as.numeric(reccpca.CETBF1.ga$grp)
v1ga <- dcast(grpvar+grp~term,data = reccpca.CETBF1.ga[,-5])
#BPCA CET B
imputeFEsfb1
imputeREsfb1
##########Comparison Models Age vs Aetiology
#Fixed effect
feccpca.CETBF2.aa <- as.data.frame(summary(mod.cc.CETBF2)$coefficients)
#Random effect
reccpca.CETBF2.aa <- as.data.frame(ranef(mod.cc.CETBF2))
reccpca.CETBF2.aa$grp <- as.numeric(reccpca.CETBF2.aa$grp)
v1aa <- dcast(grpvar+grp~term,data = reccpca.CETBF2.aa[,-5])
#BPCA CET B
imputeFEsfb2
imputeREsfb2
