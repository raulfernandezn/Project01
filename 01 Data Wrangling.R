#Project 01
#First define a directory with CET file
#Libraries
library(VIM)
library(psych)
library(GPArotation)
library(corrplot)
library(FactoMineR)
library(factoextra)
#Upload data
DF <- read.csv("CET.csv",sep = ",")
rownames(DF) <- DF$Psychology.number
#Isolate groups in CET with index
#CET A
i1 <- which(grepl("CET.A.Q",names(DF)))
#CET B
i2 <- which(grepl("CET.B.Q",names(DF)))
#########################Complete case scenario##################################
#Filter columns CET
s1 <- DF[,i1]
s2 <- DF[,i2]
cs1 <- s1[complete.cases(s1),]
cs2 <- s2[complete.cases(s2),]
########################################MCA Complete Case#####################################
#Factor conversion
rcs1 <- cs1
rcs2 <- cs2
for(i in 1:9)
{
  rcs1[,i] <- factor(rcs1[,i])
  rcs2[,i] <- factor(rcs2[,i])
}
#Factor analysis
ores.mca <- MCA(rcs1,ncp = 3, graph = FALSE)
#Eigen values
oeig.val <- get_eigenvalue(ores.mca)
#Variance
ovar <- get_mca_var(ores.mca)
#Contributions
ovar$eta2
#Scores for Individuals
oind <- get_mca_ind(ores.mca)
ocoorda <- as.data.frame(oind$coord)
ocoorda$id <- rownames(ocoorda)
#Replicate for CET B
#Factor analysis
ores.mcab <- MCA(rcs2,ncp = 3, graph = FALSE)
#Eigen values
oeig.valb <- get_eigenvalue(ores.mcab)
#Variance
ovarb <- get_mca_var(ores.mcab)
#Contributions
ovarb$eta2
#Score for Individuals
oindb <- get_mca_ind(ores.mcab)
ocoordab <- as.data.frame(oindb$coord)
ocoordab$id <- rownames(ocoordab)
################################################################################################
#####################################PCA Other Scores###########################################
#Homologation for other scores
is <- c(1:5,30:45)
OtherS <- DF[,is]
#Principal component analysis for CET A
oPCA.CETA<-principal(cs1,4,rotate="promax",scores=TRUE)
#Extract score for individuals
oScores.CETA=data.frame(oPCA.CETA$scores)
#Indicators and components
oInd.Com.CETA=cbind(cs1,oScores.CETA)
#Correlation Indicators and Components
oCorCETA=cor(oInd.Com.CETA)
oCorCETA<-oCorCETA[1:9,10:13]
#Variance accounted
oPCA.CETA$Vaccounted
#Principal component analysis for CET B
oPCA.CETB<-principal(cs2,4,rotate="promax",scores=TRUE)
#Extract score for inividuals
oScores.CETB=data.frame(oPCA.CETB$scores)
#Indicators and components
oInd.Com.CETB=cbind(cs2,oScores.CETB)
#Correlation Indicators and Components
oCorCETB=cor(oInd.Com.CETB)
oCorCETB<-oCorCETB[1:9,10:13]
#Imputation and PCA for other scores
v1 <- names(OtherS)[-c(1,2,4,5)]
OtherSImp <- hotdeck(OtherS,variable = v1, domain_var = c("Aetiology"))
vr1 <- which(grepl("_imp",names(OtherSImp)))
OtherSImp <- OtherSImp[,-vr1]
#PCA for scores
oPCA.OtherS<-principal(OtherSImp[,-c(1:5)],4,rotate="promax",scores=TRUE)
#Extract scores for individuals
oScores.OtherS=data.frame(oPCA.OtherS$scores)
#Indicators and components
oInd.Com.OtherS=cbind(OtherSImp,oScores.OtherS)
#Correlation Indicators and Components
oCorOtherS=cor(oInd.Com.OtherS[,-c(1:5)])
oCorOtherS<-oCorOtherS[1:16,17:20]
#Format variables
oInd.Com.OtherS <- oInd.Com.OtherS[,-c(2:4)]
names(oInd.Com.OtherS)[c(19:22)] <- paste0(names(oInd.Com.OtherS)[c(19:22)],".OtherS")
oInd.Com.OtherS$Psychology.number <- trimws(as.character(oInd.Com.OtherS$Psychology.number))
#Join
oInd.Com.CETA$id <- trimws(rownames(oInd.Com.CETA))
oInd.Com.CETB$id <- trimws(rownames(oInd.Com.CETB))
rownames(oInd.Com.CETA) <- NULL
rownames(oInd.Com.CETB) <- NULL
FullCETA <- merge(oInd.Com.CETA,oInd.Com.OtherS[,c(1,19:22)],by.x="id",by.y="Psychology.number")
FullCETB <- merge(oInd.Com.CETB,oInd.Com.OtherS[,c(1,19:22)],by.x="id",by.y="Psychology.number")
#Format demographic information
Dem <- OtherSImp[,c(1:5)]
Dem$Range.Age <- cut(Dem$Age,breaks = c(-Inf,29,39,49,59,69,Inf),include.lowest = T,right = T,dig.lab = 10)
Dem$Range.Education <- cut(Dem$Years.of.Edu, breaks = c(-Inf,11,15,22),include.lowest = T,right = T,dig.lab = 10)
Dem$Gender <- factor(Dem$Gender)
Dem$Aetiology <- factor(Dem$Aetiology)
Dem$Psychology.number <- trimws(as.character(Dem$Psychology.number))
Dem <- Dem[,c(1,4:7)]
#Final join
FullCETA <- merge(FullCETA,Dem,by.x="id",by.y="Psychology.number")
FullCETB <- merge(FullCETB,Dem,by.x="id",by.y="Psychology.number")
#Add score
DF$Psychology.number <- trimws(as.character(DF$Psychology.number))
FullCETA <- merge(FullCETA,DF[,c(1,26)],by.x="id",by.y="Psychology.number")
FullCETB <- merge(FullCETB,DF[,c(1,27)],by.x="id",by.y="Psychology.number")
#Remove items
ria <- which(grepl("CET.A.Q",names(FullCETA)))
rib <- which(grepl("CET.B.Q",names(FullCETB)))
FullCETA <- FullCETA[,-ria]
FullCETB <- FullCETB[,-rib]
#Modification for factors
RFullCETA <- FullCETA[,-c(2:5)]
RFullCETB <- FullCETB[,-c(2:5)]
RFullCETA <- merge(RFullCETA,ocoorda,by.x="id",by.y = "id",sort = F)
RFullCETB <- merge(RFullCETB,ocoordab,by.x="id",by.y = "id",sort = F)
########################Bootstrapping version#################################################
#Preprocess
ns1 <- DF[,c(1,5,i1)]
ns2 <- DF[,c(1,5,i2)]
#PCA CET A list objects
set.seed(1)
ldlist <- list(0)
evlist <- list(0)
scoreslist <- list(0)
corlist <- list(0)
#MCA CET A list objects
eiglist <- list(0)
etalist <- list(0)
contriblist <- list(0)
coordlist <- list(0)
#Loop
vars <- names(ns1)[which(grepl("CET",names(ns1)))]
for(i in c(1:1000))
{
  #Hot deck imputation
  DFIA <- hotdeck(ns1,variable = vars, domain_var = c("Aetiology"))
  vri <- which(grepl("_imp",names(DFIA)))
  DFIA <- DFIA[,-vri]
  #Isolate for factor analysis
  FADFIA <- DFIA[,-c(1,2)]
  #Conversion
  for(n in 1:9)
  {
    FADFIA[,n] <- factor(FADFIA[,n])
  }
  #PCA
  PCA.CETA<-principal(scale(DFIA[,-c(1,2)]),4,rotate="promax",scores=TRUE)
  #Loadings
  ld <- as.matrix(unclass(PCA.CETA$loadings))
  #Explained var
  ev <- as.matrix(PCA.CETA$Vaccounted)
  valceta <- ev[3,4]
  #Scores
  Scores.CETA=as.matrix(PCA.CETA$scores)
  dimnames(Scores.CETA)[[1]] <- DFIA$Psychology.number
  #Indicators and components
  Ind.Com.CETA=cbind(DFIA,Scores.CETA)
  #Correlation Indicators and Components
  CorCETA=cor(Ind.Com.CETA[,-c(1,2)])
  CorCETA<-CorCETA[1:9,10:13]
  #Save PCA
  ldlist[[i]] <- ld
  evlist[[i]] <- valceta
  scoreslist[[i]] <- Scores.CETA
  corlist[[i]] <- CorCETA
  #MCA
  #Model
  res.mca <- MCA(FADFIA,ncp = 3, graph = FALSE)
  #Eigen values
  eig.val <- as.matrix(get_eigenvalue(res.mca))
  eig.val <- eig.val[1:3,1:2]
  #Variance
  var <- get_mca_var(res.mca)
  #Eta contribution
  etacont <- as.matrix(var$eta2)
  #Percentage
  mcacontrib <- as.matrix(var$contrib)
  #Individuals
  ind <- get_mca_ind(res.mca)
  facoord <- as.matrix(ind$coord)
  dimnames(facoord)[[1]] <- DFIA$Psychology.number
  #Save MCA
  eiglist[[i]] <- eig.val
  etalist[[i]] <- etacont
  contriblist[[i]] <- mcacontrib
  coordlist[[i]] <- facoord
}
save(ldlist,evlist,scoreslist,corlist,file="Boot.CETA.RData")
save(eiglist,etalist,contriblist,coordlist,file="Boot.CETA.Factor.RData")
#Conglomeration of results
#Loadings PCA
avgloadings <- Reduce('+', ldlist)
avgloadings <- avgloadings/1000
#Explained variance PCA
avgexpvar <- Reduce('+', evlist)
avgexpvar <- avgexpvar/1000
zp <- do.call(c,evlist)
sd(zp)
#Scores PCA
avgscore <- Reduce('+', scoreslist)
avgscore <- avgscore/1000
#Correlation PCA
avgcor <- Reduce('+', corlist)
avgcor <- avgcor/1000
#Conglomeration for MCA
#eigen values
avgeig <- Reduce('+', eiglist)
avgeig <- avgeig/1000
#eta
avgeta <- Reduce('+', etalist)
avgeta <- avgeta/1000
#percentages
avgcontrib <- Reduce('+', contriblist)
avgcontrib <- avgcontrib/1000
#avg coord
avgcoord <- Reduce('+', coordlist)
avgcoord <- avgcoord/1000
#additional format
avgcoord <- as.data.frame(avgcoord)
avgcoord$id <- trimws(rownames(avgcoord))
names(avgcoord)[-4] <- paste0(names(avgcoord)[-4],".MCA")
#Replicate for CET B
rm(i,n)
set.seed(1)
#PCA CET B lists
ldlistb <- list(0)
evlistb <- list(0)
scoreslistb <- list(0)
corlistb <- list(0)
#MCA lists
eiglistb <- list(0)
etalistb <- list(0)
contriblistb <- list(0)
coordlistb <- list(0)
#vars
varsb <- names(ns2)[which(grepl("CET",names(ns2)))]
for(i in c(1:1000))
{
  #Hot deck imputation
  DFIB <- hotdeck(ns2,variable = varsb, domain_var = c("Aetiology"))
  vrib <- which(grepl("_imp",names(DFIB)))
  DFIB <- DFIB[,-vrib]
  #Isolate for factor analysis
  FADFIB <- DFIB[,-c(1,2)]
  #Conversion
  for(n in 1:9)
  {
    FADFIB[,n] <- factor(FADFIB[,n])
  }
  #PCA
  PCA.CETB<-principal(scale(DFIB[,-c(1,2)]),4,rotate="promax",scores=TRUE)
  #Loadings
  ldb <- as.matrix(unclass(PCA.CETB$loadings))
  #Explained var
  evb <- as.matrix(PCA.CETB$Vaccounted)
  valcetb <- evb[3,4]
  #Scores
  Scores.CETB=as.matrix(PCA.CETB$scores)
  dimnames(Scores.CETB)[[1]] <- DFIB$Psychology.number
  #Indicators and components
  Ind.Com.CETB=cbind(DFIB,Scores.CETB)
  #Correlation Indicators and Components
  CorCETB=cor(Ind.Com.CETB[,-c(1,2)])
  CorCETB<-CorCETB[1:9,10:13]
  #Save PCA
  ldlistb[[i]] <- ldb
  evlistb[[i]] <- valcetb
  scoreslistb[[i]] <- Scores.CETB
  corlistb[[i]] <- CorCETB
  #MCA
  res.mcab <- MCA(FADFIB,ncp = 3, graph = FALSE)
  #Eigen values
  eig.valb <- as.matrix(get_eigenvalue(res.mcab))
  eig.valb <- eig.valb[1:3,1:2]
  #Variance
  varb <- get_mca_var(res.mcab)
  #Eta contribution
  etacontb <- as.matrix(varb$eta2)
  #Percentage
  mcacontribb <- as.matrix(varb$contrib)
  #Individuals
  indb <- get_mca_ind(res.mcab)
  facoordb <- as.matrix(indb$coord)
  dimnames(facoordb)[[1]] <- DFIB$Psychology.number
  #Save MCA
  eiglistb[[i]] <- eig.valb
  etalistb[[i]] <- etacontb
  contriblistb[[i]] <- mcacontribb
  coordlistb[[i]] <- facoordb
}
save(ldlistb,evlistb,scoreslistb,corlistb,file="BootCETB.RData")
save(eiglistb,etalistb,contriblistb,coordlistb,file="Boot.CETB.Factor.RData")
#Conglomeration PCA CET B
#Loadings
avgloadingsb <- Reduce('+', ldlistb)
avgloadingsb <- avgloadingsb/1000
#Explained variance
avgexpvarb <- Reduce('+', evlistb)
avgexpvarb <- avgexpvarb/1000
zpb <- do.call(c,evlistb)
sd(zpb)
#Scores
avgscoreb <- Reduce('+', scoreslistb)
avgscoreb <- avgscoreb/1000
#Conglomeration MCA CET B
#eigen values
avgeigb <- Reduce('+', eiglistb)
avgeigb <- avgeigb/1000
#eta
avgetab <- Reduce('+', etalistb)
avgetab <- avgetab/1000
#percentage
avgcontribb <- Reduce('+', contriblistb)
avgcontribb <- avgcontribb/1000
#avg coord
avgcoordb <- Reduce('+', coordlistb)
avgcoordb <- avgcoordb/1000
#Additional format
avgcoordb <- as.data.frame(avgcoordb)
avgcoordb$id <- trimws(rownames(avgcoordb))
names(avgcoordb)[-4] <- paste0(names(avgcoordb)[-4],".MCA")

##################################Dataset merge#############################################
#Conversion
avgscore <- as.data.frame(avgscore)
avgscore$id <- rownames(avgscore)
avgscoreb <- as.data.frame(avgscoreb)
avgscoreb$id <- rownames(avgscoreb)
#ns1 CETA and ns2 CETB
PreCETA <- merge(ns1,avgscore,by.x = "Psychology.number",by.y="id",sort = F)
PreCETB <- merge(ns2,avgscoreb,by.x = "Psychology.number",by.y="id",sort = F)
#Demographic
GCETA <- merge(PreCETA,Dem[,-3],by.x="Psychology.number",by.y="Psychology.number",sort = F)
GCETB <- merge(PreCETB,Dem[,-3],by.x="Psychology.number",by.y="Psychology.number",sort = F)
#Add Total score
GCETA <- merge(GCETA,DF[,c(1,26)],by.x="Psychology.number",by.y="Psychology.number",sort = F)
GCETB <- merge(GCETB,DF[,c(1,27)],by.x="Psychology.number",by.y="Psychology.number",sort = F)

##############################Bootstrapped PCA for Other scores#################################
#Select variables
vo <- names(OtherS)[-c(1,2,3,4,5)]
OtherS <- OtherS[,-c(2,3,4)]
set.seed(1)
#List objects
ldlisto <- list(0)
evlisto <- list(0)
scoreslisto <- list(0)
corlisto <- list(0)
rm(i)
for(i in c(1:1000))
{
  #Hot deck imputation
  DFIO <- hotdeck(OtherS,variable = vo, domain_var = c("Aetiology"))
  vrio <- which(grepl("_imp",names(DFIO)))
  DFIO <- DFIO[,-vrio]
  #PCA
  PCA.CETO<-principal(scale(DFIO[,-c(1,2)]),4,rotate="promax",scores=TRUE)
  #Loadings
  ldo <- as.matrix(unclass(PCA.CETO$loadings))
  #Explained var
  evo <- as.matrix(PCA.CETO$Vaccounted)
  valceto <- evo[3,4]
  #Scores
  Scores.CETO=as.matrix(PCA.CETO$scores)
  dimnames(Scores.CETO)[[1]] <- DFIO$Psychology.number
  #Indicators and components
  Ind.Com.CETO=cbind(DFIO,Scores.CETO)
  #Correlation Indicators and Components
  CorCETO=cor(Ind.Com.CETO[,-c(1,2)])
  CorCETO<-CorCETO[1:16,17:20]
  #Save
  ldlisto[[i]] <- ldo
  evlisto[[i]] <- valceto
  scoreslisto[[i]] <- Scores.CETO
  corlisto[[i]] <- CorCETO
}
save(ldlisto,evlisto,scoreslisto,corlisto,file="BootCETO.RData")
#Conglomeration
#Loadings
avgloadingso <- Reduce('+', ldlisto)
avgloadingso <- avgloadingso/1000
#Explained variance
avgexpvaro <- Reduce('+', evlisto)
avgexpvaro <- avgexpvaro/1000
zpo <- do.call(c,evlisto)
sd(zpo)
#Scores
avgscoreo <- Reduce('+', scoreslisto)
avgscoreo <- avgscoreo/1000
#Format conversion
avgscoreo <- as.data.frame(avgscoreo)
avgscoreo$id <- rownames(avgscoreo)
#Set names
names(avgscoreo)[-5] <- paste0(names(avgscoreo)[-5],".OtherS")
#Definitive frames
TotalCETA <- merge(GCETA,avgscoreo,by.x="Psychology.number",by.y="id",sort = F)
TotalCETB <- merge(GCETB,avgscoreo,by.x="Psychology.number",by.y="id",sort = F)
#Remove items
ia <- which(grepl("CET.A.Q",names(TotalCETA)))
ib <- which(grepl("CET.B.Q",names(TotalCETB)))
TotalCETA <- TotalCETA[,-ia]
TotalCETB <- TotalCETB[,-ib]
TotalCETA$Aetiology <- factor(TotalCETA$Aetiology)
TotalCETB$Aetiology <- factor(TotalCETB$Aetiology)
#Format MCA variables
TotalCETA.Factors <- TotalCETA[,-c(3:6)]
TotalCETB.Factors <- TotalCETB[,-c(3:6)]
#Merge
TotalCETA.Factors <- merge(TotalCETA.Factors,avgcoord,by.x="Psychology.number",by.y="id",sort = F)
TotalCETB.Factors <- merge(TotalCETB.Factors,avgcoordb,by.x="Psychology.number",by.y="id",sort = F)
##Isolate elements to save
save(TotalCETA,TotalCETB,TotalCETA.Factors,TotalCETB.Factors,FullCETA,FullCETB,
     RFullCETA,RFullCETB,file="CleanDF.RData")
save.image("Project01.Data.Total.RData")