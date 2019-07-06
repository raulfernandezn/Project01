library(ggcorrplot)
source("addgrids3d.r")
#This script consumes output from Script 01 Data Wrangling
#####Process for complete case PCA CET A
sum1a <- as.data.frame(oPCA.CETA$Vaccounted)
sum2a <- as.data.frame(round(unclass(oPCA.CETA$loadings),3))
#Loadings
IndLoadCETA <- as.data.frame(unclass(oPCA.CETA$loadings))
#3d plot
p0 <- scatterplot3d::scatterplot3d(IndLoadCETA$RC4,IndLoadCETA$RC1,IndLoadCETA$RC3,pch=16,color="steelblue",
                                   grid=FALSE, box=FALSE,
                                   xlab="Principal Component 1",
                                   ylab="Principal Component 2",
                                   zlab="Principal Component 3",
                                   main="PCA Complete Case CET A")
addgrids3d(IndLoadCETA[, 1:3], grid = c("xy", "xz", "yz"))
p0.coord <- p0$xyz.convert(IndLoadCETA$RC4,IndLoadCETA$RC1,IndLoadCETA$RC3)
text(p0.coord$x,p0.coord$y,labels = rownames(IndLoadCETA),pos = 2,offset = 0.5,adj = c(0,0),cex=0.7,srt=45)
#Process for complete case PCA CET B
sum1a <- as.data.frame(oPCA.CETA$Vaccounted)
sum2a <- as.data.frame(round(unclass(oPCA.CETA$loadings),3))
#Loadings
IndLoadCETB <- as.data.frame(unclass(oPCA.CETB$loadings))
#3d plots
q0 <- scatterplot3d::scatterplot3d(IndLoadCETB$RC1,IndLoadCETB$RC2,IndLoadCETB$RC3,color="steelblue",pch=16,
                                   grid=FALSE, box=FALSE,
                                   xlab="Principal Component 1",
                                   ylab="Principal Component 2",
                                   zlab="Principal Component 3",
                                   main="PCA Complete Case CET B")
addgrids3d(IndLoadCETB[, 1:3], grid = c("xy", "xz", "yz"))
q0.coord <- q0$xyz.convert(IndLoadCETB$RC1,IndLoadCETB$RC2,IndLoadCETB$RC3)
text(q0.coord$x,q0.coord$y,labels = rownames(IndLoadCETB),pos = 2,offset = 0.5,adj = c(0,0),cex=0.7,srt=45)
#Process for one trial PCA Other scores
sum1o <- as.data.frame(oPCA.OtherS$Vaccounted)
sum2o <- as.data.frame(round(unclass(oPCA.OtherS$loadings),3))
#Loadings
IndLoadCETOthers <- as.data.frame(unclass(oPCA.OtherS$loadings))
#3d plots
r0 <- scatterplot3d::scatterplot3d(IndLoadCETOthers$RC1,IndLoadCETOthers$RC3,IndLoadCETOthers$RC2,color="steelblue",pch=16,
                                   grid=FALSE, box=FALSE,
                                   xlab="Principal Component 1",
                                   ylab="Principal Component 2",
                                   zlab="Principal Component 3",
                                   main="PCA Other Score Measures")
addgrids3d(IndLoadCETOthers[, 1:3], grid = c("xy", "xz", "yz"))
r0.coord <- r0$xyz.convert(IndLoadCETOthers$RC1,IndLoadCETOthers$RC3,IndLoadCETOthers$RC2)
text(r0.coord$x,r0.coord$y,labels = rownames(IndLoadCETOthers),pos = 2,offset = 0.5,adj = c(0,0),cex=0.8,srt=45)
#Process for complete case MCA CET
oeig.val <- as.data.frame(oeig.val)
oeig.val <- oeig.val[1:3,]
oeta <- as.data.frame(ovar$eta2)
#3d plots
png(filename = "FPlot01.png")
t0 <- scatterplot3d::scatterplot3d(oeta$`Dim 1`,oeta$`Dim 2`,oeta$`Dim 3`,color="purple",pch=16,
                                   grid=FALSE, box=FALSE,
                                   xlab="Factor 1",
                                   ylab="Factor 2",
                                   zlab="Factor 3",
                                   main="MCA Complete Case CET A")
addgrids3d(oeta[, 1:3], grid = c("xy", "xz", "yz"))
t0.coord <- t0$xyz.convert(oeta$`Dim 1`,oeta$`Dim 2`,oeta$`Dim 3`)
text(t0.coord$x,t0.coord$y,labels = rownames(oeta),pos = 2,offset = 0.7,adj = c(0,0),cex=0.8,srt=45)
dev.off()
#Process for complete case MCA CET B
oeig.valb <- as.data.frame(oeig.valb)
oeig.valb <- oeig.valb[1:3,]
oetab <- as.data.frame(ovarb$eta2)
#3d plots
u0 <- scatterplot3d::scatterplot3d(oetab$`Dim 1`,oetab$`Dim 2`,oetab$`Dim 3`,color="purple",pch=16,
                                   grid=FALSE, box=FALSE,
                                   xlab="Factor 1",
                                   ylab="Factor 2",
                                   zlab="Factor 3",
                                   main="MCA Complete Case CET B")
addgrids3d(oetab[, 1:3], grid = c("xy", "xz", "yz"))
u0.coord <- u0$xyz.convert(oetab$`Dim 1`,oetab$`Dim 2`,oetab$`Dim 3`)
text(u0.coord$x,u0.coord$y,labels = rownames(oetab),pos = 2,offset = 0.7,adj = c(0,0),cex=0.7,srt=45)

#Bootstrapped PCA CET A 
#Process average loadings
avgloadingsformat <- as.data.frame(avgloadings)
##Plot CET A
p1 <- scatterplot3d::scatterplot3d(avgloadingsformat$RC1,avgloadingsformat$RC3,avgloadingsformat$RC4,color="steelblue",pch=16,
                                   grid=FALSE, box=FALSE,
                                   xlab="Principal Component 1",
                                   ylab="Principal Component 2",
                                   zlab="Principal Component 3",
                                   main="PCA Imputed Data CET A")
addgrids3d(avgloadingsformat[, 1:3], grid = c("xy", "xz", "yz"))
p1.coord <- p1$xyz.convert(avgloadingsformat$RC1,avgloadingsformat$RC3,avgloadingsformat$RC4)
text(p1.coord$x,p1.coord$y,labels = rownames(avgloadingsformat),pos = 2,offset = 0.7,adj = c(0,0),cex=0.7,srt=45)
#Average explained variance
avgexpvar
#Montecarlo error CET.A
sd(zp)
#Process for bootstrapped vals CET B
avgloadingsformatb <- as.data.frame(avgloadingsb)
##Plot CET B
q1 <- scatterplot3d::scatterplot3d(avgloadingsformatb$RC2,avgloadingsformatb$RC1,avgloadingsformatb$RC4,color="steelblue",pch=16,
                                   grid=FALSE, box=FALSE,
                                   xlab="Principal Component 1",
                                   ylab="Principal Component 2",
                                   zlab="Principal Component 3",
                                   main="PCA Imputed Data CET B")
addgrids3d(avgloadingsformatb[, 1:3], grid = c("xy", "xz", "yz"))
q1.coord <- q1$xyz.convert(avgloadingsformatb$RC2,avgloadingsformatb$RC1,avgloadingsformatb$RC4)
text(q1.coord$x,q1.coord$y,labels = rownames(avgloadingsformatb),pos = 2,offset = 0.7,adj = c(0,0),cex=0.7,srt=45)
#Cum var
avgexpvarb
#Process for bootstrapped vals Other Scores
avgloadingsformato <- as.data.frame(avgloadingso)
##Plot Others
r1 <- scatterplot3d::scatterplot3d(avgloadingsformato$RC1,avgloadingsformato$RC3,avgloadingsformato$RC2,color="steelblue",pch=16,
                                   grid=FALSE, box=FALSE,
                                   xlab="Principal Component 1",
                                   ylab="Principal Component 2",
                                   zlab="Principal Component 3",
                                   main="PCA Imputed Data Other Score Measures")
addgrids3d(avgloadingsformato[, 1:3], grid = c("xy", "xz", "yz"))
r1.coord <- r1$xyz.convert(avgloadingsformato$RC1,avgloadingsformato$RC3,avgloadingsformato$RC2)
text(r1.coord$x,r1.coord$y,labels = rownames(avgloadingsformato),pos = 2,offset = 0.5,adj = c(0,0),cex=0.7,srt=25)
#Bootstrap MCA CET A
#Eigen vals
avgeig
#Eta
DFFA <- as.data.frame(avgeta)
#Plot
s1 <- scatterplot3d::scatterplot3d(DFFA$`Dim 1`,DFFA$`Dim 2`,DFFA$`Dim 3`,color="purple",pch=16,
                                   grid=FALSE, box=FALSE,
                                   xlab="Factor 1",
                                   ylab="Factor 2",
                                   zlab="Factor 3",
                                   main="MCA Imputed Data CET A")
addgrids3d(DFFA[, 1:3], grid = c("xy", "xz", "yz"))
s1.coord <- s1$xyz.convert(DFFA$`Dim 1`,DFFA$`Dim 2`,DFFA$`Dim 3`)
text(s1.coord$x,s1.coord$y,labels = rownames(DFFA),pos = 2,offset = 0.5,adj = c(0,0),cex=0.7,srt=45)
#Bootstrap MCA CET B
#Eigen vals
avgeigb
#Eta
DFFB <- as.data.frame(avgetab)
#Plot
s2 <- scatterplot3d::scatterplot3d(DFFB$`Dim 1`,DFFB$`Dim 2`,DFFB$`Dim 3`,color="purple",pch=16,
                                   grid=FALSE, box=FALSE,
                                   xlab="Factor 1",
                                   ylab="Factor 2",
                                   zlab="Factor 3",
                                   main="MCA Imputed Data CET B")
addgrids3d(DFFB[, 1:3], grid = c("xy", "xz", "yz"))
s2.coord <- s2$xyz.convert(DFFB$`Dim 1`,DFFB$`Dim 2`,DFFB$`Dim 3`)
text(s2.coord$x,s2.coord$y,labels = rownames(DFFB),pos = 2,offset = 0.5,adj = c(0,0),cex=0.7,srt=45)
dev.off()