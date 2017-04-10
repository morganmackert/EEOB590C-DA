# Multivariate GLM: MANOVA, Multivariate Regression, MANCOVA

         #Packages: geomorph, vegan

rm(list=ls())
library(geomorph) 
library(vegan)

bumpus<-read.csv("Lab/Data/bumpus.csv",header=T)
bumpus.data<-log(as.matrix(bumpus[,(5:13)])) # matrix of linear measurements
sex<-as.factor(bumpus[,2])
surv<-as.factor(bumpus[,4])
TotalLength<-bumpus.data[,1]
Y<-bumpus.data[,-1]

#__________________________________________________________________________#
#Describe your data
(100*apply(bumpus.data,2,sd)) /apply(bumpus.data,2,mean)  # CV: Coefficient of Variation

cor.bumpus<-cor(bumpus.data)
cor.bumpus
pairs(bumpus.data)
vcv.bumpus<-var(bumpus.data)
vcv.bumpus
var(scale(bumpus.data))
dist(bumpus.data)

#__________________________________________________________________________#
#single factor MANOVA
model1<-lm(bumpus.data~sex)
summary(model1)	#yields a set of univariate analyses

summary(manova(model1))	#does multivariate test (using Pillai's)
summary(manova(model1),test="Wilks")	#does multivariate test (using Wilks)

##### Permutation-MANOVA
procD.lm(bumpus.data~sex,iter=9999) #in geomorph

#__________________________________________________________________________#
#Factorial MANOVA
model2<-lm(bumpus.data~sex*surv)
summary(manova(model2))

# VARIANCE COMPONENTS: SSCP matrices
summary(manova(model2))$SS[[1]]	#SSCP matrices for each term in model
summary(manova(model2))$SS[[2]]
summary(manova(model2))$SS[[3]]
summary(manova(model2))$SS[[4]]

#Permutation Factorial MANOVA
procD.lm(bumpus.data~sex*surv)  #NOTE: RRPP USED HERE

#NOTE compare to: 
procD.lm(bumpus.data~sex*surv,RRPP=FALSE)  #'standard' Y-value permutation: not optimal for factorial models

#Pairwise comparisons based on Distances among LS.means (D_Euclid in this case)
advanced.procD.lm(bumpus.data~sex*surv,~sex+surv,groups=~sex*surv)


#__________________________________________________________________________#
### Multivariate Regression
summary(manova(Y~TotalLength))

procD.lm(Y~TotalLength)  #permutation-regression

### Visualizing multivariate regression with summary variables
#Below is the 'by hand' version: functions in geomorph

# 1: Visualize Regression PC1
PC1<-prcomp(Y)$x[,1]
plot(TotalLength, PC1, pch=21, cex=1.25, bg="black")

# 2: Visualize with Regression scores: Drake and Klingenberg 2008
B<-(coef(lm(Y~TotalLength)))
s<-Y%*%B[2,]%*%sqrt(solve(t(B[2,])%*%B[2,])) #slope portion only
plot(TotalLength, s, pch=21, cex=1.25, bg="black")

# 3: Stylized Regression lines: Adams and Nistri 2010
yhat.mancova<-predict(lm(Y~TotalLength))
P<-prcomp(yhat.mancova)$x[,1] #1st dim of predicteds: Adams and Nistri 2010
##NOTE: Adams and Nistri had multiple groups and used: shape~group*size)
plot(TotalLength, P, pch=21, cex=1.5, bg="black")

#__________________________________________________________________________#
#MANCOVA	
summary(manova(lm(Y~TotalLength*sex*surv)))
# FIT COMMON SLOPE
summary(manova(lm(Y~TotalLength+sex*surv)))

#permutation-MANCOVA
procD.lm(Y~TotalLength*sex)

#another implementation with pairwise slope tests included
allom.res<-procD.allometry(Y~TotalLength, ~sex,logsz=FALSE)
summary(allom.res)

### Visualizing MANCOVA
plot(allom.res,method="PredLine")
plot(allom.res,method="RegScore")

#__________________________________________________________________________#
###  Distance Matrix MANOVA: where Y-data input as distance matrix
gdf<-geomorph.data.frame(morph.dist=dist(bumpus.data),sex=sex,surv=surv)
procD.lm(morph.dist~sex*surv,data=gdf) #note: same SS, MS, etc. as above with raw data!


## USEFUL for MANOVA for non-Euclidean data
data(sipoo)  #birds on island data from vegan
dist.sipoo<-vegdist(sipoo, method="jaccard")
sipoo.area <-  c(1.1, 2.1, 2.2, 3.1, 3.5, 5.8, 6, 6.1, 6.5, 11.4, 13,
                 14.5, 16.1 ,17.5, 28.7, 40.5, 104.5, 233) 
names(sipoo.area)<-row.names(sipoo)
island.sz<-ifelse(sipoo.area>10,1,0); island.sz<-as.factor(island.sz)
sipoo.gdf<-geomorph.data.frame(dist.sipoo=dist.sipoo,island.sz=island.sz,sipoo.area=sipoo.area)

#MANOVA comparing community composition in small-large islands
procD.lm(dist.sipoo~island.sz,data=sipoo.gdf)

#### NOTE: the vegan function 'adonis' also performs ANOVA/REGRESSION, BUT ONLY using full permutation
