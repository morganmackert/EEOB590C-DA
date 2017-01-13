# Various regression models

    #packages: geomorph, lmodel2, scatterplot3d

rm(list=ls())
bumpus<-read.csv("Data/bumpus.csv",header=T)
bumpus.data<-log(as.matrix(bumpus[,(5:13)])) # matrix of log-linear measurements
sex<-as.factor(bumpus[,2])
surv<-as.factor(bumpus[,4])
SexBySurv<-as.factor(paste(sex,surv))
Y<-as.matrix(bumpus.data[,1])
X1<-bumpus.data[,2]
X2<-bumpus.data[,3]

#__________________________________________________________________________#
#simple linear regression  (Model I Regression)
lm(Y~X1)

#more information
model1<-lm(Y~X1)
summary(model1)	#provides regression coefficients	
anova(model1)	#provides model term tests
plot(X1,Y,pch=21, bg="black", cex=2)
abline(model1,lwd=2,col="red")


##### Distance-based MANOVA + randomization
library(geomorph)
procD.lm(Y~X1)

#Additional Diagnostic plots
plot(model1)
	#plot 1: fitted vs. residuals (structure in plot BAD)
	#plot 2: residual quantiles (near line GOOD)
	#plot 3: fitted vs. sqrt resid (in a triangle BAD)
	#plot 4: cooks distance (identify points with high leverage)

#__________________________________________________________________________#
#### Model II Regression 
#install.packages("lmodel2")
library(lmodel2)
lmodel2(Y~X1,nperm=999)
  RMA<-lmodel2(Y~X1)
plot(RMA, pch=21,cex=2, bg="black")
abline(model1,lwd=2,col="blue")

#__________________________________________________________________________#
#multiple regression
summary(lm(Y~X1+X2))
anova(lm(Y~X1+X2))

# Plot for multiple regression
#install.packages("scatterplot3d")
library(scatterplot3d)
plot<-scatterplot3d(X1,X2,Y)
plot$plane3d(lm(Y~X1+X2))

#__________________________________________________________________________##polynomial regression
summary(lm(Y~X1+I(X1^2)))	#NOTE: 'I' identity matrix and 'allows' function sqr to be used 
anova(lm(Y~X1+I(X1^2)))

plot(X1,Y)
abline(model1,col="red")
  z<-(order(X1))  
  X1new<-X1[z]
  Ynew<-Y[z]
  polymodel<-predict(lm(Ynew~X1new+I(X1new^2)))
lines(X1new,polymodel)

#__________________________________________________________________________#
###ANCOVA
anova(lm(Y~X2*SexBySurv))

#Implies slopes for M and F statistically equivalent, so drop and re-run
anova(lm(Y~X2+SexBySurv))
model.ancova<-lm(Y~X2+SexBySurv)

colo<-rep("blue",nrow(bumpus.data)); colo[which(SexBySurv == 'f TRUE')]<-"red";
	colo[which(SexBySurv == 'f FALSE')]<-"pink"; colo[which(SexBySurv == 'm FALSE')]<-"lightblue";

plot(X2,Y,pch=21,cex=2,bg=colo)
abline(model.ancova$coefficients[1],model.ancova$coefficients[2])
abline((model.ancova$coefficients[1]+model.ancova$coefficients[3]),model.ancova$coefficients[2])
abline((model.ancova$coefficients[1]+model.ancova$coefficients[4]),model.ancova$coefficients[2])
abline((model.ancova$coefficients[1]+model.ancova$coefficients[5]),model.ancova$coefficients[2])
	#note: parameters are intercept, slope, deviations of int. for groups

#__________________________________________________________________________#
## SIMPLE RESAMPLING EXAMPLE
F.obs<-anova(lm(Y~X1))$F[[1]]  #Find Test value and save
permute<-999
P.Ftest<-1
F.rand.vec<-NULL
F.rand.vec<-rbind(F.rand.vec,F.obs)
for(i in 1:permute){
  ###Shuffle Data
	y.rand<-sample(Y)	#Resample vector 
	F.rand<-anova(lm(y.rand~X1))$F[[1]]  
	F.rand.vec<-rbind(F.rand.vec,F.rand)
	P.Ftest<-if(F.rand>=F.obs) (P.Ftest+1) else P.Ftest
}  
P.Ftest<-P.Ftest/(permute+1)
P.Ftest
####Plot
hist(F.rand.vec,40,freq=T,col="gray")
segments(F.obs, 0, F.obs, 50)  ##Plot Observed value

rm(list=ls())
