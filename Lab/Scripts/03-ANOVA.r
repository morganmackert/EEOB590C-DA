# Various ANOVA models
   #Packages:   geomorph


rm(list=ls())
bumpus<-read.csv("Data/bumpus.csv",header=T)
bumpus.data<-log(as.matrix(bumpus[,(5:13)])) # matrix of log-linear measurements
sex<-as.factor(bumpus[,2])
surv<-as.factor(bumpus[,4])
SexBySurv<-as.factor(paste(sex,surv))
 TotalLength<-bumpus.data[,1]


#  Some components of the linear model
model1<-lm(TotalLength~sex)
anova(model1)  	#Statistical summary
predict(model1)		#Predicted values
resid(model1)		#Residuals

## Multi-group histogram
sexdata<-split(TotalLength,sex)
hist(sexdata$m,col="blue",freq=TRUE)
hist(sexdata$f,col="red",freq=TRUE,add=TRUE)
#ANOTHER WAY, SHOWING OVERLAP

hist(sexdata$m, col=rgb(0,0,0,1),freq=TRUE)
hist(sexdata$f, col=rgb(0.8,0.8,0.8,0.5), freq=TRUE,add=T)


#__________________________________________________________________________#
#Single-Factor ANOVA
anova(lm(TotalLength~sex))
plot(sex,TotalLength,ylab="Total Length")
  summary(aov(TotalLength~sex))    #another way to do univariate ANOVA

##### Distance-based MANOVA + randomization
library(geomorph) #DCA's package for morphometrics: some functions useful for other data

procD.lm(TotalLength~sex)

### Single-Factor ANOVA (4 groups)
anova(lm(TotalLength~SexBySurv))

#post-hoc tests
pairwise.t.test(TotalLength, SexBySurv, p.adj = "none")  #standard

#distance-based with randomization
advanced.procD.lm(TotalLength~SexBySurv, ~ 1, groups = ~SexBySurv)

#__________________________________________________________________________#
###Factorial ANOVA
anova(lm(TotalLength~sex+surv+sex:surv))

###SHORTHAND for FULL Factorial
anova(lm(TotalLength~sex*surv))
procD.lm(TotalLength~sex*surv)

###NOTE: pairwise comparison interpretations can change when using different
  #null models
advanced.procD.lm(TotalLength~sex*surv, ~ sex + surv, groups = ~sex*surv)
advanced.procD.lm(TotalLength~sex*surv, ~ 1, groups = ~sex*surv)
#  First uses RRPP for sequential models, 2nd uses grand mean (Y~1)
  #latter pairwise comparisons equivalent to merging factors A & B
advanced.procD.lm(TotalLength~SexBySurv, ~ 1, groups = ~SexBySurv)

#__________________________________________________________________________#
## Type I vs Type III SS:   Balanced Design
A<-gl(2,20)
B<-gl(2,1,40)
Y<-rnorm(40)
anova(lm(Y~A*B))
anova(lm(Y~B*A))  #for balanced designs order doesn't matter

### Now for unbalanced designs
##Order of variables matters with Type I SS
anova(lm(TotalLength~sex+surv+sex:surv))
anova(lm(TotalLength~surv+sex+sex:surv))

#NOTE: type III (marginal) SS are the same (found by dropping 1 term from the model)
drop1(aov(TotalLength~sex+surv+sex:surv),.~.)
drop1(aov(TotalLength~surv+sex+sex:surv),.~.)

### BUT, Sum of SS type III != total SSmodel!
anova(lm(TotalLength~paste(sex,surv))) # total SSModel: A, B, A:B
  #type I SS
sum(anova(lm(TotalLength~sex*surv))[[2]][1:3])   
  #type III SS
sum(drop1(aov(TotalLength~sex*surv),.~.)[[2]][2:4]) 

#__________________________________________________________________________#
###Nested ANOVA (survival nested within sex)
anova(lm(TotalLength~sex/surv))
###NOTE: SS correct, but significance tests of effects run against error term.
	#Not correct (sex should be tested against nested effecgt, surv)
  anova(lm(TotalLength~sex*surv))  #NOTE the SSB.nested = SSB + SSA:B

##FUNCTION to Manually obtain F and P for nested ANOVA
	# FOR SIMPLE NESTED MODEL ONLY (where order is known)
nested.anova <-function(anova.output){
  nr <- nrow(anova.output) # number of rows in ANOVA table
  ms <- anova.output[,3] # mean squares
  df <- anova.output[,1] # degrees of freedom
  F <- ms[-nr] / ms[-1] # new F ratios
  pval <- 1-pf(F, df[-nr], df[-1]) # new p-values
  anova.output[-nr,4] <- F; anova.output[-nr,5] <- pval
  anova.output
}

anova(lm(TotalLength~sex/surv))
nested.anova(anova(lm(TotalLength~sex/surv)))

#__________________________________________________________________________#
##test effects by comparing models
model6<-(lm(TotalLength~sex+surv))
model7<-(lm(TotalLength~sex))
  anova(model7,model6,test="F")


#__________________________________________________________________________#
## SIMPLE RESAMPLING EXAMPLE
F.obs<-anova(lm(TotalLength~sex))[1,4]  #Find Test value and save
permute<-999
P.Ftest<-1
F.rand.vec<-NULL
F.rand.vec<-rbind(F.rand.vec,F.obs)
for(i in 1:permute){
  ###Shuffle Data
	y.rand<-sample(TotalLength)  #SHUFFLE FIRST COL
	F.rand<-anova(lm(y.rand~sex))[1,4]  
	F.rand.vec<-rbind(F.rand.vec,F.rand)
	P.Ftest<-if(F.rand>=F.obs) (P.Ftest+1) else P.Ftest
}  
F.obs
P.Ftest/(permute+1)

hist(F.rand.vec,40,freq=T,col="gray")
segments(F.obs, 0, F.obs, 50)  ##Plot Observed value
#__________________________________________________________________________#

rm(list=ls())
