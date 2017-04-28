#HOMEWORK 1: Lectures 1-4
#Bryan Juarez, Nick Lyon, Ashley St. Clair, Morgan Mackert

#Question 1_____________________________________________________________________#
rm(list=ls())
library(mvtnorm) #load library to generate data
#Create data set with no relationship
sigma <- matrix(c(1, 0, 0, 1), nrow = 2) #choose covariance matrix (correlation matrix since variance on diagonal = 1)
r.2 <- NULL #set vector to store r.squared for each of 100 datasets
#generate 100 datasets with 100 observations and calculate r.squared each time
for(i in 1:100){ 
	dat <- rmvnorm(1000, rep(0,  nrow(sigma)),  sigma = sigma)
	r.2 <- c(r.2, summary(lm(dat[, 1] ~ dat[, 2]))$r.squared)
}
mean(sqrt(r.2)) #confirm average r.2 is about 0
plot(dat)
Y <-dat[,1]
X <- dat[,2]
#Create data set with postive relationship
sigma.cor <- matrix(c(1, .95, .95, 1), nrow = 2) #same as above except we changed the covariance
r.2.cor <- NULL
for(i in 1:100){
	dat.cor <- rmvnorm(1000, rep(0,  nrow(sigma.cor)), sigma = sigma.cor)
	r.2.cor <- c(r.2.cor, summary(lm(dat.cor[, 1] ~ dat.cor[, 2]))$r.squared)
}
mean(sqrt(r.2.cor)) #confirm average r.squared is about .95
plot(dat.cor)
Y.cor <-dat.cor[,1]
X.cor<-dat.cor[,2]


#First permutation method
lm(Y~X)
model1<-lm(Y~X)
summary(model1)	#provides regression coefficients	
anova(model1)	#provides model term tests
summary(model1)
perm.snail <- function(permute, model){
	Y <- model$model$Y
	X <- model$model$X
	F.obs<-anova(model)$F[[1]]  #Find Test value and save
	F.rand.vec <- NULL
	F.rand.vec[permute + 1] <- F.obs
	for(i in 1:permute){
  	###Shuffle Data
  		Y.rand <- sample(Y)	#Resample vector 
  		F.rand.vec <- c(F.rand.vec, anova(lm(Y.rand ~ X))$F[[1]])  
	}  
	P.Ftest <- length(which(F.rand.vec >= F.rand.vec[permute + 1])) / (permute + 1) 
	return(P.Ftest)
}

perm.snail(999, model1)

#Second permutation method

perm.lightning <- function(permute, model){
	Y <- model$model$Y #get Y and X from model
	X <- model$model$X
	F.obs <- anova(model)$F[[1]] #get observed from model
	rand <- replicate(permute, sample(Y)) #generate permuted Ys
	p.val <- length(which(F.obs > apply(rand, 2, function(x) {anova(lm(x ~ X))$F[[1]]} )) + 1) / (permute + 1) #vectorized loop
	return(p.val) 
	}
	
perm.lightning(999, model1)

microbenchmark(perm.snail(100, model1), perm.lightning(100, model1))




