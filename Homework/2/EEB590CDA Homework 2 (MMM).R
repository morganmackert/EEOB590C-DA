#-----------------------------------------------#
#                   EEB 590C-DA                 #
#                   Homework 2                  #
#                  Lectures 6-10                #
#-----------------------------------------------#

#Assignment: 
#1: Select one of the two datasets (HW2.dat1.csv or HW2.dat2.csv).  Each contain a multivariate dataset and several independent (X) variables.  Using the methods learned in weeks 6-10, examine patterns in the dataset. You may use one or more (or all) of the X-variables, and a variety of methods to describe the patterns.

#You must use at least one method from the material learned in: 
#Weeks 6-7
#Week 8
#Week 9
#Week 10  

#USE COMMENTS IN THE R CODE to describe what the patterns you find represent.

#-----------------------------------------------#
#                    Week 6                     #
# Multivarariate Data and General Linear Models #
#-----------------------------------------------#

#Clear environment
rm(list=ls())

#Load libraries
library(ggcorrplot)
library(geomorph)
library(vegan)

#Read in data
dat1 <- read.csv("Homework/2/HW2.dat1.csv",header=T)
dat1.data <- as.matrix(dat1[,(4:9)])
X1<-as.factor(dat1[,1])
X2<-as.factor(dat1[,2])
X3<-dat1.data[,3]
Y1<-dat1.data[,1]
Y2<-dat1.data[,-1]

#-----------------------------------------------#
#Describe data
(100*apply(dat1.data,2,sd)) /apply(dat1.data,2,mean)
#Gives coefficient of variation; the variation relative to the mean of the trait.

cor <- cor(dat1.data)
cor
ggcorrplot(cor)
vcv.dat1 <- var(dat1.data)
vcv.dat1
var(scale(dat1.data))
dist(dat1.data)
dist(dat1.data[1:3,])

#-----------------------------------------------#
#Permutation MANOVA
procD.lm(dat1.data ~ X1,iter=9999)

#-----------------------------------------------#
#Factorial MANOVA
model2 <- lm(dat1.data ~ X1*X2)
summary(model2)
summary(manova(model2))

#-----------------------------------------------#
#                    Week 8                     #
#                  Ordination                   #
#-----------------------------------------------#
pca.dat1 <- prcomp(dat1.data) #uses SVD for computations (mean-centered is default)
summary(pca.dat1)
PC.scores <- pca.dat1$x
pca.dat1$rotation[,1] #1st PC axis only
plot(PC.scores, xlab = "PC I", ylab = "PC II", asp = 1, pch = 21, cex = 1.5)

#-----------------------------------------------#
#Eigenvalue analysis via broken stick model
screeplot(pca.dat1,bstick = TRUE) #Implies 2 PCAs sufficient graphical representation
plot(dist(dat1.data),dist(PC.scores[,1:2])) #We see teardrop spread

#-----------------------------------------------#
#Biplot with PCA
biplot(pca.dat1)

#-----------------------------------------------#
#                    Week 9                     #
#                  Clustering                   #
#-----------------------------------------------#
#PCoA of dat1 data
#dat2 <- dat1[-c(10:136), ]
dat1.dist <- dist(dat1)
dat1.PCoA <- cmdscale(dat1.dist)
plot(dat1.PCoA, pch=21, bg= 'black', cex=1.5, asp=1)
text(PCoA[,1] + 1.5, PCoA[,2], row.names(dat1.dist))

#Clustering methods
dat1.single <- hclust(dat1.dist, method = "single")       #Single-link
dat1.complete <- hclust(dat1.dist, method = "complete")   #Complete-link
dat1.upgma <- hclust(dat1.dist, method = "average")       #UPGMA = average-link
dat1.upgmc <- hclust(dat1.dist, method = "centroid")      #UPGMC
dat1.wpgma <- hclust(dat1.dist, method = "mcquitty")      #WPGMA
dat1.wpgmc <- hclust(dat1.dist, method = "median")        #WPGMC
dat1.wards <- hclust(dat1.dist, method = "ward.D")         #Ward's

#Plot different clustering methods
plot(dat1.single, hang = -1, lwd = 2)
plot(as.dendrogram(dat1.single),horiz = TRUE, lwd = 4, xlim = c(16,-1))    #single-link
plot(as.dendrogram(dat1.complete), horiz = TRUE, lwd = 4, xlim = c(16,-1))  #complete-link
plot(as.dendrogram(dat1.upgma), horiz = TRUE, lwd = 4, xlim = c(16,-1))     #UPGMA
plot(as.dendrogram(dat1.wpgma), horiz = TRUE, lwd = 4, xlim = c(16,-1))     #WPGMA
plot(as.dendrogram(dat1.upgmc),horiz=TRUE,lwd=4,xlim=c(16,-1))     #UPGMC
plot(as.dendrogram(dat1.wpgmc),horiz=TRUE,lwd=4,xlim=c(16,-1))     #WPGMC
plot(as.dendrogram(dat1.wards),horiz=TRUE,lwd=4,xlim=c(16,-1))     #Ward's

plot(dat1.dist,cophenetic(dat1.upgma))

#-----------------------------------------------#
#                   Week 10                     #
#       Multivarariate Association and          #
#            Canonical Ordination               #
#-----------------------------------------------#
#From perANOVA, we know that  X3 and X1 significantly affect Y1
procD.lm(Y1 ~ X1 + X2 + X3, data = dat1)

#It is of value to consider how those variables are associated with Y1, and such a question can be asked via partial least squares analysis
dat.pls <- two.b.pls(dat1$X3, dat1$Y1)
summary(dat.pls)
plot(dat.pls)
#This analysis show that X3 and Y1 are significantly correlated

#It is now of value to consider if the association of Y1 with X3 varies between the two levels of X1 when X1 = 0, X3 and Y1 are significantly correlated
dat.x1.0 <- subset(dat1, dat1$X1 == 0)
x1.0.pls <- two.b.pls(dat.x1.0$X3, dat.x1.0$Y1)
summary(x1.0.pls)

#When X1 = 1, X3 and Y1 are also correlated
dat.x1.1 <- subset(dat1, dat1$X1 == 1)
x1.1.pls <- two.b.pls(dat.x1.1$X3, dat.x1.1$Y1)
summary(x1.1.pls)

plot(dat.pls) #Association of Y1 with X3 (not considering X1)
plot(x1.0.pls) #Association of Y1 with X3 (X1 = 0)
plot(x1.1.pls) #Association of Y1 with X3 (X1 = 1)

#perANOVA further supports that the relationship of X3 with Y1 does not vary as a function of X1
procD.lm(Y1 ~ X3 * X1, data = dat1)
#Nonsignificant interaction term























