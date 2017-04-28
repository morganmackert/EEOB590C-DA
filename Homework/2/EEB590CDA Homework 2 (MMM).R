#-----------------------------------------------#
#                   EEB 590C-DA                 #
#                   Homework 2                  #
#                  Lectures 6-10                #
#          Morgan Mackert, Bryan Juarez,        #
#          Ashely St. Clair, Nick Lyon          #
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
library(fpc)

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
summary(princomp(dat1[, 4:9])) #Summarize pc scores
var(dat1[, 4:9]) #Check variances of original variables
PC.scores <- prcomp(dat1[, 4:9])$x #Calculate pc scores for Y data and plot
plot(PC.scores, asp = 1)

#UPGMA
dist <- dist(PC.scores) #Obtain distance matrix for PCs
upgma <- hclust(dist, method = "average") #UPGMA cluster analysis on distance matrix 
plot(as.dendrogram(upgma), horiz = TRUE, lwd = 4)  #Plot cluster analysis dendrogram

#K-means
#Calculate Calinski-Harabasz index within:between group variance index for different number of clusters
CI <- array(dim = c(1, 10)) 
for(i in 2:12){
  kclusters <- kmeans(PC.scores, i)
  CI[i - 1] <- calinhara(dist, kclusters$cluster) 
}
plot(2:12, CI) #Plot index for different number of clusters

plot(PC.scores[, 1:2], col = kmeans(PC.scores, 3)$cluster, asp = 1) #Plot clusters on PCs
prcomp(dat1[, 4:9])$rotation #Check which linear combination of original varibales may help distinguish among groups

#Here we formed clusters based on the PCs of the original Y data in dataset 1. We note that the variance of Y1 is much larger than the variance of any of the other Y variables. Depending on the aim of the study, we might want to scale our variables, or not. However, we proceeded with cluster analysis on these data without scaling it. We then calculated distances between individuals and performed a UPGMA cluster analysis and formed a dendrogram. Visual examination of dendrogram shows that 2-5 clusters may be useful. We then calculated the Calinski Harabasz index for each step of the cluster analysis (2-12 clusters). This analysis showed that 3, 5, 6, and, particularly 8 clusters result in the highest within:between group variance. For simplicity, here we analyze the 3-cluster solution. Labeling points in PC space according to group membership in the 3-cluster solution shows that the difference between groups is well-summarized by differences in PC1, but not so well by differences in PC2. We found that this pattern is due to the fact that Y1 loads highly on PC1 whereas the rest of the variables do not. Therefore, PC1 is a contrast between Y1 and all other variables. In summary, the simplest interpretation of this dataset is that there are about 3 groups of observations in our data which are easily identifiable using PC1. Should we want to analyze this dataset while accounting for the large variance in Y1, we must simply rescale our original variables or even perform principal components analysis on the correlation matrix of our variables.

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