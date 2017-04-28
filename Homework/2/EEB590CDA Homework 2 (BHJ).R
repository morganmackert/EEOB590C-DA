

####Chapter 9
library(fpc) #load library and data
dat1 <- read.csv("Homework/2/HW2.dat1.csv")

summary(princomp(dat1[, 4:9])) #summarize pc scores
var(dat1[, 4:9]) #check variances of original variables
PC.scores <- prcomp(dat1[, 4:9])$x #calculate pc scores for Y data and plot
plot(PC.scores, asp = 1)

##UPGMA
dist <- dist(PC.scores) #obtain distance matrix for PCs
upgma <- hclust(dist, method = "average") #UPGMA cluster analysis on distance matrix 
plot(as.dendrogram(upgma), horiz = TRUE, lwd = 4)  #plot cluster analysis dendrogram

#K-means
#Calculate Calinski-Harabasz index within:between group variance index for different number of clusters
CI <- array(dim = c(1, 10)) 
for(i in 2:12){
	kclusters <- kmeans(PC.scores, i)
	CI[i - 1] <- calinhara(dist, kclusters$cluster) 
}
plot(2:12, CI) #plot index for different number of clusters

plot(PC.scores[, 1:2], col = kmeans(PC.scores, 3)$cluster, asp = 1) #plot clusters on PCs
prcomp(dat1[, 4:9])$rotation #Check which linear combination of original varibales may help distinguish among groups

# Here we formed clusters based on the PCs of the original Y data in dataset 1. We note that the variance of Y1 is much larger than the variance of any of the other Y variables. Depending on the aim of the study, we might want to scale our variables, or not. However, we proceeded with cluster analysis on these data without scaling it. We then calculated distances between individuals and performed a UPGMA cluster analysis and formed a dendrogram. Visual examination of dendrogram shows that 2-5 clusters may be useful. We then calculated the Calinski Harabasz index for each step of the cluster analysis (2-12 clusters). This analysis showed that 3, 5, 6, and, particularly 8 clusters result in the highest within:between group variance. For simplicity, here we analyze the 3-cluster solution. Labeling points in PC space according to group membership in the 3-cluster solution shows that the difference between groups is well-summarized by differences in PC1, but not so well by differences in PC2. We found that this pattern is due to the fact that Y1 loads highly on PC1 whereas the rest of the variables do not. Therefore, PC1 is a contrast between Y1 and all other variables. In summary, the simplest interpretation of this dataset is that there are about 3 groups of observations in our data which are easily identifiable using PC1. Should we want to analyze this dataset while accounting for the large variance in Y1, we must simply rescale our original variables or even perform principal components analysis on the correlation matrix of our variables.