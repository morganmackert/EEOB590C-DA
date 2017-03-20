###Some Ordination Methods
      #Packages: geomorph [for data], vegan, MASS

rm(list=ls())

########################
#	PCA: Principal Components Analysis
library(geomorph)
   data(plethodon)
   Y.gpa<-gpagen(plethodon$land, print.progress=FALSE)   
   pleth.data<-scale(two.d.array(Y.gpa$coords),scale=FALSE)  #center the data
   gp<-factor(paste(plethodon$species,plethodon$site))
   
pca.pleth<-prcomp(pleth.data)  #uses SVD for computations (mean-centered is default)
summary(pca.pleth) #Last four are redundant; explain negligible variation
PC.scores<-pca.pleth$x
pca.pleth$rotation[,1]  #1st PC axis only

plot(PC.scores,xlab="PC I", ylab="PC II",asp=1,pch=21,bg=gp,cex = 1.5) #asp=1 is very important! Means aspect ratio = 1, without this will distort the plot.

#Eigenvalue analysis via broken stick model
library(vegan)
screeplot(pca.pleth,bstick = TRUE)  #implies 2 PCAs sufficient graphical representation

##Plot of actual (full data set) vs. PC1-2 (reduced data set) distances 
plot(dist(pleth.data),dist(PC.scores[,1:2]))

### PCA via svd
svd.res<-svd(pleth.data)  
svd.res$d^2/sum(svd.res$d^2)   #same % variation per PC axis
  pc.scores.svd<-svd.res$u%*%diag(svd.res$d)  #PCA scores
plot(pc.scores.svd,asp=1)

#### PCA "by hand" via eigen-analysis
vcv.pleth<-var(pleth.data)	#Calculate PC axes
pc.pleth<-eigen(vcv.pleth) 
pc.pleth$values/sum(pc.pleth$values)   #same % variation per PC axis
pc.scores<-pleth.data%*%pc.pleth$vectors	#Projection
plot(pc.scores,xlab="PC I", ylab="PC II",asp=1,pch=21,bg=gp,cex = 1.5) #colored by group
    #note, axes can be reflected and still retain correct information  

#############Biplot with PCA
biplot(pca.pleth)

###########################
# 	PCoA
pleth.dist<-dist(pleth.data)
PCoA<-cmdscale(pleth.dist)   #from vegan
plot(PCoA,pch=21,bg=gp,cex=1.5,asp=1)

#############################
#	NMDS
pleth.dist<-dist(pleth.data)
pleth.nmds <- metaMDS(pleth.dist, autotransform=FALSE, k=2)
    #A nice function. Runs 20 times with different starting points; tries to avoid local optima
plot(pleth.nmds)
plot(pleth.nmds,type='t')

  #Plot of actual to plot distances (often curved)
plot(pleth.dist, dist(scores(pleth.nmds, display='sites'), method='eucl'), xlab='D.obs', ylab='D.plot')

data(dune)  #from vegan
dune
dune.dist<-vegdist(dune)  #default = Bray-Curtis distance
dune.nmds <- metaMDS(dune.dist, autotransform=FALSE, k=2)
#A nice function. Runs 20 times with different starting points; tries to avoid local optima
plot(dune.nmds)
plot(dune.nmds,type='t')

#Plot of actual to plot distances (often curved)
plot(dune.dist, dist(scores(dune.nmds, display='sites'), method='eucl'), xlab='D.obs', ylab='D.plot')


#################
#Correspondence Analysis
dune.cca<-cca(dune)  #from vegan
plot(dune.cca)

#Detrended Correspondence Analysis
dune.dca<-decorana(dune)
plot(dune.dca)
