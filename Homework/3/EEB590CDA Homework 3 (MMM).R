#-----------------------------------------------#
#                   EEB 590C-DA                 #
#                   Homework 3                  #
#                 Lectures 11-13                #
#          Morgan Mackert, Bryan Juarez,        #
#          Ashley St. Clair, Nick Lyon          #
#-----------------------------------------------#

#Assignment: 
#Obtain a dataset. This may be one of your own, a dataset from DRYAD, or some other source.  Identify hypotheses for this dataset and analyze patterns in the data. You may use any methods learned during the semester, but at least one analysis must come from material learned in:
#Week 11
#Week 12
#Week 13

#USE COMMENTS IN THE R CODE to describe what the patterns you find represent.

#-----------------------------------------------#
#                   Week 11                     #
#               Model Selection                 #
#-----------------------------------------------#

#Clear environment
rm(list=ls())

#Load libraries and data
library(vegan)
data(varespec)
#Lichen sp with the greatest cover is Cladstel
data(varechem)
#Of the metals on which data is collected, Aluminum is the most abundant
testdata <- cbind(varespec, varechem)
str(testdata)

#Want to test relationship of "important" chemicals with the most abundant lichen cover spp.
#Nitrogen is known to be very important to plant survival and is often an arbiter of which plant species are most successful in a given environment.
#Aluminum however is typically much more rare in plant habitats and due to that rarity can inform species abundances.
#Of additional importance for lichen is the amount of hummus (sort of vegetative matter) that they are growing on.
#This sort of structural component for habitat can have implications for lichen abundances outside of nutrient availability.

#Q: Do common nutrients (i.e. N), rare nutrients (i.e. Al), or habitat stucture (i.e. hummus depth) have distinct or interactive effects on lichen abundance?

#Need to include the null model just to see if any of these variables explain more variation in the data than nothing
lm.null <- lm(Cladstel ~ 1 , data = testdata); summary(lm.null)

#H0: N does not vary significantly with Cladstel cover
#HA: N does explain Cladstel cover
lm.n <- lm(Cladstel ~ N, data = testdata); summary(lm.n) 
#N = sig

#H0: Al does not vary significantly with Cladstel cover
#HA: Al does explain Cladstel cover
lm.al <- lm(Cladstel ~ Al, data = testdata); summary(lm.al) 
#Al = sig

#H0: Hummus depth does not vary significantly with Cladstel cover
#HA: Hummus does explain Cladstel cover
lm.hum <- lm(Cladstel ~ Humdepth, data = testdata); summary(lm.hum) 
#Hum = *almost* significant

#Included together?
lm.noint <- lm(Cladstel ~ N + Al + Humdepth, data = testdata); summary(lm.noint)
#N = sig, Al = *almost*, Hum = NS

#With interactive effects?
lm.int <- lm(Cladstel ~ N * Al * Humdepth, data = testdata); summary(lm.int) 
#all = NS

aics <- AIC(lm.null, lm.n, lm.al, lm.hum, lm.noint, lm.int); aics
sort(aics[,2])
#AIC scores
#lm.noint < lm.al < lm.n < lm.hum < lm.null < lm.int
#∆AIC is <4 between lm.noint and lm.al, but > 4 between lm.noint and all other models

#Conclusions:
#The model without interaction terms is the model with the lowest AIC score and should be the selected one.
#The interaction model actually has a slightly higher AIC than the null model and *definitely* should not be used.
#Each of the hypothesized important variables explain more data than simply the null (though not always significant).
#BUT a model informed with all three variables as discrete effects explains the most variation in the data.

#-----------------------------------------------#
#         Principle Components Analysis         #
#                 and Clustering                #
#-----------------------------------------------#

#Load libraries and data
library(geiger)
library(geomorph)
data("geospiza") #Darwins finch data on 5 morphological variable for 14 Species
geospiza #View the data
geo.dat <- geospiza$geospiza.data #Create data file

#PCA: Principal Components Analysis
pca.geo <- prcomp(geo.dat)  #uses SVD for computations (mean-centered is default)
summary(pca.geo) #PC1 explains 90% of variance
PC.scores <- pca.geo$x  
pca.geo$rotation  
plot(PC.scores,xlab = "PC I", ylab = "PC II", asp = 1, pch = 21, bg = 'magenta', cex = 1.5)
#By looking at the rotation matrix of the PCA we are able to see the loadings of all the variables on the PCs. If we do this we find that beakD and gonysW explain the most variation in PC1. BeakD explaining 72% and gonysW explaining 57%. It is no surprise that these two explain variation because beak and gonys are related; they are both aspects of the mandibles and bills of birds. PC2 represents a strong contrast between CulmenL which loads the most negatively and beakD which loads the most positively. The next three PCs also represent contrasts but interpreting the information becomes a lot more hazy. PC3 is a contrast between BeakD which loads most negatively and tarsusL which loads most positively. 

#Eigenvalue analysis via broken stick model
screeplot(pca.geo,bstick = TRUE) #PC1 explains most of the variation in dataset, but can also interpret up to about 3 PCs
plot(dist(geo.dat), dist(PC.scores[,1:2])) #The PC space reliably represents distances between our original objects, there is a strong linear relationship.

#K-means clustering
#Compare TESS of the k-means PC scores
TESS <- array(NA,6)
for (i in 1:6){
  TESS[i]<-kmeans(PC.scores,i)$tot.withinss
}
plot(TESS)
#seems to bottom out at 3 groups so we will use kclusters with 3 groups

#K-means groups plotted on PC space (PC1 vs PC2)
kclusters3 <- kmeans(PC.scores,3)
plot(PC.scores[,1:2],col=kclusters3$cluster)
points(kclusters3$centers, col = 1:3, pch = 8, cex=2)

#-----------------------------------------------#
#                   Week 12                     #
#        Phylogenetic Comparative Methods       #
#-----------------------------------------------#
dat <- treedata(geospiza$geospiza.tree, geospiza$geospiza.data) #Get data, prune tree and data
gdf <- geomorph.data.frame(beak= dat$data[, 4], gonys = dat$data[, 5], wing = dat$data[, 1])
beak.gonys <- procD.pgls(beak~gonys, dat = gdf, phy = dat$phy, rep = 999) #Test whether beak and gonys are related as expected from PCA analysis
beak.gonys[[1]][1, 7] #p-value < 0.001
beak.wing <- procD.pgls(beak~wing, dat = gdf, phy = dat$phy, rep = 999) #Test whether birds with larger wings which can presumably carry larger loads have large beaks which potentially allows them to eat larger prey items.
beak.wing[[1]][1, 7] #p-value < 0.001

#As expected, when acknowledging that species share common ancestry, there is a significant positive relationship between beak length and gonys size, as well as between wing length and beak length. However, these analyses ignore the positive relationship between the relevant variables and body size.