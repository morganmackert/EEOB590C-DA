# Some Phylogenetic Comparative Methods: PCMs
   ## packages: nlme, ape, geiger, phytools, geomorph

library(nlme)	# Contains GLS 
library(ape)	# Many phylogenetics-related analyses (see Paradis, Analysis of Phylogenetics and Evolution with R)
library(geiger) # Many evolutionary analyses
library(phytools) 
library(geomorph) #contains multivariate PCMs

## Read data, phylogeny, and match the two
#NOTE: a critical issue in R is that the species names in the data matrix [names(Y) or rownames(Y) 
    #depending on input type] match those on the phylogeny [phy$tip.label]: 
anolis<-read.table("Data/Lab-12-anolis.SVL.txt", sep=",", header=T)
anolis<-na.exclude(anolis)       
anolis[,5:6]<-log10(anolis[,5:6])
    head(anolis)
tree<-read.tree("Data/Lab-12-anolis.tree.tre")  #read tree
plot(tree,  cex=0.5)

# Prune tree to match taxa
droplist<-setdiff(tree$tip, rownames(anolis))
tree<-drop.tip(tree,droplist)	
plot(tree, cex=0.7) # dev.off()  #for fixing plot window if needed
#NOTE: the function 'treedata' in geiger prunes both the data matrix and the tree to match

anolis<-anolis[tree$tip,] #re-order data to match tree order (necessary for some functions)
ml<-anolis$Male_SVL; names(ml)<-rownames(anolis)
fml<-anolis$Female_SVL; names(fml)<-rownames(anolis)
ecomorph<-as.factor(anolis$geo_ecomorph); names(ecomorph)<-rownames(anolis)
Y<-cbind(anolis$Male_SVL,anolis$Female_SVL);rownames(Y)<-rownames(anolis)
gdf<-geomorph.data.frame(ml=ml,fml=fml,ecomorph=ecomorph,Y=Y,tree=tree)

  #another dataset (multivariate)
data(plethspecies) 
Y.gpa<-gpagen(plethspecies$land, print.progress = FALSE)    #GPA-alignment
land.gps<-c("A","A","A","A","A","B","B","B","B","B","B") 
gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
names(gp.end)<-plethspecies$phy$tip

###################
### Analyses
   #Regression: non-phylogenetic
anova(lm(anolis$Female_SVL~anolis$Male_SVL))  
summary(lm(anolis$Female_SVL~anolis$Male_SVL))
plot(anolis$Female_SVL~anolis$Male_SVL)
abline(lm(anolis$Female_SVL~anolis$Male_SVL))

summary(gls(Female_SVL~Male_SVL, data=anolis))  #NOTE: non-phylogenetic done another way

# 1: PhylogeneticRegression (PGLS)
   #using GLS
bm.gls<-gls(Female_SVL~Male_SVL, correlation=corBrownian(phy=tree), data=anolis)
summary(bm.gls)  #Here the correlation structure of the phylogeny is used
anova(bm.gls)
  #using D-PGLS: same
pgls.res<-procD.pgls(fml~ml,phy=tree,data=gdf, print.progress = FALSE) 
summary(pgls.res)
  #using Independent Contrasts: same
picM<-pic(anolis$Male_SVL, tree)
picF<-pic(anolis$Female_SVL, tree)
cor.test(picM, picF)
summary(lm(picF~picM - 1)) #Contrasts run through origin: see Garland et al. 1992
plot(picF~picM)
abline(lm(picF~picM - 1))
   #Phylogenetic anova
procD.pgls(fml~ecomorph,phy=tree,data=gdf, print.progress = FALSE)
anova(gls(Female_SVL~geo_ecomorph, correlation=corBrownian(phy=tree), data=anolis))  #same

  #multivariate phy-anova/regression (even when p>N)
procD.pgls(Y~ecomorph,phy=tree,data=gdf, print.progress = FALSE)  #Here Y is multivariate Y~X|phy. 

# 2: Phylogenetic PLS: multivariate
IT<- phylo.integration(Y.gpa$coords,partition.gp=land.gps,phy=plethspecies$phy, print.progress = FALSE)
summary(IT) ; plot(IT)

# 3: Phylogenetic ordination
plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords,ancStates = FALSE)

# 4: Phylogenetic signal
phylosig(tree, anolis$Female_SVL, method="K", test=T, nsim=1000)  #phytools
physignal(gdf$fml,gdf$tree, print.progress = FALSE)   #geomorph
res<-physignal(gdf$Y,gdf$tree)   #multivariate example
summary(res); plot(res)

# 5: Compare evolutionary rates
  #univariate [geomorph]: provide data, tree and groups
compare.evol.rates(A=fml,gp=ecomorph,phy = tree, print.progress = FALSE)
  #multivariate net evolutionary rates [geomorph]
ER<-compare.evol.rates(A=Y.gpa$coords, phy=plethspecies$phy,gp=gp.end, print.progress = FALSE)
summary(ER); plot(ER)
   #NOTE: [phytools] can compare univariate rates, and rate matrices

# 6: Comparing evolutionary models

## Compare BM vs. OU models using GEIGER
geo=get(data(geospiza))      #a smaller dataset (for time reasons)
tmp=treedata(geo$phy, geo$dat)
phy=tmp$phy; dat=tmp$data[,"tarsusL"] 
plot(phy)

bm.M<-fitContinuous(phy, dat, model="BM")
ou.M<-fitContinuous(phy, dat, bounds = list(alpha = c(min = exp(-500), max = exp(2))),model="OU")
  bm.M
  ou.M

# Compare models: LRT & AIC
LRT.M<-(2*(ou.M$opt$lnL-bm.M$opt$lnL))
prob.M<-pchisq(LRT.M, 1, lower.tail=FALSE)
LRT.M
prob.M

bm.M$opt$aic
ou.M$opt$aic #OU does not provide a better fit: use simpler model (BM)

# NOTE: The package OUCH, MVSLOUCH, OUwie, and otherscan compare more complicated models: BM1, BMM, OU1, OUM, etc.).