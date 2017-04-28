####Some Spatial Statistics
   #Packages: spatstat, vegan, spdep, nlme

# Point patterns
library("spatstat")

#1: Ripley's K
data(cells)
plot(cells)
plot(Kest(cells))   #shows K for various models: isotropic, poisson, etc.
E<-envelope(cells,Kest,nsim=100,rank=2)
plot(E) # plot shows that observed are underdispersed at small spatial scales


# Generate Spatially Autocorrelated Data  (one could also read some in!)
lat<-runif(50,0,5); long<-runif(50,0,5)
g<-cbind(lat,long)  #create an XY spatial grid
y <- sqrt(diag(g%*%t(g))) + rnorm(nrow(g))  #a spatially-autocorrelated variable
        # Valuee is associated with distance from origin
t <- sample(rep(c(1,2),length(y)/2))  #2 ecological 'groups' constrained by spatial contingency
sc <- 0.0   #spatial contingency: 0->1
for(i in 1:length(y)){
  z = scale(y,scale=sd(y))
  crit = 1-sc
  crit =c(-crit/2,crit/2)*3
  if(z[i]<=min(crit)) t[i]=1
  if(z[i]>=max(crit)) t[i]=2 
}

plot(g, pch=21, bg=t, cex=y, asp=1, main="Species diversity proprotional to circle size 
Color designates ecological type")

#2: Covariation of Geography and Data
library(vegan)
mantel(dist(g), dist(y), permutations = 9999)  # Mantel association
mantel.partial(dist(t), dist(y), dist(g), permutations = 9999)  #3-way Mantel holding group constant

#3: Spatial Autocorrelation
library(spdep)
W<-tri2nb(g) #weights with Delauney tesselation
moran.test(y,nb2listw(W))   #positive autocorrelation

#4: Account for Spatial Non-Independence
ols.fit <- lm(y~t, x=T)  #spatial proximity not considered
summary(ols.fit)
anova(ols.fit)

# GLS models with different spatial autocorrelations
  # Warning!  False convergeneces possible
t <- as.factor(t)
geo=data.frame(g,t,y)

library(nlme)
gls.fit.exp = gls(y~t, data=geo, correlation=corExp(form=~lat+long)) #exponential
gls.fit.gaus = gls(y~t, data=geo,  correlation=corGaus(form=~lat+long)) #gaussian
gls.fit.spher = gls(y~t, data=geo,  correlation=corSpher(form=~lat+long)) #spherical
gls.fit.lin = gls(y~t, data=geo,  correlation=corLin(form=~lat+long)) #linear

# look at coeffcients: VERY DIFFERENT when spatial non-independence considered
ols.fit
gls.fit.exp 
gls.fit.gaus 
gls.fit.spher 
gls.fit.lin

# model comparisons
AIC(gls.fit.exp, gls.fit.gaus,gls.fit.lin) #Careful when y is multivariate, see Model Selection lecture
   ## Exponential decay of spatial dependence is best model
# Anova
anova(gls.fit.exp)
anova(gls.fit.gaus)
anova(gls.fit.spher)
anova(gls.fit.lin)

#  NOTE: one could also incorporate spatial dependency 'by hand' using Gower's 
  #distance approach to find spatial covariance matrix
#Pseudo-code
## PCOA BY HAND
### Principal Coordinate Analysis
#  D =   #some distance matrix
#  n = #number of points (y)
#  A = -.5*as.matrix(D^2)
#  ones = matrix(1,n)
#  I = diag(1,n)
#  G = (I-1/n*ones%*%(t(ones)))%*%A%*%(I-1/n*ones%*%(t(ones)))

## GLS fit by hand while accounting for spatial covariance
#  B = ginv(t(X)%*%ginv(G^2)%*%X)%*%t(X)%*%ginv(G^2)%*%y

##NOTE: one should use a 'cut-off' beyond which points are considered independent
  #This is what the exponential, gaussian etc. models above do
  #see Cressie's work on spatial stats
  


