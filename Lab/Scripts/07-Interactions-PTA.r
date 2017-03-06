#####  Interaction Term exploration: Trajectory Analysis	
rm(list=ls())
library(geomorph)

##Interaction terms and sub-components
data(plethodon)
Y.gpa<-gpagen(plethodon$land,print.progress=FALSE)    #GPA-alignment
gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, site = plethodon$site)

# Example of a test of a factor interaction, plus pairwise comparisons 
procD.lm(coords ~ site*species, data = gdf,print.progress=FALSE)
advanced.procD.lm(coords ~ site*species, ~ site + species, groups = ~site*species, 
                  iter=499, data = gdf,print.progress=FALSE)

# Example of a test of a factor interaction, plus pairwise comparisons, 
# accounting for a common allometry  
procD.lm(coords ~ Csize+site*species, data = gdf,print.progress=FALSE)
advanced.procD.lm(coords ~ Csize + site*species, 
                  ~ log(Csize) + site + species, 
                  groups = ~ site*species, iter = 499, 
                  data = gdf,print.progress=FALSE)

## Example of a comparison of slopes
procD.lm(coords ~ Csize+site*species, data = gdf,print.progress=FALSE)
advanced.procD.lm(coords ~ Csize + site*species, 
                  ~ log(Csize) + site + species, 
                  groups = ~ site*species, sslope = ~ log(Csize), angle.type = "deg", iter = 499, 
                  data = gdf,print.progress=FALSE)

# Example of a test of homogeneity of slopes, plus pairwise slopes comparisons
gdf$group <- factor(paste(gdf$species, gdf$site, sep="."))
procD.lm(coords ~ Csize*group, data = gdf,print.progress=FALSE)
advanced.procD.lm(coords ~ log(Csize) + group, 
                  ~ log(Csize) * group, 
                  groups = ~ group, 
                  slope = ~ log(Csize), angle.type = "deg", 
                  iter = 499, data = gdf,print.progress=FALSE)


##### Trajectory Analysis
#1: Estimate trajectories from LS means in factorial mode
data(plethodon) 
Y.gpa <- gpagen(plethodon$land,print.progress=FALSE)   
Y.dist<-dist(two.d.array(Y.gpa$coords))
gdf <- geomorph.data.frame(Y.dist=Y.dist,Y.gpa, species = plethodon$species, site = plethodon$site)

TA <- trajectory.analysis(coords ~ species*site, data=gdf, iter=199,print.progress=FALSE)
summary(TA, angle.type = "deg")
plot(TA)

#Run using distances between objects in Y-data space
TA2 <- trajectory.analysis(Y.dist ~ species*site, data=gdf, iter=199,print.progress=FALSE)
summary(TA2, angle.type = "deg")
plot(TA2)


#2: Add covariate to analysis
TA <- trajectory.analysis(f1 = coords ~ species*site, f2 = ~ Csize, 
                          data=gdf, iter=199,print.progress=FALSE)
summary(TA, angle.type = "deg")
plot(TA)

#3: Compare groups of pre-existing trajectories
data(motionpaths)

gdf <- geomorph.data.frame(trajectories = motionpaths$trajectories,
                           groups = motionpaths$groups)
TA <- trajectory.analysis(f1 = trajectories ~ groups, 
                          traj.pts = 5, data=gdf, iter=199,print.progress=FALSE)
summary(TA)
plot(TA)


