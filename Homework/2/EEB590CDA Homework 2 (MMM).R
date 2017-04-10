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

#Read in data
dat <- read.csv("Homework/2/HW2.dat1.csv",header=T)
dat1 <- as.matrix(dat[,(3:9)])
X3 <- as.factor(dat1[,1])

#Load libraries

#-----------------------------------------------#
#Describe data
cor <- cor(dat1)
cor
pairs(dat1)
var <- var(dat1)
var

#-----------------------------------------------#
#Single factor MANOVA
datmodel <- lm(dat1)
summary(datmodel)
summary(manova(datmodel))
summary(manova(datmodel), test = "Wilks")









