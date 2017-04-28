##  ------------------------------------------------------------------------------------------------------------------------------  ##
                                    # Homework 1 -- Lectures 1 - 4
##  ------------------------------------------------------------------------------------------------------------------------------  ##
# This assignment is due prior to class in week 6 (Feb 20).

# Guidelines
  ## You are to self-select and work in groups: 3-4 in a group.
  ## For the assignment below submit one R-script.
  ## Annotations via comments are highly encouraged. 
  ## The script should run!

# Required libraries
library(mvtnorm)

##  ------------------------------------------------------------------------------------------------------------------------------  ##  
                                            # Question 1
##  ------------------------------------------------------------------------------------------------------------------------------  ##
# Select some form of linear model containing a single dependent variable (continuous) and at least 1 independent variable.


# Next, simulate two datasets: the first with no relationship between X and Y
cor1 <- 0
data1 <- rmvnorm(n = 100, mean = c(4, 10), sigma = matrix(c(1, cor1, cor1, 1), 2, 2) )
plot(data1, xlab = "X", ylab = "Y", main = "Correlation of ~0")
cor(data1)

# and the second with some positive association between X & Y. 
cor2 <- 0.85
data2 <- rmvnorm(n = 100, mean = c(4, 10), sigma = matrix(c(1, cor2, cor2, 1), 2, 2) )
plot(data2, xlab = "X", ylab = "Y", main = "Correlation of ~0.8")
cor(data2)


# Perform 100 simulations under each condition.
  ## Custom function for simulations
simdat <- function(numsims, avg, corel) {
  output <- NULL
  output$length <- 1:numsims
  
  for(i in output$length) {
    require(mvtnorm)
    intermed <- rmvnorm(100, avg, sigma = matrix(c(1, corel, corel, 1), 2, 2) )
    output$stat <- summary(lm(intermed[,2] ~ intermed[,1]))$r.squared
  }
  
  return(as.data.frame(output))
}

# Dataset one (r = 0)
dat1 <- simdat(100, c(4, 10), 0)

# Dataset two (r = 0.85)
dat2 <- simdat(100, c(4, 10), 0.85)


# Run the linear models on all datasets to confirm that on average,
  ## the patterns for condition 1 (no relationship) and condition 2 (some relationship) are met.
  ## (HINT: this requires determining an appropriate summary measure extracted from the linear model). 

# Average R^2 for simulations
mean(dat1$stat)

# Average R^2 for simulations
mean(dat2$stat)

##  ------------------------------------------------------------------------------------------------------------------------------  ##
                                              # Question 2
##  ------------------------------------------------------------------------------------------------------------------------------  ##
# Devise a permutation procedure to evaluate the above linear model.
# Write code for this permutation procedure.
permtest <- function(df, numperms){
  notscram <- df[,1]
  forscram <- df[,2]
  output <- NULL
  output$length <- 1:numperms
  r2.rand <- NULL
  r2.rand[numperms + 1] <- summary(lm(forscram ~ notscram))$r.squared
  
  for(i in output$length){
    stopgap <- sample(forscram, replace = F)
    r2.rand[i] <- summary(lm(stopgap ~ notscram))$r.squared
  }

  return(hist(r2.rand, main = "Histogram of R^2 Values", xlab = "R^2"))
}

# permutations for dataset1 (correlation = ~ 0)
permtest(data1, 1000); abline(v = as.numeric(summary(lm(data1[,2]~data1[,1]))$r.squared))

# permutations for dataset 2 (correlation = ~0.8)
permtest(data2, 1000); abline(v = as.numeric(summary(lm(data2[,2]~data2[,1]))$r.squared))


# Next, devise a SECOND implementation of the same permutation procedure (ie, code the procedure in a different manner).



##########################
# Y'all are on your own for this one, I used up all my bright ideas in the first method
##########################


# For a single dataset compare the two implementations for their computational performance.
library(microbenchmark)

microbenchmark(permtest(data1, 5))
  # permtest = kind gross, but functional, need to specify dataframe and number of permutations a priori

  ## Summarize your findings via comments in the code (e.g., which approach was faster?  Any thoughts as to why?). 

