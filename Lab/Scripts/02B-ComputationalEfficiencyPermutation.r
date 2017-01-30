##  Some Initial Comments on Computational Efficiency

#As one improves in their R-skills, computational efficiency becomes more critical.  In R one can 
#find numerous ways to do the same thing, but these can differ wildly in how long they take to complete. 

#A useful way to compare implementations is with system.time, or better still microbenchmark (the latter is
   #much more precise and can be used to isolate the slow lines of code in a function you are writing)

       #Packages: microbenchmark

#For example: 
library(microbenchmark)

mymean <- function(x){
  n <- length(x)
  tmp <- 0
  for (i in 1:n){
    tmp <- tmp+x[i]
  }
  mn <- tmp/n
  return(mn)
}
x<-rnorm(100)
n<-length(x)

microbenchmark(mean(x),mymean(x),sum(x)/length(x))
#NOTES on this simple experiment
 #1: 'mymean' is very slow, so while correct this is not a great implementation for obtaining the mean!   
  #2: the base function 'mean' is not the fastest either.  using (sum(x)/length(x)) is 10X faster!
      #NOTE: some reasons why: may functions need to do data I/O checks first (e.g., dimensionality).
     #This slows R down.  But if you know the these inputs, you can rewrite and speed things up.



######A larger example: Different loops, pulling out computations, and lapply/map functions
library(microbenchmark)
quick.pls <- geomorph:::quick.pls
perm.index <- geomorph:::perm.index
n <- 50
p <- 50

a <- factor(round(runif(n)))
b <- factor(round(runif(n)))
form <- ~a*b
Terms <- terms(form)

# Example sub-design matrices
Xa <- model.matrix(Terms[1])
Xb <- model.matrix(Terms[1:2])
Xab <- model.matrix(Terms[1:3])
Xnull <- as.matrix(Xa[,1])

iter <- 999
ind <- perm.index(n, iter)
Y <- matrix(rnorm(p*n), n, p)

# loop versions of lm.fit for all models.  Must be universal for any 
# number of X matrices

# First, a bad but intuitively reasonable for loop
Result <- list(iter+1)
SSEnull <- sum(lm.fit(Xnull, Y)$residuals^2)
k <- ncol(attr(Terms, "factors"))

for(i in 1:(iter+1)){
  SSE.result <- array(1:k)
  y <- Y[ind[[i]],]
  for(j in 1:k){
    X <- model.matrix(Terms[1:j])
    fit<- lm.fit(X,y)
    SSE <- sum(fit$residuals^2)
    SSE.result[j] <- SSE
  }
  Result[[i]] <- c(SSEnull,SSE.result)
}

Result <- simplify2array(Result)
Result[,1:10]


# Second, the same for loop with one important change
# before the loop...
Xs <- list(Xnull = Xnull, Xa = Xa, Xb = Xb, Xab = Xab) # this is the change

Result <- list(iter+1)

for(i in 1:(iter+1)){
  SSE.result <- array(1:length(Xs))
  y <- Y[ind[[i]],]
  for(j in 1:length(Xs)){
    X <- Xs[[j]]# this is the change followed through
    fit<- lm.fit(X,y)
    SSE <- sum(fit$residuals^2)
    SSE.result[j] <- SSE
  }
  Result[[i]] <- SSE.result
}

Result <- simplify2array(Result)
Result[,1:10]


# Map version (uses the logic of the second for loop without so much
# cumbersome code.)

SSEs <- function(Y, Xs) Map(function(x) sum(lm.fit(x,Y)$residuals^2), Xs)

Result <- simplify2array(Map(function(j) SSEs(Y[j,], Xs), ind))
Result[,1:10]

### Time comparison --------------------------------------

microbenchmark(
  { # loop 1
    Result <- list(iter+1)
    SSEnull <- sum(lm.fit(Xnull, Y)$residuals^2)
    k <- ncol(attr(Terms, "factors"))
    
    for(i in 1:(iter+1)){
      SSE.result <- array(1:k)
      y <- Y[ind[[i]],]
      for(j in 1:k){
        X <- model.matrix(Terms[1:j])
        fit<- lm.fit(X,y)
        SSE <- sum(fit$residuals^2)
        SSE.result[j] <- SSE
      }
      Result[[i]] <- c(SSEnull,SSE.result)
    }
    
    Result <- simplify2array(Result)
  },  
  
  { # loop 2
    Result <- list(iter+1)
    
    for(i in 1:(iter+1)){
      SSE.result <- array(1:length(Xs))
      y <- Y[ind[[i]],]
      for(j in 1:length(Xs)){
        X <- Xs[[j]]
        fit<- lm.fit(X,y)
        SSE <- sum(fit$residuals^2)
        SSE.result[j] <- SSE
      }
      Result[[i]] <- SSE.result
    }
    Result <- simplify2array(Result)
  },  
  
  Result <- simplify2array(Map(function(j) SSEs(Y[j,], Xs), ind)), 
  
  times=5)

