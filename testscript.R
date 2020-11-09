################################################################################################################
########################## Test script for obtaining the sparse shooting S estimator ###########################
##### "Sparse regression for large data sets with outliers" by Lea Bottmer, Christophe Croux and Ines Wilms ####
################################################################################################################
rm(list=ls())     # Clean memory
# Load libraries
library(Rcpp) 
library(pcaPP) 
library(robustbase) 
library(robustHD) 
library(mvtnorm)

# Source functions
source("sparseshootS.R")
sourceCpp("sparseshootS.cpp")

# Data example
p <- 50 # number of predictors
n <- 100 # sample size
epsilon <- 0.01 # level of cellwise contamination
K <- 0.1 # Percentage of parameters that are non-zero
beta <- c(rep(1, K*p), rep(0,(1-K)*p))
lambda_grid_points <- 100 # number of values in lambda-grid
cutvalue <- 3 # cut-off value for outlier flaging

# Generate data
Sigma <- matrix(NA, p, p)
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j] = 0.5^(abs(i-j))
  }
}
X <- rmvnorm(n, sigma = Sigma)
e <- rnorm(n, sd = 1)
y <- X%*%beta+e
if(epsilon>0){
  ## Contamination
  vec <- as.vector(X)
  percent <- round(epsilon*n*p)
  cell.variable <- sample(1:length(vec), percent, replace = FALSE)
  vec[cell.variable] = rnorm(percent,50,1)
  X = matrix(vec,n,p)
}

# Sparse Shooting S
fitshootsparse <- sparseshooting(x = X, y = y,  wvalue = cutvalue, nlambda = lambda_grid_points)
fitshootsparse$coef

# Shooting S
fitshoot <- shooting(x = X, y = y,  wvalue = cutvalue)
fitshoot$coef