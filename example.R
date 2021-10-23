
library(nortest)
library(MASS)
library(splines)
n <- 200; p <- 10; rho <- 0.5
sigma <- toeplitz(rho^seq(0, p-1))
X <- mvrnorm(n, mu=rep(0, p), sigma) 
error <- rnorm(n)
#example1
y <- X[,1] + X[,2] + error
condsdr(X, y, trans=1, seed=1, h=10)
x <- X[,1]; z <- X[,-1]
condcor(x, y, z, trans=1, seed=1, h=10)
#example2
y <- 3*sin(X[,1])+3*sin(X[,p])+0.1*error
condsdr(X, y, trans=1, seed=1, h=10)
x <- X[,1]; z <- X[,-1]
condcor(x, y, z, trans=1, seed=1, h=10)
x <- X[,3]; z <- X[,-3]
condcor(x, y, z, trans=1, seed=1, h=10)
#example3
y=sign(X[,1]+X[,p])*exp(X[,2]+X[,p-1])+0.1*error
condsdr(X, y, trans=1, seed=1, h=10)
x <- X[,1]; z <- X[,-1]
condcor(x, y, z, trans=1, seed=1, h=10)
x <- X[,3]; z <- X[,-3]
condcor(x, y, z, trans=1, seed=1, h=10)

