# Simulated example 2
# X1, independent from X3, 
# X2 correlated with X1 and X3
# Y ~ X1 + X2 + X3
#
# linear model to be fitted:
# y ~ x1+...+x_{p1} + epsilon
#
library(mgcv)
library(ggplot2)
library(grid)
library(maptools)# For pointLabel



source("relev.ghost.var.R")
source("relev.rand.perm.R")


printing <- TRUE # FALSE # 

set.seed(123456)
#
n1 <- 2000 # size of the training sample 
n2 <- 1000 # size of the test sample

sigma.1 <- 1 # sd for x1
sigma.2 <- 1 # sd for x2
sigma.3 <- 1 # sd for x3

sigma.eps <- 1 # residual sd for defining y

rho <- .6 # correlation between x2,x3 and X2,X1: must be rho <= sqrt(2)/2= 0.7071068

beta1 <- 1 # coef. of y=x_1+...+x_{p1}
beta2 <- 1 # coef. of y=x_1+...+x_{p2}
beta3 <- 1 # coef. of y=x_1+...+x_{p2}

# Generating variables x2 and x3
X1 <- sigma.1 * matrix(rnorm(n1+n2),ncol=1)

# Generating variables x2 and x3
Sigma <- matrix(rho, nrow=3, ncol=3)
diag(Sigma) <- 1
Sigma[1,3] <- 0
Sigma[3,1] <- 0
eig.Sigma <- eigen(Sigma)
sqrt.Sigma <- eig.Sigma$vectors %*% diag(eig.Sigma$values^.5) %*% t(eig.Sigma$vectors)
X <- matrix(rnorm((n1+n2)*3),ncol=3) %*% sqrt.Sigma %*%diag(c(sigma.1,sigma.2,sigma.3))
X1<-X[,1]
X2<-X[,2]
X3<-X[,3]
# defining the response variable
y <- beta1*X1 + beta2*X2 + beta3*X3 + rnorm(n1+n2,sd=sigma.eps)
colnames(X) <- paste0("x",1:3)
yX <- as.data.frame(cbind(y,X))
colnames(yX) <- c("y",paste0("x",1:3))

# Training sample:
tr.sample <- (1:n1)
# Test sample:-tr.sample

# Fitting the linear model
lm.tr <- lm(y ~ ., data=yX, subset = tr.sample)
(sum.lm.tr <- summary(lm.tr))

# Predicting in the test sample
y.hat.ts <- as.numeric( predict(lm.tr,newdata = yX[-tr.sample,]) )

# variable relevance matrix 
# by ghost variables

relev.ghost.out <- relev.ghost.var(model=lm.tr, 
                                   newdata = yX[-tr.sample,], 
                                   func.model.ghost.var= lm)

if (printing) pdf("Ex6_Relev_GH.pdf", width = 8, height = 6) 
plot.relev.ghost.var(relev.ghost.out, n1=n1, resid.var=sum.lm.tr$sigma^2,
                     sum.lm.tr=sum.lm.tr)
if (printing) dev.off()

relevance2cluster(relev.ghost.out$V, sub="",xlab="Relevance by ghost variables")

#######
# variable relevance matrix 
# by random permutation

relev.rand.out <- relev.rand.perm(model=lm.tr, 
                                  newdata = yX[-tr.sample,], 
                                  func.model.ghost.var= lm)

if (printing) pdf("Ex6_Relev_RP.pdf", width = 8, height = 6) 
plot.relev.rand.perm(relev.rand.out, relev.ghost=relev.ghost.out$relev.ghost,
                     sum.lm.tr=sum.lm.tr)
if (printing) dev.off() 

relevance2cluster(relev.rand.out$V, sub="",xlab="Relevance by random permutations")