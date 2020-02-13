# Simulated example 1
# X1, independent from X2 and X3, that are correlated.
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


set.seed(123456)# for the graphics in the slides
#
n1 <- 2000 # size of the training sample 
n2 <- 1000 # size of the test sample

sigma.1 <- 1 # sd for x1
sigma.2 <- 1 # sd for x2
sigma.3 <- 1 # sd for x3

sigma.eps <- 1 # residual sd for defining y

rho <- .95 # correlation between x2 and x3

beta1 <- 1 # coef. of y=x_1+...+x_{p1}
beta2 <- 1 # coef. of y=x_1+...+x_{p2}
beta3 <- 1 # coef. of y=x_1+...+x_{p2}

# Generating variables x2 and x3
X1 <- sigma.1 * matrix(rnorm(n1+n2),ncol=1)

# Generating variables x2 and x3
Sigma.2 <- matrix(rho, nrow=2, ncol=2)
diag(Sigma.2) <- 1
eig.Sigma.2 <- eigen(Sigma.2)
sqrt.Sigma.2 <- eig.Sigma.2$vectors %*% diag(eig.Sigma.2$values^.5) %*% t(eig.Sigma.2$vectors)
X23 <- matrix(rnorm((n1+n2)*2),ncol=2) %*% sqrt.Sigma.2 %*%diag(c(sigma.2,sigma.3))

X2<-X23[,1]
X3<-X23[,2]

# defining the response variable
y <- beta1*X1 + beta2*X2 + beta3*X3 + rnorm(n1+n2,sd=sigma.eps)

X <- cbind(X1,X2,X3)
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

if (printing) pdf("Ex5_Relev_GH.pdf", width = 8, height = 6) 
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

if (printing) pdf("Ex5_Relev_RP.pdf", width = 8, height = 6) 
plot.relev.rand.perm(relev.rand.out, relev.ghost=relev.ghost.out$relev.ghost,
                     sum.lm.tr=sum.lm.tr)
if (printing) dev.off() 

relevance2cluster(relev.rand.out$V, sub="",xlab="Relevance by random permutations")

#### Sobre el signo de la corr. entre las columnas de la matriz A
#### que son residuos de las regresiones de una var. explicativa
#### sobre las otras:
if (1==0){
  cor(X)
  A<-relev.ghost.out$A
  cor(A)
  lm.2.3 <- lm(X[,2] ~ X[,3])
  lm.3.2 <- lm(X[,3] ~ X[,2])
  
  par(mfrow=c(1,1))
  plot(X[,2], X[,3], col=8)
  abline(a=0,b=1,col=2,lwd=2)
  abline(lm.2.3$coefficients,col=4,lwd=2)
  abline(a=-lm.3.2$coefficients[1]/lm.3.2$coefficients[2],b=1/lm.3.2$coefficients[2],col=3,lwd=2)
  j<-7
  i <- sample(1:dim(X)[1],j)
  points(X[i,2], X[i,3],pch=19,col=1:j,cex=1.5)
  points(lm.2.3$fitted.values[i], 
         lm.3.2$fitted.values[i],pch=15,col=1:j,cex=1.5)
  legend("bottomright",c("Original values", 
                         "Selected original values",
                         "Corresponding fitted values"),
         pch=c(1,19,15),cex=c(1.5,1.5,1.5),col=c(8,2,2))
  
  op <- par(mfrow=c(2,2))
  plot(X[,2], X[,3],col=round(abs(X[,2]-X[,3]))+1)
  plot(X[,2], X[,3],col=round(abs(X[,2]+X[,3]))+1)
  plot(lm.2.3$fitted.values, lm.3.2$fitted.values, col=round(abs(X[,2]-X[,3]))+1)
  plot(lm.2.3$fitted.values, lm.3.2$fitted.values, col=round(abs(X[,2]+X[,3]))+1)
  
  plot(X[,2], X[,3],col=round(abs(X[,2]+X[,3]))+1)
  plot(lm.2.3$residuals, lm.3.2$residuals, col=round(abs(X[,2]+X[,3]))+1)
  plot(X[,2], X[,3],col=round(abs(X[,2]-X[,3]))+1)
  plot(lm.2.3$residuals, lm.3.2$residuals, col=round(abs(X[,2]-X[,3]))+1)
  par(op)
}