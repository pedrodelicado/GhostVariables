# Variation of simul_rel_vars_Ejemplo_7.R
# - y is generated according to a gam model, and estimated with a gam
# 

# Simulated example
# X1,X2 independent from X3,X4,X5,
# X1 and X2 correlated,
# X3, X4 and X5 correlated,
#

library(mgcv)
library(ggplot2)
library(grid)
library(maptools)# For pointLabel

source("relev.ghost.var.R")
source("relev.rand.perm.R")

printing <-  FALSE # TRUE #FALSE # 
set.seed(1234)
#
n1 <- 2000 # size of the training sample 
n2 <- 1000 # size of the test sample

rho.12  <- .95 # correlation between x1 and x2
rho.345 <- .7 # correlation between x3, X4 and x5

sigma.1 <- 1 # sd for x1
sigma.2 <- 1 # sd for x2
sigma.3 <- 1 # sd for x3
sigma.4 <- 1 # sd for x4
sigma.5 <- 1 # sd for x5

sigma.eps <- 1/2 # residual sd for defining y

beta1 <- 1 # coef. of y=x_1+...+x_{p1}
beta2 <- 0.5#0#0.5 #.25 # entre -1 i 1. Marca la quantitat d'aditivitat del terme en (X2,X3): 
             # quan m?s lluny de 0, m?s aditivitat
             # 0, .25, 1
            # For beta2=0, do beta2.2=1
            # For beta2==.5, beta2.2=1
beta2.2 <- 1
#if (beta2==0) beta2.2 <- 2/3
beta3 <- 1 # coef. of y=x_1+...+x_{p2}
beta4 <- 1 # coef. of y=x_1+...+x_{p2}
beta5 <- 1 # coef. of y=x_1+...+x_{p2}

# Generating variables x1 and x2
Sigma.12 <- matrix(rho.12, nrow=2, ncol=2)
diag(Sigma.12) <- 1
eig.Sigma.12 <- eigen(Sigma.12)
sqrt.Sigma.12 <- eig.Sigma.12$vectors %*% diag(eig.Sigma.12$values^.5) %*% t(eig.Sigma.12$vectors)
X12 <- matrix(rnorm((n1+n2)*2),ncol=2) %*% sqrt.Sigma.12 %*%diag(c(sigma.1,sigma.2))

X1<-X12[,1]
X2<-X12[,2]


# Generating variables x3, x4 and X5
Sigma.345 <- matrix(rho.345, nrow=3, ncol=3)
diag(Sigma.345) <- 1
eig.Sigma.345 <- eigen(Sigma.345)
sqrt.Sigma.345 <- eig.Sigma.345$vectors %*% diag(eig.Sigma.345$values^.5) %*% t(eig.Sigma.345$vectors)
X345 <- matrix(rnorm((n1+n2)*3),ncol=3) %*% sqrt.Sigma.345 %*%diag(c(sigma.3,sigma.4,sigma.5))

X3<-X345[,1]
X4<-X345[,2]
X5<-X345[,3]

# defining the response variable
#y <- beta1*X1 + beta2*X2 + beta3*X3 + beta4*X4 + beta5*X5 + rnorm(n1+n2,sd=sigma.eps)
#y <- beta1*cos(X1) + beta2*cos(X2) + beta3*cos(X3) + beta4*cos(X4) + beta5*cos(X5) + rnorm(n1+n2,sd=sigma.eps)
#y <- beta1*X1 + (beta2*(X2+X3) + sqrt(1-beta2^2)*X2*X3) + beta4*X4 + beta5*X5 + rnorm(n1+n2,sd=sigma.eps)

#y <- beta1*cos(X1) + beta2.2*(beta2*(cos(X2)+cos(X3)) + sqrt(1-beta2^2)*X2*X3) + 
#  beta4*cos(X4) + beta5*cos(X5) + rnorm(n1+n2,sd=sigma.eps)

y <- beta1*sin(X1) + beta2.2*(beta2*(cos(X2)+cos(X3)) + (1-beta2)*X2*X3) + 
  beta4*cos(X4) + beta5*sin(X5) + rnorm(n1+n2,sd=sigma.eps)

X <- cbind(X1,X2,X3,X4,X5)
colnames(X) <- paste0("x",1:5)
yX <- as.data.frame(cbind(y,X))
colnames(yX) <- c("y",paste0("x",1:5))

# Training sample:
tr.sample <- (1:n1)
# Test sample:-tr.sample

# Fitting the linear model
#gam.tr <- gam(y ~ s(x1)+s(x2)+s(x3)+s(x4)+s(x5), data=yX, subset = tr.sample)
gam.tr <- gam(y ~ s(x1)+s(x2,x3)+s(x4)+s(x5), data=yX, subset = tr.sample)
(sum.gam.tr <- summary(gam.tr))

if (printing) pdf(file="Ex7_gam_model_variaciones.pdf", width=8,height=6)
plot(gam.tr,residuals=TRUE,pages=1)
if (printing) dev.off()

# Predicting in the test sample
y.hat.ts <- as.numeric( predict(gam.tr,newdata = yX[-tr.sample,]) )

# variable relevance matrix 
# by ghost variables

relev.ghost.out <- relev.ghost.var(model=gam.tr, newdata = yX[-tr.sample,], func.model.ghost.var= lm)

if (printing) pdf("Ex7_gam_Relev_GH_variaciones.pdf", width = 8, height = 6) 
plot.relev.ghost.var(relev.ghost.out, #resid.var=gam.tr$sig2, 
                     n1=n1)
if (printing) dev.off()


#######
# variable relevance matrix 
# by random permutation

relev.rand.out <- relev.rand.perm(model=gam.tr, newdata = yX[-tr.sample,], func.model.ghost.var= lm)

if (printing) pdf("Ex7_gam_Relev_RP_variaciones.pdf", width = 8, height = 6) 
plot.relev.rand.perm(relev.rand.out, relev.ghost=relev.ghost.out$relev.ghost)
if (printing) dev.off()

