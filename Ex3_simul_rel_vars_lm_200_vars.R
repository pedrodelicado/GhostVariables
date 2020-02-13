# Simulated example 3
# For j=1,...,4{
#   (xj.1,...,xj.pj) normals with sd sigma.j and correlations rho.j
#   zj = xj.1+...+xj.pj
# }
# where rho.1 = rho.3 = 0, so variables in blocks x1 and x3 are independent.
# Moreoven blocks x1, x2, x3 and x4 are independent. 
# So z1, z3, z4, z2 are independent.
#
# Response variable:
# y = beta1*z1 + beta2*z2 + beta3*z1 + beta4*z2 + epsilon2
# with epsilon ~ N(0,sigma.eps^2)
# and beta3=beta4=0, so variables in x3 and x4 are irrelevant for y.
#
# linear model to be fitted:
# y ~ x1.1+...+x1.p1 + x2.1+...+x2.p2 + x3.1+...+x3.p3 + x4.1+...+x4.p4
#

printing <- FALSE # TRUE # 
printing.paper <- FALSE # TRUE # 
set.seed(123456) # to reproduce results

n1 <- 2000 # size of the training sample 
n2 <- 1000 # size of the test sample
p1 <- 50 # number of uncorrelated variables relevant for y
p2 <- 50 # number of correlated variables relevant for y
p3 <- 50 # number of uncorrelated variables irrelevant for y
p4 <- 50 # number of correlated variables irrelevant for y

sigma.1 <- 1 # sd for the p1 variables x1.1,...,x1.p1
sigma.2 <- 1 # sd for the p2 variables x2.1,...,x2.p2
sigma.3 <- 1 # sd for the p3 variables x3.1,...,x3.p3
sigma.4 <- 1 # sd for the p4 variables x4.1,...,x4.p4

sigma.eps <- 1 # residual sd for defining y

# rho.1 = rho.3 = 0
rho.2 <- .95 # correlation between p2 variables
rho.4 <- .95 # correlation between p4 variables

beta1 <- .5 # coef. of z1=x1.1+...+x1.p1
beta2 <- 1  # coef. of z2=x2.1+...+x2.p2
beta3 <- 0  # coef. of variables in X3
beta4 <- 0  # coef. of variables in X4

# Generating the p1 variables
X1 <- sigma.1 * matrix(rnorm((n1+n2)*p1),ncol=p1)
z1 <- apply(X1,1,sum)

# Generating the p2 variables
Sigma.2 <- matrix(rho.2, nrow=p2, ncol=p2)
diag(Sigma.2) <- 1
eig.Sigma.2 <- eigen(Sigma.2)
sqrt.Sigma.2 <- eig.Sigma.2$vectors %*% diag(eig.Sigma.2$values^.5) %*% t(eig.Sigma.2$vectors)
X2 <- sigma.2 * matrix(rnorm((n1+n2)*p2),ncol=p2) %*% sqrt.Sigma.2
z2 <- apply(X2,1,sum)
#z2 <- as.numeric(X2 %*% (1:p2))/p2

# Generating the p3 variables
X3 <- sigma.3 * matrix(rnorm((n1+n2)*p3),ncol=p3)
z3 <- apply(X3,1,sum)

# Generating the p4 variables
Sigma.4 <- matrix(rho.4, nrow=p4, ncol=p4)
diag(Sigma.4) <- 1
eig.Sigma.4 <- eigen(Sigma.4)
sqrt.Sigma.4 <- eig.Sigma.4$vectors %*% diag(eig.Sigma.4$values^.5) %*% t(eig.Sigma.4$vectors)
X4 <- sigma.4 * matrix(rnorm((n1+n2)*p4),ncol=p4) %*% sqrt.Sigma.4
z4 <- apply(X4,1,sum)


# defining the response variable
y <- beta1*z1 + beta2*z2 + beta3*z3 + beta4*z4 + rnorm(n1+n2,sd=sigma.eps)

X <- cbind(X1,X2,X3,X4)
colnames(X) <- c( paste0("x1.",1:p1), paste0("x2.",1:p2) , paste0("x3.",1:p3), paste0("x4.",1:p4) )
yX <- as.data.frame(cbind(y,X))

# Training sample:
tr.sample <- (1:n1)
# Test sample:-tr.sample

# Fitting the linear model
lm.tr <- lm(y ~ ., data=yX, subset = tr.sample)
(sum.lm.tr<-summary(lm.tr))

# Predicting in the test sample
y.hat.ts <- as.numeric( predict(lm.tr,newdata = yX[-tr.sample,]) )

# variable relevance matrix
library(mgcv)
library(ggplot2)
library(grid)
library(maptools)# For pointLabel
source("relev.ghost.var.R")
source("relev.rand.perm.R")

par(mfrow=c(1,1))
relev.ghost.out <- relev.ghost.var(model=lm.tr, 
                             newdata = yX[-tr.sample,],
                             func.model.ghost.var= lm)

matplot(relev.ghost.out$eig.V$vectors[51:100,c(173,179,180)],type="l");abline(h=0,col=8)
matplot(relev.ghost.out$eig.V$vectors[,c(173,179,180)],type="l");abline(h=0,col=8)

if (printing) pdf("Ex_3_200_vars_relev_ghost.pdf",height = 16, width = 12)
plot.relev.ghost.var(relev.ghost.out, n1=n1, resid.var=sum.lm.tr$sigma^2,
                     sum.lm.tr=sum.lm.tr,
                     vars=c(1,2,50,51,99,100,173,179,180))
if (printing) dev.off()

par(mfrow=c(1,1))
if (printing) pdf("Ex_3_200_cluster_ghost.pdf",height = 8, width = 12)
aux <- relevance2cluster(relev.ghost.out$V, sub="",
                  xlab="Relevance by ghost variables",
                  labels=substr(colnames(X),start=2,stop=2),
                  ret.dist = TRUE)
if (printing) dev.off()
#mds.ghost <- cmdscale(d=aux,k=4)
#pairs(mds.ghost,col=rep(c(1,2,3,4),each=50))

relev.perm.out <- relev.rand.perm(model=lm.tr, newdata = yX[-tr.sample,], func.model.ghost.var= lm)
if (printing) pdf("Ex_3_200_vars_relev_perm.pdf",height = 20, width = 12)
plot.relev.rand.perm(relev.perm.out, relev.ghost=relev.ghost.out$relev.ghost,
                     vars=c(1,2,40,45,100,101,102,140,147,160, 199, 200))
if (printing) dev.off()

par(mfrow=c(1,1))
if (printing) pdf("Ex_3_200_cluster_perm.pdf",height = 8, width = 12)
relevance2cluster(relev.perm.out$V, sub="",
                  xlab="Relevance by random permutations",
                  labels=substr(colnames(X),start=2,stop=2))
if (printing) dev.off()

# Figures for the paper
resid.var <- sum.lm.tr$sigma^2
F.transformed <- resid.var*sum.lm.tr$coefficients[-1,3]^2/n1
alpha <- .01
p <- dim(X)[2]
F.critic.transformed <- resid.var*qf(1-alpha,1,n1-p-1)/n1
relev.ghost <- relev.ghost.out$relev.ghost
relev.rp <- relev.perm.out$relev.rp

eig.vals.V <- relev.ghost.out$eig.V$values
eig.vals.V.tilde <- relev.perm.out$eig.V$values

if (printing.paper) pdf("Ex_3_200.pdf", height = 8, width = 8)
op<-par(mfrow=c(3,2))

# 1
plot(F.transformed,relev.ghost,
     xlim=c(0,max(c(F.transformed,relev.ghost))),
     ylim=c(0,max(c(F.transformed,relev.ghost))),
     xlab=expression(paste("F-statistics*",hat(sigma)^2/n[1])), 
     ylab="Relev. by ghost variables",
     col=substr(colnames(X),start=2,stop=2),
     pch=19)
pointLabel(F.transformed,relev.ghost, substr(colnames(X),start=2,stop=2), 
           col=substr(colnames(X),start=2,stop=2))
abline(a=0,b=1,col=8,lwd=2)
abline(v=F.critic.transformed,h=F.critic.transformed,lty=2,col="blue",lwd=2)

# 2
plot(relev.rp,relev.ghost,
     xlim=c(0,max(relev.rp)),
     ylim=c(0,max(relev.ghost)),
     xlab="Relev. by rand.perm.", 
     ylab="Relev. by ghost variables",
     col=substr(colnames(X),start=2,stop=2),
     pch=19)
pointLabel(relev.rp,relev.ghost, substr(colnames(X),start=2,stop=2), 
           col=substr(colnames(X),start=2,stop=2))

# 3
par(xaxp=c(1,p,min(p,5)))
plot(eig.vals.V, ylim=c(0,max(eig.vals.V)),
     main=expression("Eigenvalues of matrix V"),
     ylab="Eigenvalues",type="b")
abline(h=0,col="red",lty=2)

# 4 
plot(eig.vals.V.tilde, ylim=c(0,max(eig.vals.V.tilde)),
     main=expression(paste("Eigenvalues of matrix ",tilde(V))),
     ylab="Eigenvalues",type="b")
abline(h=0,col="red",lty=2)

# 5 
relevance2cluster(relev.ghost.out$V, sub="",
                         xlab="Relevance by ghost variables",
                         labels=substr(colnames(X),start=2,stop=2),
                  axes=FALSE, ylab="")

# 6
relevance2cluster(relev.perm.out$V, sub="",
                  xlab="Relevance by random permutations",
                  labels=substr(colnames(X),start=2,stop=2),
                  axes=FALSE, ylab="")

par(op)
if (printing.paper) dev.off()

