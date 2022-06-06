### Computing the case-variable relevance case, 
### using the Ghost Variables strategy
#
### Version for glmnet models fitted with library 'glmnet'
# Ghost variables are fitted with glmnet
relev_GhVar_glmnet <- function(glmnet_model, X_ts, y_ts, 
                               func.model.ghost.var="glmnet",
                               nfolds=3, nlambda=10){  
  require(glmnet)

  n2 <- dim(X_ts)[1]
  p <- dim(X_ts)[2]

  y.hat.ts <- predict(glmnet_model, newx = X_ts)
  MSPE.ts <- sum((y_ts-y.hat.ts)^2)/n2
  
  A <- matrix(0,nrow=n2, ncol=p)
  colnames(A) <- colnames(as.matrix(as.data.frame(X_ts)))
  GhostX <- A
  if (func.model.ghost.var!="glmnet"){# then "lm"
    for (j in (1:p)){
      xj.hat <- lm(X_ts[,j]~X_ts[,-j])$fitted.values
      xj.true <- X_ts[,j]
      X_ts[,j] <- xj.hat
      y.hat.ts.j <- predict(glmnet_model, newx = X_ts)
      X_ts[,j] <- xj.true
      A[,j] <- y.hat.ts - y.hat.ts.j
      GhostX[,j] <- xj.hat
    }
  }else{
    X_aux <- scale(X_ts)
    attr_X_aux <- attributes(X_aux)
    mx <- attr_X_aux$`scaled:center`
    sx <- attr_X_aux$`scaled:scale`
    for (j in (1:p)){
      xj.hat <- mx[j] + sx[j]*predict(
        cv.glmnet(x=X_aux[,-j], y=X_aux[,j],standardize=FALSE,
                  nfolds=nfolds, nlambda=nlambda),
        newx=X_aux[,-j]
      )
      xj.true <- X_ts[,j]
      X_ts[,j] <- xj.hat
      y.hat.ts.j <- predict(glmnet_model, newx = X_ts)
      X_ts[,j] <- xj.true
      A[,j] <- y.hat.ts - y.hat.ts.j
      GhostX[,j] <- xj.hat
    }
  }
  # V=(1/n2)*t(A)%*%A
  V=(1/n2)*t(A)%*%A / MSPE.ts
  relev.ghost <- diag(V)
  eig.V <- eigen(V)
  return(list(A=A, V=V, GhostX=GhostX, 
              relev.ghost=relev.ghost, 
              eig.V=eig.V,
              y.hat.test=y.hat.ts,
              MSPE.test=MSPE.ts)
  )
}
