### Computing the case-variable relevance matrix, 
### using the random permutation strategy

relev.rand.perm <- function(model,newdata=model$call$data,
                           nperm=1, 
                           permut.gam.vars =NULL, ...){
  #func.model <- eval(parse(text=class(model)[1])) # What kind of model has been fitted 
  #data.tr <- model$call$data[model$call$subset,] # data used for training the model
  #n <- dim(data.tr)[1]
  attr.model <- attributes(model$terms) #getting the varaible names in the model
  term.labels <- attr.model$term.labels #getting the varaible names in the model
  # re-ordering variables in gam models to recover the order desired by the user 
  if (!is.null(permut.gam.vars)) term.labels <- term.labels[permut.gam.vars]
  
  p <- length(term.labels)

  n2 <- dim(newdata)[1] # newdata is the test sample
  # Response variable in the test sample
  name.y <- attr.model$variables[[2]]
  y.ts <- newdata[[ name.y ]]
  # Predicting in the test sample
  y.hat.ts <- as.numeric( predict(model,newdata = newdata) )
  
  MSPE.ts <- sum((y.ts-y.hat.ts)^2)/n2

  A <- matrix(0,nrow=n2, ncol=p)
  colnames(A) <- term.labels
  diag.cov.X <- numeric(p)
  yX.ts <- newdata
  for (perm in 1:nperm){
    for (j in (1:p)){
      yX.ts.aux <- yX.ts
      xj.perm <- sample(yX.ts[,term.labels[j]])
      diag.cov.X[j] <- var(xj.perm)
      yX.ts.aux[,term.labels[j]] <- xj.perm
      y.hat.ts.j <- as.numeric( predict(model,newdata = yX.ts.aux) )
      A[,j] <- A[,j] + y.hat.ts - y.hat.ts.j
    }
  }
  A <- A/nperm
  # V=(1/n2)*t(A)%*%A
  V=(1/n2)*t(A)%*%A / MSPE.ts
  relev.rp <- diag(V)
  eig.V <- eigen(V)
  return(list(A=A, V=V, 
              relev.rp=relev.rp, 
              eig.V=eig.V,
              y.hat.test=y.hat.ts,
              MSPE.test=MSPE.ts,
              diag.cov.X=diag.cov.X)
         )
}

plot.relev.rand.perm <- function(relev.rand.out,
                                 vars=NULL, sum.lm.tr=NULL, relev.ghost=NULL,
                                 alpha=.001, ncols.plot=3){
  A <- relev.rand.out$A
  V <- relev.rand.out$V
  eig.V <- relev.rand.out$eig.V
  relev.rp <- relev.rand.out$relev.rp
  diag.cov.X <- relev.rand.out$diag.cov.X
  
  p  <- dim(A)[2]
  
  if (ncols.plot<3){
    ncols.plot<-3 
    warning("The number of plot columns must be at least 3")
  }
  max.plots <- 4*ncols.plot
  if (is.null(vars)){
    vars <- 1:min(max.plots,p)
  }else{
    if (length(vars)>max.plots){
      vars <- vars[1,max.plots]
      warning(
        paste("Only the first", max.plots, "selected variables in 'vars' are used"))
    }
  }
  n.vars <- length(vars)
  nrows.plot <- 1 + n.vars%/%ncols.plot + (n.vars%%ncols.plot>0)
  
  if (!is.null(sum.lm.tr)){
    F.transformed <- 2*sum.lm.tr$coefficients[-1,1]^2*diag.cov.X
    #F.critic.transformed <- qf(1-alpha,1,n1-p-1)*2*sum.lm.tr$coefficients[-1,2]^2*diag.cov.X
  #}else{
    #F.critic.transformed <- numeric(p)
  }
  
  rel.Gh <- data.frame(relev.rp=relev.rp)
  rel.Gh$var.names <- colnames(A)
  
  plot.rel.Gh <- ggplot(rel.Gh) +
    #geom_bar(aes(x=reorder(var.names, relev.univ.rp), y=relev.univ.rp),
    geom_bar(aes(x=reorder(var.names,X=length(var.names):1), y=relev.rp), 
             stat="identity", fill="darkgray") +
    ggtitle("Relev. by rand.permut.") +
    theme(axis.title=element_blank())+
    theme_bw()+
    ylab("Relevance")+
    xlab("Variable name") + 
    coord_flip()

  # if (!is.null(sum.lm.tr)){
  #   plot.rel.Gh <- plot.rel.Gh +
  #   geom_hline(aes(yintercept = F.critic.transformed),color="blue",size=1.5,linetype=2)
  # }    
  # plot.rel.Gh <- plot.rel.Gh + coord_flip()
  
  # plot.rel.Gh.pctg <- ggplot(rel.Gh) +
  #   #geom_bar(aes(x=reorder(var.names, relev.univ.rp), y=relev.univ.rp),
  #   geom_bar(aes(x=var.names, y=100*relev.rp/sum(relev.rp)),
  #            stat="identity") +
  #   coord_flip() +
  #   ggtitle("Relev. by rand.permut. (% of total relevance)") +
  #   theme(axis.title=element_blank())
  
  # eigen-structure
  # eig.V <- eigen(V)
  eig.vals.V <- eig.V$values
  eig.vecs.V <- eig.V$vectors
  
  expl.var <- round(100*eig.vals.V/sum(eig.vals.V),2)
  cum.expl.var <- cumsum(expl.var)
  
  # op <-par(mfrow=c(2,2))
  # plot(eig.vals.V, main="Eigenvalues of matrix V",ylab="Eigenvalues", type="b")
  # for (j in (1:p)){
  #   plot(eig.V$vectors[,j],main=paste("Eigenvector",j,", Expl.Var.:",expl.var[j],"%"))
  #   abline(h=0,col=2,lty=2)
  # }
  # par(op)
  
  
  eig.V.df <- as.data.frame(eig.V$vectors)
  names(eig.V.df) <- colnames(A)
  
  op <-par(mfrow=c(nrows.plot,ncols.plot))
  plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
  
  if (!is.null(relev.ghost)){
    plot(relev.rp,relev.ghost,
         xlim=c(0,max(relev.rp)),
         ylim=c(0,max(relev.ghost)),
         xlab="Relev. by rand.perm.", 
         ylab="Relev. by ghost variables")
    pointLabel(relev.rp,relev.ghost, colnames(A))
  }else{
    if (!is.null(sum.lm.tr)){
      plot(F.transformed,relev.rp,
           xlim=c(0,max(c(F.transformed,relev.rp))),
           ylim=c(0,max(c(F.transformed,relev.rp))),
           xlab="2*beta*Var(x_i)", ylab="Relev. by rand.perm.")
      pointLabel(F.transformed,relev.rp, colnames(A))
      abline(a=0,b=1,col=2)
      #abline(v=F.critic.transformed,h=F.critic.transformed,lty=2,col="blue",lwd=2)
    }else{
      plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
    }
  }
  par(xaxp=c(1,p,min(p,5)))
  plot(eig.vals.V, ylim=c(0,max(eig.vals.V)),
       main=expression(paste("Eigenvalues of matrix ",tilde(V))),
       ylab="Eigenvalues",type="b")
  abline(h=0,col="red",lty=2)
  
  par(op)
  
  pushViewport(viewport(layout = grid.layout(nrows.plot, ncols.plot))) #package grid
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  print(plot.rel.Gh,vp = vplayout(1,1))
  jj<-0
  for (j in vars){
    jj<-jj+1
    print(
      ggplot(eig.V.df) +
#       geom_bar(aes(x=var.names, y=eig.V.df[,j]),
        geom_bar(aes(x=reorder(names(eig.V.df),X=length(names(eig.V.df)):1), 
                     y=eig.V.df[,j]),
                 stat="identity") +
        geom_hline(aes(yintercept=0),color="red",linetype=2,size=1) +
        ylim(min(eig.V.df[,j])-.5,max(eig.V.df[,j])+.5) +
        coord_flip() +
        ggtitle(paste0("Eig.vect.",j,", Expl.Var.: ",expl.var[j],"%")) +
        theme(axis.title=element_blank(),plot.title = element_text(size = 12)),
      vp = vplayout(2+(jj-1)%/%ncols.plot, 1+(jj-1)%%ncols.plot)
    )
  }
}