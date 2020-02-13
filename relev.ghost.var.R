### Computing the case-variable relevance case, 
### using the Ghost Variables strategy

relev.ghost.var <- function(model,newdata=model$call$data,
                           func.model.ghost.var=gam, ...){
  #func.model <- eval(parse(text=class(model)[1])) # What kind of model has been fitted 
  #data.tr <- model$call$data[model$call$subset,] # data used for training the model
  #n <- dim(data.tr)[1]
  attr.model <- attributes(model$terms) #getting the varaible names in the model
  term.labels <- attr.model$term.labels #getting the varaible names in the model
  if (formals(func.model.ghost.var)[[2]]==formals(gam)[[2]]){
    s.term.labels <- paste0("s(",term.labels,")")
  }else{
    s.term.labels <- term.labels
  }
  p <- length(term.labels)

  n2 <- dim(newdata)[1] # newdata is the test sample
  # Predicting in the test sample
  y.hat.ts <- as.numeric( predict(model,newdata = newdata) )
  
  A <- matrix(0,nrow=n2, ncol=p)
  colnames(A) <- term.labels
  GhostX <- A
  yX.ts <- newdata
  for (j in (1:p)){
    yX.ts.aux <- yX.ts
    form.j <- as.formula(paste0(term.labels[j],"~",paste(s.term.labels[-j],collapse="+")))
    xj.hat <- func.model.ghost.var(form.j, data=yX.ts[,term.labels],...)$fitted.values
    yX.ts.aux[,term.labels[j]] <- xj.hat
    y.hat.ts.j <- as.numeric( predict(model,newdata = yX.ts.aux) )
    A[,j] <- y.hat.ts - y.hat.ts.j
    GhostX[,j] <- xj.hat
  }
  V=(1/n2)*t(A)%*%A
  relev.ghost <- diag(V)
  eig.V <- eigen(V)
  return(list(A=A, V=V, GhostX=GhostX, 
              relev.ghost=relev.ghost, 
              eig.V=eig.V,
              y.hat.test=y.hat.ts)
         )
}

plot.relev.ghost.var <- function(relev.ghost.out, n1, resid.var,
                                 vars=NULL, sum.lm.tr=NULL,
                                 alpha=.01, ncols.plot=3){
  A <- relev.ghost.out$A
  V <- relev.ghost.out$V
  eig.V <- relev.ghost.out$eig.V
  GhostX <- relev.ghost.out$GhostX
  relev.ghost <- relev.ghost.out$relev.ghost

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
    F.transformed <- resid.var*sum.lm.tr$coefficients[-1,3]^2/n1
  }
  F.critic.transformed <- resid.var*qf(1-alpha,1,n1-p-1)/n1

  rel.Gh <- data.frame(relev.ghost=relev.ghost)
  rel.Gh$var.names <- colnames(A)
  
  plot.rel.Gh <- ggplot(rel.Gh) +
    geom_bar(aes(x=reorder(var.names,X=length(var.names):1), y=relev.ghost), 
             stat="identity", fill="darkgray") +
    ggtitle("Relev. by ghost variables") +
    geom_hline(aes(yintercept = F.critic.transformed),color="blue",size=1.5,linetype=2)+
    theme(axis.title=element_blank())+
    theme_bw()+
    ylab("Relevance")+
    xlab("Variable name") +
    coord_flip()
    
  plot.rel.Gh.pctg <- ggplot(rel.Gh) +
    geom_bar(aes(x=reorder(var.names,X=length(var.names):1), 
                 y=100*relev.ghost/sum(relev.ghost)), 
             stat="identity", fill="darkgray") +
    coord_flip() +
    ggtitle("Relev. by ghost variables (% of total relevance)") +
    theme(axis.title=element_blank())
  
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
  eig.V.df$var.names <- colnames(A)
  
  op <-par(mfrow=c(nrows.plot,ncols.plot))
  plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
  
  if (!is.null(sum.lm.tr)){
    plot(F.transformed,relev.ghost,
         xlim=c(0,max(c(F.transformed,relev.ghost))),
         ylim=c(0,max(c(F.transformed,relev.ghost))),
         xlab=expression(paste("F-statistics*",hat(sigma)^2/n[1])), 
         ylab="Relev. by ghost variables")
    pointLabel(F.transformed,relev.ghost, colnames(A))
    abline(a=0,b=1,col=2)
    abline(v=F.critic.transformed,h=F.critic.transformed,lty=2,col="blue",lwd=2)
  }else{
    plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
  } 
  
  par(xaxp=c(1,p,min(p,5)))
  plot(eig.vals.V, ylim=c(0,max(eig.vals.V)),
       main=expression("Eigenvalues of matrix V"),
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
        geom_bar(aes(x=reorder(eig.V.df$var.names,X=length(eig.V.df$var.names):1), 
                     y=eig.V.df[,j]), stat="identity") +
        geom_hline(aes(yintercept=0),color="red",linetype=2,size=1) +
        ylim(min(eig.V.df[,j])-.5,max(eig.V.df[,j])+.5) +
        coord_flip() +
        ggtitle(paste0("Eig.vect.",j,", Expl.Var.: ",expl.var[j],"%")) +
        theme(axis.title=element_blank(),plot.title = element_text(size = 12)),
      vp = vplayout(2+(jj-1)%/%ncols.plot, 1+(jj-1)%%ncols.plot)
    )
  }
}

relevance2cluster <- function(V,method="ward.D2",ret.dist=FALSE,...){
  V2<-V^2
  dV2 <- diag(V2)
  p<- length(dV2)
  one <- 1+numeric(p)
  # Distance matrix from similarity matrix V2:
  W <- as.dist(sqrt(one %*% t(dV2) + dV2 %*% t(one) -2*V2))
  hcl.W <- hclust(W, method = method)
  plot(hcl.W,...)
  if (ret.dist) return(W)
}