#' @test check.constr(x=matrix(c(.1,.2,.5,.8),ncol=1),function(x)(x-0.5))
#' @test check.constr(x=matrix(c(.1,.2,.5,.8),ncol=1),function(x)(x<=0.5))
#' @test check.constr(x=matrix(c(.1,.2,.5,.8),ncol=2),function(x)(sum(x)<=.7))
#' @test check.constr(x=list(.45,.4),function(x)(sum(x)<=.8))
#' @test check.constr(x=c(.45,.4),function(x)(sum(x)<=.8))
#' @test check.constr(x=matrix(c(.45,.4),ncol=2),function(x)(sum(x)<=.8))
#' @test check.constr(x=as.matrix(expand.grid(seq(0,1,,10),seq(0,1,,10))),constraint.fun=function(x){(sum((x))>1)})
#' @test check.constr(x=as.matrix(expand.grid(seq(0,1,,10),seq(0,1,,10))),constraint.fun=NULL)
check.constr <- function(x,constraint.fun) {
  is.checked = Vectorize(function(x) {ifelse(is.logical(x),x,x<=0)})
  
  if (is.null(constraint.fun)) return(rep(TRUE,nrow(x)))
  
  if (is.list(x)) x=unlist(x)
  if (!is.matrix(x)) x = t(as.matrix(x))
  
  if (is.list(constraint.fun)) {
    constr = is.checked( apply(X=x,MARGIN=1,FUN=constraint.fun[[1]]))  # is.checked( ifelse( !is.matrix(x),constraint.fun(x),apply(X=x,MARGIN=1,FUN=constraint.fun[[1]]) ))
    if (length(constraint.fun)>1) {
      for (i in 2:length(constraint.fun)) {
        constr = constr * is.checked(apply(X=x,MARGIN=1,FUN=constraint.fun[[i]])) #is.checked( ifelse( !is.matrix(x),constraint.fun(x),apply(X=x,MARGIN=1,FUN=constraint.fun[[i]]) ))
      }
    }
  } else {
    constr = is.checked(apply(X=x,MARGIN=1,FUN=constraint.fun)) #is.checked( ifelse( !is.matrix(x),constraint.fun(x),apply(X=x,MARGIN=1,FUN=constraint.fun) ))
  }
  return(constr)
}

#' @test X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); constrEI(seq(0,1,,100),kmi,constraint.fun=function(x)(x<0.49))
#' @test X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)constrEI(x,kmi,constraint.fun=function(x)(x<0.49)),dim=1)
#' @test X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); constrEI(as.matrix(expand.grid(seq(0,1,,10),seq(0,1,,10))),kmi,constraint.fun=function(x){(sum((x))>1)})
#' @test X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)constrEI(x,kmi,constraint.fun=function(x){(sum((x))>1)}),dim=2)
constrEI <- function(x, model, constraint.fun=NULL, plugin=NULL, type="UK", envir=NULL) {
  #cat("      constrEI\n")
  d=model@d
  
  if (!is.matrix(x)) x = matrix(x,ncol=d)
  #if (!is.null(print.constrEI)) print.constrEI(x)
  #     if (is.null(constraint.fun)) return(EI(x,model=model,type=type, envir=envir))
  # 
  #     if (is.list(constraint.fun)) {
  #     const = apply(X=x,MARGIN=1,FUN=constraint.fun[[1]])<0
  #     if (length(constraint.fun)>1) {
  #       for (i in 2:length(constraint.fun)) {
  #         const = const * (apply(X=x,MARGIN=1,FUN=constraint.fun[[i]])<0)
  #       }
  #     }
  #     }else{
  #       const = apply(X=x,MARGIN=1,FUN=constraint.fun)<0
  #     }
  #     return(EI(x,model=model,plugin=plugin,type=type, envir=envir) * const)
  
  constr<-check.constr(x,constraint.fun)
  if(length(which(constr))>0)
    ei=EI(x[which(constr),],model=model,plugin=plugin,type=type, envir=envir)
  else 
    ei=0
  ei.constr=rep(0,nrow(x))
  ei.constr[which(constr)]=ei
  
  #cat("EI");print(ei)
  #cat("constr");print(constr)
  #cat("EI*constr");print(ei*constr)
  #cat("      /constrEI\n")
  return(ei.constr)
}

#   constrEI.grad = function(x) {
#     if (!is.matrix(x)) x = matrix(x,ncol=d)
#     if (length(constraint.fun)>1) {
#       const = (constraint.fun[1](x)<0)
#       for (c in 2:length(constraint.fun)) {
#         const = const * (constraint.fun[i](x)<0)
#       }
#     } else 
#       const = const * (constraint.fun(x)<0)
#     return(EI.grad(x,model=model,envir=EI.envir) * const)
#   }


#' @test X=matrix(runif(10),ncol=1); y=-sin(2*pi*X); kmi <- km(design=X,response=y); EI(seq(0,1,,100),model=kmi)
#' @test X=matrix(runif(10),ncol=1); y=-sin(2*pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)EI(x,kmi),dim=1)
#' @test X=matrix(runif(10),ncol=1); y=-sin(2*pi*X); kmi <- km(design=X,response=y);  sectionview.km(kmi,xlim=c(0,1)); DiceView::sectionview.fun(function(x)EI(x,kmi,plugin=function(x)localmin.fun(x,kmi@X,kmi@y)),dim=1,xlim=c(0,1)); DiceView::sectionview.fun(function(x)EI(x,kmi),col='red',dim=1,add=T);
#' @test X=matrix(runif(20),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)EI(x,kmi),dim=2); points(X)
#' @test set.seed(1); X=matrix(runif(100),ncol=2); y=branin(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)EI(x,kmi,plugin=function(x)localmin.fun(x,kmi@X,kmi@y)),dim=2); DiceView::contourview.fun(function(x)EI(x,kmi,plugin=NULL),dim=2,col='red',add=T); points(X); points(X[which.localmin.delaunay(X,y),],col='red',pch=20);
EI <- function (x, model, plugin=NULL, type="UK", envir=NULL) {
  
  if (is.null(plugin)){ plugin <- min(model@y) }
  if (!is.function(plugin)) m=function(x) rep(plugin,nrow(x)) else  m <- plugin
  
  ########################################################################################
  # Convert x in proper format(s)
  if (!is.matrix(x)) x <- matrix(x,ncol= model@d)
  d <- ncol(x)
  if (d != model@d){ stop("x does not have the right number of columns (",d," instead of ",model@d,")") }
  newdata <- x
  colnames(newdata) = colnames(model@X) 
  
  ########################################################################################
  #cat("predict...")
  if (!is.numeric(newdata)) 
    newdata = as.numeric(newdata)
  if(is.matrix(newdata)) { 
    if (ncol(newdata)!=d) newdata = as_x(newdata)
  } else {
    #if (length(newdata)!=d) 
    newdata = as_x(newdata)
  }
  predx <- predict.km(object=model, newdata=newdata, type=type, checkNames = FALSE)
  #cat(" done.\n")
  kriging.mean <- predx$mean
  kriging.sd   <- predx$sd
  
  xcr <- (m(x) - kriging.mean)/kriging.sd
  
  #if (kriging.sd/sqrt(model@covariance@sd2) < 1e-06) 
  #{ res <- 0
  #  xcr <- xcr.prob <- xcr.dens <- NULL
  #} else  {
  xcr.prob <- pnorm(xcr)
  xcr.dens <- dnorm(xcr)	        
  res <- (m(x) - kriging.mean) * xcr.prob + kriging.sd * xcr.dens
  #}
  too.close = which(kriging.sd/sqrt(model@covariance@sd2) < 1e-06)
  res[too.close] <- pmax(0,m(x)[too.close] - kriging.mean[too.close])
  
  if (!is.null(envir)) 
  { 
    xcr[too.close] <- NaN
    xcr.prob[too.close] <- NaN
    xcr.dens[too.close] <- NaN
    assign("xcr", xcr, envir=envir)
    assign("xcr.prob", xcr.prob, envir=envir)
    assign("xcr.dens", xcr.dens, envir=envir)
    assign("kriging.sd", kriging.sd, envir=envir)
    assign("c", predx$c, envir=envir)
    assign("Tinv.c", predx$Tinv.c, envir=envir)
    
  }
  return(res)
}




EI.grad <- function(x, model, plugin=NULL, type="UK", envir=NULL){ 
  
  ########################################################################################
  if (is.null(plugin)){ plugin <- min(model@y) }
  m <- plugin
  
  ########################################################################################
  # Convert x in proper format(s)
  d <- length(x)
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  ########################################################################################
  # Get quantities related to the model
  T <- model@T
  X <- model@X
  z <- model@z
  u <- model@M
  covStruct <- model@covariance
  
  # Get quantities related to the prediction
  if (is.null(envir))
  {  
    predx <- predict(object=model, newdata=newdata, type=type, checkNames = FALSE,se.compute=TRUE,cov.compute=FALSE)
    kriging.mean <- predx$mean
    kriging.sd <- predx$sd
    v <- predx$Tinv.c
    c <- predx$c
    
    xcr <- (m - kriging.mean)/kriging.sd
    xcr.prob <- pnorm(xcr)
    xcr.dens <- dnorm(xcr)    
  } else
  {  # If uploaded through "envir", no prediction computation is necessary 
    toget <- matrix(c("xcr", "xcr.prob", "xcr.dens", "kriging.sd", "c", "Tinv.c"), 1, 6)
    apply(toget, 2, get, envir=envir)
    xcr        <- envir$xcr
    xcr.prob   <- envir$xcr.prob
    xcr.dens   <- envir$xcr.dens
    kriging.sd <- envir$kriging.sd
    c          <- envir$c
    v          <- envir$Tinv.c
  }
  
  F.newdata <- model.matrix(model@trend.formula, data=newdata)
  
  ########################################################################################
  # Pursue calculation only if standard deviation is non-zero
  if ( kriging.sd/sqrt(model@covariance@sd2) < 1e-06) 
  { ei.grad <- rep(0,d)
  } else 
  { # Compute derivatives of the covariance and trend functions
    dc <- covVector.dx(x=newdata.num, X=X, object=covStruct, c=c)  
    f.deltax <- trend.deltax(x=newdata.num, model=model)
    
    # Compute gradients of the kriging mean and variance
    W <- backsolve(t(T), dc, upper.tri=FALSE)
    kriging.mean.grad <- t(W)%*%z + t(model@trend.coef%*%f.deltax)
    
    if (type=="UK")
    { tuuinv <- solve(t(u)%*%u)
      kriging.sd2.grad <-  t( -2*t(v)%*%W +
                                2*(F.newdata - t(v)%*%u )%*% tuuinv %*%
                                (f.deltax - t(t(W)%*%u) ))
    } else
    { kriging.sd2.grad <-  t( -2*t(v)%*%W) }
    
    kriging.sd.grad <- kriging.sd2.grad / (2*kriging.sd)
    
    # Compute gradient of EI
    ei.grad <- - kriging.mean.grad * xcr.prob + kriging.sd.grad * xcr.dens
  }
  ########################################################################################  
  return(ei.grad)
}