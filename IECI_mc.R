#' @test X=matrix(runif(100),50,2);any.row.equals(X,X[1,])
#' @test X=matrix(runif(100),50,2);any.row.equals(X,X[1,]+1)
any.row.equals <- function(X,x) {
  #X<<-X
  #x<<-x
  if (is.null(x)) return(FALSE)
#   for (i in 1:nrow(X)) {
#     eq=TRUE
#     for (j in 1:ncol(X)) {
#       if (X[i,j]!=x[j]) {eq=FALSE; break;}
#     }
#     if (eq)    {  	
#       return(TRUE)
#     }
#   }
#   return(FALSE)
return(any(rowSums(abs(matrix(as.numeric(X),nrow(X),ncol(X)) - matrix(as.numeric(x),nrow(X),ncol(X),byrow=TRUE)))==0))
}

# set.seed(1);x=matrix(runif(1000000*3),1000000);system.time(any(rowSums(subset(x==x[1,]))==3))
# utilisateur     système      écoulé 
# 0.082       0.000       0.083

# set.seed(1);x=matrix(runif(1000000*3),1000000);system.time(any.row.equals(x,x[1,]))
# utilisateur     système      écoulé 
# 0.228       0.000       0.228

# rowMatch <- function(A,B) {
#   # Rows in A that match the rows in B
#   # The row indexes correspond to A
#   f <- function(...) paste(..., sep=":")
#   #if(!is.matrix(B)) B <- matrix(B, 1, length(B))
#   a <- do.call("f", as.data.frame(A))
#   b <- do.call("f", as.data.frame(B)) 
#   match(b, a)
# } 
# set.seed(1);x=matrix(runif(1000000*3),1000000);system.time(!is.na(rowMatch(x,x[1,])))

library(randtoolbox)

#' Compute "Expected Conditional Improvement" (ECI_mc) criterion
#' 
#' ECI_mc(x0,x) = E_y0~Y(x0)[ EI(x) |  Y(x0)=y0 ]
#'           = E_y0~Y(x0)[ E_y~Y(x)[ ( min(Y(X0)) - y )+ | Y(X0)=Y0 ] | Y(x0)=y0 ]
#' 
#' @param x0 point where a virtual experiment is added, with a response following Y process
#' @param x  point where the EI is evaluated
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)ECI_mc(x,x0=0.5,model=kmi),dim=1,npoints=101); abline(v=X); abline(v=0.5,col='red')
ECI_mc <- function(x, x0, model, y0.sample,...) {
  return(ECI_mc_vec.o2(x, x0, model, y0.sample=100, ...))
}

#' unoptimized version
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)ECI_mc.o0(x,x0=0.5,model=kmi),dim=1,npoints=20); abline(v=X); abline(v=0.5,col='red')
ECI_mc.o0 <- function(x, x0, model,  y0.sample=10, mc=1) {
  x  = t(as.matrix(x))
  x0 = t(as.matrix(x0))
  
  Yx0 = predict_nobias_km( object = model, newdata = x0,  cov.compute = FALSE,  se.compute = TRUE) 
  my0 = Yx0$mean
  sdy0 = Yx0$sd
  Fy0 = Yx0$F.newdata
  
  EIx_y0 = function(y0) {
    if (length(model@noise.var) != 0) {
      new.noise.var <- c(model@noise.var,0)
    } else {
      new.noise.var <- rep(0,model@n+1)            
    }
    k_y0 = km(formula=model@trend.formula, design=rbind(model@X, as.matrix(x0)),response=rbind(model@y, as.matrix(y0)),
              covtype=model@covariance@name,coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance),coef.var=model@covariance@sd2, 
              noise.var=new.noise.var)
    
    Yx.update = predict_nobias_km(object = k_y0, newdata = x,  cov.compute = FALSE, se.compute = TRUE)
    
    minYn <- min(y0,min(k_y0@y))      
    yimp <- (minYn - Yx.update$mean) / Yx.update$sd
    p.yimp <- pnorm(yimp)
    d.yimp <- dnorm(yimp)          
    EI <- (minYn - Yx.update$mean) * p.yimp + Yx.update$sd * d.yimp
    
    #         points(x=x0,y=-y0)
    #         X=seq(from=0,to=1,l=100)
    #         Y=predict.km(k_y0,newdata=X)$mean
    #         lines(X,-Y,col=rgb(0,0,0,.1))
    
    return(EI)
  }
  
  if (length(y0.sample)==1) 
    y0.sample = my0+sdy0*qnorm(sobol(n=y0.sample,init=T))#rnorm(n=y0.sample, my0, sdy0)
  else 
    y0.sample = y0.sample*sdy0 + my0
  
  if(mc<=1) {
    sampleEI_Y0 = unlist(lapply( FUN=EIx_y0, X=y0.sample ))        
  } else {
    library(multicore)
    sampleEI_Y0 = unlist(mclapply( FUN=EIx_y0, X=y0.sample, mc.cores=as.integer(mc) ))
  }
  
  return(mean(sampleEI_Y0))
}

library(KrigInv)

#' level 1 optimized version : use update.km to build updated kriging
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)ECI_mc.o1(x,x0=0.5,model=kmi),dim=1,npoints=20); abline(v=X); abline(v=0.5,col='red')
ECI_mc.o1 <- function(x, x0, model,  y0.sample=10, mc=1) {
  
  x  = t(as.matrix(x))
  x0 = t(as.matrix(x0))
  
  Yx0 = predict_nobias_km(object = model, newdata = x0,  cov.compute = FALSE, se.compute = TRUE) 
  my0 = Yx0$mean
  sdy0 = Yx0$sd
  Fy0 = Yx0$F.newdata
  
  EIx_y0 = function(y0) {
    k_y0 = update_km(model = model, NewX = x0, NewY =  y0, F.newdata = Fy0)
    
    Yx.update = predict_nobias_km(object = k_y0, newdata = x,  cov.compute = FALSE, se.compute = TRUE)
    
    minYn <- min(y0,min(k_y0@y))      
    yimp <- (minYn - Yx.update$mean) / Yx.update$sd
    p.yimp <- pnorm(yimp)
    d.yimp <- dnorm(yimp)          
    EI <- (minYn - Yx.update$mean) * p.yimp + Yx.update$sd * d.yimp
    
    #         points(x=x0,y=-y0)
    #             X=seq(from=0,to=1,l=100)
    #             Y=predict.km(k_y0,newdata=X)$mean
    #             lines(X,-Y,col=rgb(0,0,0,.1))
    
    return(EI)
  }
  
  if (length(y0.sample)==1) 
    y0.sample = my0+sdy0*qnorm(sobol(n=y0.sample,init=T))#rnorm(n=y0.sample, my0, sdy0)
  else 
    y0.sample = y0.sample*sdy0 + my0
  
  if(mc<=1) {
    sampleEI_Y0 = unlist(lapply( FUN=EIx_y0, X=y0.sample ))        
  } else {
    library(multicore)
    sampleEI_Y0 = unlist(mclapply( FUN=EIx_y0, X=y0.sample, mc.cores=as.integer(mc) ))
  }
  
  return(mean(sampleEI_Y0))
}

library(KrigInv)

#' level 2 optimized version : use Emery's formulas to update kriging
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)ECI_mc.o2(x,x0=0.5,model=kmi),dim=1,npoints=20); abline(v=X); abline(v=0.5,col='red')
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); ECI_mc.o2(x=matrix(runif(2),ncol=2),x0=matrix(.5,ncol=2),model=kmi)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); apply(FUN=ECI_mc.o2,X=matrix(runif(20),ncol=2),MARGIN=1,x0=matrix(.5,ncol=2),model=kmi)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)ECI_mc.o2(x,x0=matrix(.5,ncol=2),model=kmi),dim=2,npoints=20); points(X); points(.5,.5,col='red'); DiceView::contourview.fun(function(x)EI(x=x,model=kmi),dim=2,npoints=20,add=TRUE,col='red')
ECI_mc.o2 <- function(x, x0, model,  y0.sample=10, mc=1, kn=NULL, precalc.data=NULL) {
  
  #if(any.row.equals(model@X,x)) return(0)
  
  if (!is.matrix(x)) x = matrix(x,nrow=1)
  
  Yx = predict_nobias_km(object = model, newdata = x,  cov.compute = FALSE, se.compute = TRUE,low.memory=TRUE) 
  #Yx = predict.km(bias.correct=TRUE,object = model, newdata = x,  cov.compute = FALSE, se.compute = TRUE,low.memory=TRUE) 
  my = Yx$mean
  sdy = Yx$sd
  
  if (!is.matrix(x0)) x0 = t(as.matrix(x0))    
  Yx0 = predict_nobias_km(object = model, newdata = x0,  cov.compute = FALSE, se.compute = TRUE) 
  my0 = Yx0$mean
  sdy0 = Yx0$sd
  Fy0 = Yx0$F.newdata
  cy0 = Yx0$c
  
  if (is.null(kn)) {
    if (is.null(precalc.data)) precalc.data <- precomputeUpdateData(model,x)
    kn <- computeQuickKrigcov(model,x,x0,precalc.data, Fy0 , cy0)
  }
  
  EIx_y0 = function(y0) {
    
    Yx.update = predict_update_km(
      newXvalue = y0, newXmean = my0, newXvar = sdy0^2, newdata.oldmean = my,  newdata.oldsd = sdy, 
      kn = kn) 
    
    minYn <- min(y0,min(model@y))      
    yimp <- (minYn - Yx.update$mean) / Yx.update$sd
    p.yimp <- pnorm(yimp)
    d.yimp <- dnorm(yimp)          
    EI <- (minYn - Yx.update$mean) * p.yimp + Yx.update$sd * d.yimp
    
    return(EI)
  }
  
  if (length(y0.sample)==1) 
    y0.sample = my0+sdy0*qnorm(sobol(n=y0.sample,init=T))#rnorm(n=y0.sample, my0, sdy0)
  else 
    y0.sample = y0.sample*sdy0 + my0
  
  if(mc<=1) {
    sampleEI_Y0 = unlist(lapply( FUN=EIx_y0, X=y0.sample ))        
  } else {
    library(multicore)
    sampleEI_Y0 = unlist(mclapply( FUN=EIx_y0, X=y0.sample, mc.cores=as.integer(mc) ))
  }
  
  return(mean(sampleEI_Y0))
}

#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); ECI_mc_vec.o2(.51,x0=0.5,model=kmi);ECI_mc_vec.o2(.49,x0=0.5,model=kmi);ECI_mc_vec.o2(c(.51,.49),x0=0.5,model=kmi);
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); plot(function(x)ECI_mc_vec.o2(x,x0=0.5,model=kmi)); abline(v=X); abline(v=0.5,col='red')
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); x=matrix(runif(2),ncol=2); ECI_mc.o2(x,x0=matrix(.5,ncol=2),model=kmi); ECI_mc_vec.o2(x,x0=matrix(.5,ncol=2),model=kmi)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); x=matrix(runif(20),ncol=2); apply(FUN=ECI_mc.o2,X=x,MARGIN=1,x0=matrix(.5,ncol=2),model=kmi); ECI_mc_vec.o2(x,x0=matrix(.5,ncol=2),model=kmi)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)ECI_mc.o2(x,x0=matrix(.5,ncol=2),model=kmi),dim=2,npoints=20); points(X); points(.5,.5,col='red'); DiceView::contourview.fun(function(x)EI(x=x,model=kmi),dim=2,npoints=20,add=TRUE,col='red')
ECI_mc_vec.o2 <- function(x, x0, model,  y0.sample=100, kn=NULL, precalc.data=NULL) {
  #cat("        ECI_mc_vec.o2\n")
  
  if (!is.matrix(x)) x = as.matrix(x)
  Yx = predict_nobias_km(object = model, newdata = x,  cov.compute = FALSE, se.compute = TRUE,low.memory=TRUE) 
  #Yx = predict.km(bias.correct=TRUE,object = model, newdata = x,  cov.compute = FALSE, se.compute = TRUE,low.memory=TRUE) 
  my = Yx$mean
  sdy = Yx$sd
  
  if (!is.matrix(x0)) x0 = t(as.matrix(x0))    
  Yx0 = predict_nobias_km(object = model, newdata = x0,  cov.compute = FALSE, se.compute = TRUE) 
  my0 = Yx0$mean
  sdy0 = Yx0$sd
  Fy0 = Yx0$F.newdata
  cy0 = Yx0$c
  
  #print(my0)
  
  if (is.null(kn)) {
    if (is.null(precalc.data)) precalc.data <- precomputeUpdateData(model,x)
    kn <- computeQuickKrigcov(model,x,x0,precalc.data, Fy0 , cy0)
  }
  
  EIx_y0 = function(y0) {
    
    Yx.update = predict_update_km(
      newXvalue = y0, newXmean = my0, newXvar = sdy0^2, newdata.oldmean = my,  newdata.oldsd = sdy, 
      kn = kn) 
    
    minYn <- min(y0,min(model@y))      
    yimp <- (minYn - Yx.update$mean) / Yx.update$sd
    p.yimp <- pnorm(yimp)
    d.yimp <- dnorm(yimp)          
    EI <- (minYn - Yx.update$mean) * p.yimp + Yx.update$sd * d.yimp
    
    return(EI)
  }
  
  if (length(y0.sample)==1) 
    y0.sample = my0+sdy0*qnorm(sobol(n=y0.sample,init=T))#rnorm(n=y0.sample, my0, sdy0)
  else 
    y0.sample = y0.sample*sdy0 + my0
  
  sampleEI_Y0 = plyr::maply( .fun=EIx_y0, .data=y0.sample )
  
  #anyrowequals_X <<- apply(FUN=function(x)any.row.equals(model@X,x),X=(as.matrix(x)),MARGIN=1)
  #print(anyrowequals_X)
  #anyrowequals_x0 <<- apply(FUN=function(xx)any.row.equals(as.matrix(x0),xx),X=(as.matrix(x)),MARGIN=1)
  #print(anyrowequals_x0)
  #if (length(x0)>1) {
  #  if (length(which(anyrowequals_x0))>0) sampleEI_Y0[,which(anyrowequals_x0)] = 0
  #  if (length(which(anyrowequals_X))>0) sampleEI_Y0[,which(anyrowequals_X)] = 0
  #}
  
  #cat("        /ECI_mc_vec.o2\n")
  if (nrow(x)==1) return(mean(sampleEI_Y0)) else return(colMeans(sampleEI_Y0))
}

#' @param x0
#' @param model 
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); IECI_mc(x0=.5,model=kmi)
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)IECI_mc(x0=x,model=kmi,lower=c(0,0),upper=c(1,1)),dim=1,npoints=100); xx=seq(f=0,t=1,l=100); lines(xx,EI(xx,model=kmi),col='red'); abline(v=X)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); IECI_mc(x0=c(.5,.5),model=kmi,lower=c(0,0),upper=c(1,1))
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)IECI_mc(x0=x,model=kmi,lower=c(0,0),upper=c(1,1)),dim=2,npoints=20); DiceView::contourview.fun(function(x)EI(x=x,model=kmi),dim=2,npoints=20,add=TRUE,col='red'); points(X);
#' @test set.seed(1);X=matrix(runif(200),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)IECI_mc(x0=x,model=kmi,lower=c(0,0),upper=c(1,1)),dim=2,npoints=30); DiceView::contourview.fun(function(x)EI(x=x,model=kmi),dim=2,npoints=30,add=TRUE,col='red'); points(X);
IECI_mc <- function(x0, model, lower=0,upper=1,y0.sample=10, ...) {
  return(IECI_mc.o1(x0, model, lower,upper,y0.sample=y0.sample, ...))
}

# #' @test sobol_inside(10,2,0,1)
# #' @test sobol_inside(10,1,0,1)
# sobol_inside <- function(n,dim,lower,upper) {
#   s = sobol(n=n,dim=dim)
#   #s * (upper-lower) - lower
#   matrix(s * (upper-lower) - lower,ncol=dim)
# }

#' @test set.seed(1);X=matrix(runif(200),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)IECI_mc.o1(x0=x,model=kmi,upper=c(1,1),lower=c(0,0)),dim=2,npoints=30); DiceView::contourview.fun(function(x)EI(x=x,model=kmi),dim=2,npoints=30,add=TRUE,col='red'); points(X);
#' @test set.seed(1);X=matrix(runif(200),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); mEI=max_EI.xstart(kmi,lower=c(0,0),upper=c(1,1))$value; sample.EI = sample.resampling(density.fun = function(x) mEI/1000+EI(x,model),lower,upper,N=10^d); DiceView::contourview.fun(function(x)IECI_mc.o1(x0=x,model=kmi,upper=c(1,1),lower=c(0,0),integration.points = sample.EI$sample,integration.weights = 1/sample.EI$sample.density),dim=2,npoints=30); DiceView::contourview.fun(function(x)EI(x=x,model=kmi),dim=2,npoints=30,add=TRUE,col='red'); points(X);
IECI_mc.o1 <- function (x0, model, lower,upper,y0.sample, plot.integration=NULL,...) {
  if (any.row.equals(model@X,x0)) {
    return(NaN)
  }
  
  opts=list(...)  
  integration.points=opts[['integration.points']]
  integration.weights=opts[['integration.weights']]
  if (is.null(integration.weights)) integration.weights=1
  
  if (!is.matrix(x0)) x0=matrix(x0,ncol=model@d)
  
  d=ncol(x0)
  if (is.null(integration.points) || !is.matrix(integration.points)){
    if (is.null(integration.points)) {
      n=10^(d+1)
    } else if (length(integration.points)==1) {
      n=integration.points
    } else n=nrow(integration.points)
    integration.points = sobol_inside(n,d,lower,upper)
  }
  
  if(length(integration.weights)==1) 
    integration.weights = rep(integration.weights,nrow(integration.points))
  if (sum(integration.weights) != 1) 
    integration.weights = integration.weights / sum(integration.weights)
  
  if (!is.null(plot.integration)) try(plot.integration(integration.points,integration.weights))
  
  precalc.data = opts[['precalc.data']]
  new_IEI = sum(ECI_mc(x=integration.points,x0=x0,model=model,precalc.data=precalc.data,y0.sample=y0.sample) * integration.weights)
  return(new_IEI)
}


#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); constrIECI_mc(x0=.5,model=kmi,constraint.fun=function(x)(x<0.49),integration.points=100)
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)constrIECI_mc(x0=x,model=kmi,constraint.fun=function(x)(x<0.5),integration.points=100),dim=1,npoints=100); DiceView::sectionview.fun(function(x)IECI_mc(x0=x,model=kmi,integration.points=100),dim=1,npoints=100,add=TRUE,col='green'); xx=seq(f=0,t=1,l=100); lines(xx,EI(xx,model=kmi),col='red'); abline(v=X)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); IECI_mc(x0=c(.5,.5),model=kmi); constrIECI_mc(x0=c(.5,.5),model=kmi,constraint.fun=NULL)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); IECI_mc(x0=c(.5,.5),model=kmi); constrIECI_mc(x0=c(.5,.5),model=kmi,constraint.fun=function(x)sum(x)<1)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)constrIECI_mc(x0=x,model=kmi,constraint.fun=function(x)sum(x)>1),dim=2,npoints=20); DiceView::contourview.fun(function(x)IECI_mc(x0=x,model=kmi),dim=2,npoints=20,add=TRUE,col='green'); DiceView::contourview.fun(function(x)EI(x=x,model=kmi),dim=2,npoints=20,add=TRUE,col='red'); points(X);
constrIECI_mc <- function(x0, model, constraint.fun, lower, upper, ...) {
  return(constrIECI_mc.o1(x0, model, constraint.fun, lower, upper, ...))
}

constrIECI_mc.o1 <- function (x0, model, constraint.fun, lower, upper, ...) {
  #cat("      constrIECI_mc.o1\n")
  if (any.row.equals(model@X,x0)) {
    return(NaN)
  }
  
  opts=list(...)  
  integration.points=opts[['integration.points']]
  integration.weights=opts[['integration.weights']]
  if (is.null(integration.points) & is.null(integration.weights)) integration.weights=1
  
  if (is.null(constraint.fun)) return(IECI_mc.o1(x0,model, lower, upper, integration.points=integration.points, integration.weights=integration.weights, ...))
  
  if (!is.matrix(x0)) x0=matrix(x0,ncol=model@d)
  
  d=length(x0)
  if (is.null(integration.points) || !is.matrix(integration.points)){
    if (is.null(integration.points)) {
      n=10^d
    } else if (length(integration.points)==1) {
      n=integration.points
    } else n=nrow(integration.points)
    integration.points = sobol_inside(n,d,lower,upper)
  }
  
  constr_not_checked = which(!check.constr(integration.points,constraint.fun))
  integration.points = integration.points[-constr_not_checked,]
  
  if(length(integration.weights)==1) 
    integration.weights = rep(integration.weights,nrow(integration.points))
  else 
    integration.weights = integration.weights[-constr_not_checked]
  if (sum(integration.weights) != 1) 
    integration.weights = integration.weights / sum(integration.weights)
  
  #print(integration.points )
  #print(integration.weights)

  precalc.data = opts[['precalc.data']]
  new_IEI = sum(ECI_mc(x=integration.points,x0=x0,model=model,precalc.data=precalc.data) * integration.weights)
  #cat("      /constrIECI_mc.o1\n")
  return(new_IEI)
}
