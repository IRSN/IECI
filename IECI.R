#' "Integerated Expected Conditional Improvement" criterion. It is the reduction over integral EI, obtained when x is added in the doe.
#' 
#'            /             /   /
#' IECI(x0) = | ECI(x,x0) = |   | EI(x | Y(x0)=y0)
#'            / x           / x / y0~N()
#'      


#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); IEI(model=kmi)
IEI <- function(model,  integration.points=NULL, integration.weights=1, lower=0, upper=1, ...) {
  d=length(ncol(model@X))
  if (is.null(integration.points) || !is.matrix(integration.points)){
    if (is.null(integration.points)) {
      n=10^d
    } else if (length(integration.points)==1) {
      n=integration.points
    } else n=nrow(integration.points)
    integration.points = sobol_inside(n,d,lower,upper)
  }  
  
  if(length(integration.weights)==1) 
    integration.weights = rep(integration.weights,nrow(integration.points))
  if (sum(integration.weights) != 1) 
    integration.weights = integration.weights / sum(integration.weights)
  
  IEI = sum(apply(FUN=EI,X=integration.points,MARGIN=1,model=model,...) * integration.weights)      
  return(IEI)
}

source("int_Phi.R")
require(KrigInv)

#' @test X=matrix(runif(100),50,2);row.equals(X,X[1,])
#' @test X=matrix(runif(100),50,2);row.equals(X,X[1:2,])
#' @test X=matrix(runif(100),50,2);row.equals(X,x=rbind(X[1,],X[3,]))
#' @test X=matrix(runif(100),50,2);row.equals(X,X[1,]+1)
.fastcheck <- function(x,matrix){
  nc <- ncol(matrix)
  rec.check <- function(r,i,id){
    id[id] <- matrix[id,i] %in% r[i]
    if(i<nc & any(id)) rec.check(r,i+1,id) else any(id)
  }
  base::apply(x,1,rec.check,1,rep(TRUE,nrow(matrix)))
}
row.equals <- function(X,x) {
  if (is.null(x)) return(FALSE)
  if (!is.matrix(x)) x= matrix(x,nrow=1)
  #return(which(rowSums(abs(matrix(as.numeric(X),nrow(X),ncol(X)) - matrix(as.numeric(x),nrow(X),ncol(X),byrow=TRUE)))==0))
  which(.fastcheck(X,x))
}


#' Compute "Expected Conditional Improvement" (ECI) criterion
#' 
#' ECI_x0[x] = E_y0~Y(x0)[ EI[x] |  Y(x0)=y0 ]
#'           = E_y0~Y(x0)[ E_y~Y(x)[ ( min(Y(X)) - y )^+ | Y(x)=y ] | Y(x0)=y0 ]
#' 
#' @param x0 point where a virtual experiment is added, with a response following Y process
#' @param x  point where the EI is evaluated
#' 
#' Compare with Monte Carlo
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); ECI(.505,.5,kmi)
#' @test source("IECI_mc.R");set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); ECI_mc(.505,.5,kmi)
#' 
#' Vectorization
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); ECI(x=.51,x0=.5,model=kmi)
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); ECI(x=matrix(c(.51,.49)),x0=.5,model=kmi); ECI(x=.51,x0=.5,model=kmi); ECI(x=.49,x0=.5,model=kmi)
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); ECI(x=.51,x0=matrix(c(.3,.4,.5)),model=kmi);  ECI(x=.51,x0=.3,model=kmi);  ECI(x=.51,x0=.4,model=kmi);  ECI(x=.51,x0=.5,model=kmi)
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); ECI(x=matrix(c(.51,.49)),x0=matrix(c(.3,.4,.5)),model=kmi); ECI(x=matrix(c(.51,.49)),x0=.3,model=kmi); ECI(x=matrix(c(.51,.49)),x0=.4,model=kmi); ECI(x=matrix(c(.51,.49)),x0=.5,model=kmi)
#'
#' Compare with EI
#' @test x0=.5; set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)EI(x,model=kmi),dim=1,npoints=501,xlim=c(.4,.6)); abline(v=X); abline(v=x0,col='red');
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); x0=.54;plot(function(x)    EI(x,model=kmi),col='blue',type='l',                                 xlim=c(.4,.6),add=F); source("IECI_mc.R");plot(function(x)ECI_mc(x,x0=x0,model=kmi,y0.sample=100),col='red',                       xlim=c(.4,.6),add=T); plot(function(x)   ECI(x=matrix(x,ncol=1),x0=x0,model=kmi),col='green',lty=2,type='l',     xlim=c(.4,.6),add=T);abline(v=X); abline(v=x0,col='black');
#'
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)EI(x,model=kmi),dim=1,npoints=501,xlim=c(.4,.6)); abline(v=X); abline(v=x0,col='red');
#' @test x0=.59; set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); plot(function(x)EI(x,model=kmi),xlim=c(.4,.6),col='green'); xx=seq(f=.4,t=.6,l=501);lines(xx,ECI_mc(xx,x0=x0,model=kmi,y0.sample=100),col='blue'); lines(xx,ECI(matrix(xx,ncol=1),x0=x0,model=kmi),col='red'); abline(v=X); abline(v=x0,col='red');
#' @test x0=.55; set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); plot(function(x)EI(x,model=kmi),xlim=c(.4,.6),col='green'); xx=seq(f=.4,t=.6,l=501);lines(xx,ECI_mc(xx,x0=x0,model=kmi,y0.sample=100),col='blue'); lines(xx,ECI(matrix(xx,ncol=1),x0=x0,model=kmi),col='red'); abline(v=X); abline(v=x0,col='red');
#' @test x0=.50; set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); plot(function(x)EI(x,model=kmi),xlim=c(.4,.6),col='green'); xx=seq(f=.4,t=.6,l=501);lines(xx,1E3*ECI_mc(xx,x0=x0,model=kmi,y0.sample=1000),col='blue'); lines(xx,1E3*ECI(matrix(xx,ncol=1),x0=x0,model=kmi),col='red'); abline(v=X); abline(v=x0,col='red');
#' @test x0=.4990; set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); plot(function(x)EI(x,model=kmi),xlim=c(.4,.6),col='green'); xx=seq(f=.4,t=.6,l=501);lines(xx,5E2*ECI_mc(xx,x0=x0,model=kmi,y0.sample=1000),col='blue'); lines(xx,5E2*ECI(matrix(xx,ncol=1),x0=x0,model=kmi),col='red'); abline(v=X); abline(v=x0,col='red');
#'
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); x=matrix(seq(f=0,t=1,l=101),ncol=1);eci=ECI(x=matrix(runif(100),ncol=1),x0=x,model=kmi)
#'
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)ECI(x,x0=0.5,model=kmi),dim=1,npoints=500,xlim=c(.4,.6)); abline(v=X); abline(v=0.5,col='red'); DiceView::sectionview.fun(function(x)ECI_mc(x,x0=0.5,model=kmi),dim=1,npoints=500,add=TRUE,col='green',xlim=c(.4,.6));
#'
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); ECI(x=c(.5,.5),x0=c(.6,.6),model=kmi)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); ECI(x=matrix(runif(10),ncol=2),x0=c(.6,.6),model=kmi)
ECI <- function(x, x0, model, kn=NULL, precalc.data=NULL,Yx=NULL) {  
  if (!is.matrix(x)) x = matrix(x,nrow=1)  
  if (!is.atomic(x0)) x0 = unlist(x0)
  if (!is.matrix(x0)) x0 = t(as.matrix(x0))    
  l0=nrow(x0)
  
  x.old=x
  x.eq.x0 = row.equals(x,x0)
  if (length(x.eq.x0)>0) {warning("skipping x[]: ",paste(collapse = " ",x.eq.x0)); x = x[-x.eq.x0,]}
  
  if (!is.matrix(x)) x = matrix(x,ncol=ncol(x.old))  
  
  if(length(x)>0) {
    l=nrow(x)
    
    if (is.null(Yx)) 
      Yx = predict_nobias_km(object = model, newdata = x,  cov.compute = FALSE, se.compute = TRUE,low.memory=TRUE)
    else 
      if (length(x.eq.x0)>0) Yx=list(mean=Yx$mean[-x.eq.x0],sd=Yx$sd[-x.eq.x0])
    
    my = Yx$mean
    sdy = Yx$sd
    
    Yx0 = predict_nobias_km(object = model, newdata = x0,  cov.compute = TRUE, se.compute = TRUE) 
    my0 = Yx0$mean
    sdy0 = Yx0$sd
    Fy0 = Yx0$F.newdata
    cy0 = Yx0$c
    #covy0 = Yx0$cov
    
    # phi := x -> exp(-x^2/2)
    #
    #             /x
    # Phi := x -> | exp(-t^2/2) dt
    #             /-oo
    #
    # EI_m(x) = 1/sqrt(2*pi) * ( (m-Y(x)$mean) * Phi( (m-Y(x)$mean)/Y(x)$sd ) + Y(x)$sd * phi( (m-Y(x)$mean)/Y(x)$sd ) )
    #
    
    m = min(model@y)
    
    # @ref Corrected Kriging update formulae for batch-sequential data assimilation, Clément Chevalier
    # Y_{X+x0}(x)$mean = Y_X(x)$mean + lambda_{X+x0}(x)' * (y0 - Y_X(x0)$mean)
    # Y_{X+x0}(x)$sd2 = Y_X(x)$sd2 - lambda_{X+x0}(x)' * [cov(x0,X)]^-1 * lambda_{X+x0}(x) 
    #
    
    if (is.null(kn)) {
      if (is.null(precalc.data)) 
        precalc.data <- precomputeUpdateData(model,x) 
      else 
        if (length(x.eq.x0)>0) 
          precalc.data=list(Kinv.c.olddata=precalc.data$Kinv.c.olddata[,-x.eq.x0],Kinv.F=precalc.data$Kinv.F,first.member=precalc.data$first.member[-x.eq.x0,])
      kn <- computeQuickKrigcov(model,x,x0,precalc.data, Fy0 , cy0)
    } else 
      if (length(x.eq.x0)>0) kn=kn[-x.eq.x0]
    
    # non vectorized version (ie. one x0):
    #lambdax.new = kn/sdy0^2
    #sdy.new = sqrt(pmax(0, sdy * sdy - lambdax.new * lambdax.new * sdy0^2))
    
    # vectorized version (ie. many x0):
    lambdax.new =  kn %*% (diag(x = 1/sdy0^2,ncol=l0,nrow=l0))
    sdy.new = matrix(sqrt(pmax(0, (matrix(sdy * sdy,nrow=l,ncol=l0) - (lambdax.new * kn)))),ncol=l0,nrow=l)
    
    sdynew.eq.0 = which(sdy.new==0)
    sdy.new[sdynew.eq.0] = 1
    
    m_lim = rep((m-my0)/sdy0,each=l)
    
    #                          /m                                                      /+oo
    # ECI[x | Y.new(x0)=Y0] := |    EI_y0.new[x | Y.new(x0)=y0] d phi((y0-my0)/sdy0) + |    EI_m.new[x | Y.new(x0)=y0] d phi((y0-my0)/sdy0) 
    #                          /-oo                                                    /m
    #
    #                       := ECI_y0                                                + ECI_m 
    
    #                             /(m-Y(x0)$mean)/Y(x0)$sd
    # ECI_y0[x] =  Y.new(x)$sd    |    ( (a+b*z) * Phi(a+b*z) + phi(a+b*z) ) * exp(-z^2/2) dz
    #                             /-oo
    #
    # a :=            ( Y(x0)$mean - Y(x)$mean ) / sqrt( Y_X(x)$sd2 - lambda_{X+x0}(x)' * [cov(x0,X)]^-1 * lambda_{X+x0}(x) )
    # b :=  Y_X(x0)$sd * ( 1 - lambda_{X+x0}(x)' ) / sqrt( Y_X(x)$sd2 - lambda_{X+x0}(x)' * [cov(x0,X)]^-1 * lambda_{X+x0}(x) )
    
    a_y0 = -outer(my,my0,"-") / sdy.new   #non-vec: ( my0 - my ) / sdy.new
    b_y0 = matrix(sdy0,l,l0,T) * ( 1 - lambdax.new ) / sdy.new   #non-vec: sdy0 * ( 1 - t(lambdax.new) ) / sdy.new
    ECI_y0_int_abPhiabphi = int_abPhiabphi(upper=m_lim,lower=-Inf,a=array(a_y0),b=array(b_y0)) 
    ECI_y0_int_phiabphi   =   int_phiabphi(upper=m_lim,lower=-Inf,a=array(a_y0),b=array(b_y0))
    ECI_y0 = ECI_y0_int_abPhiabphi + ECI_y0_int_phiabphi
    ECI_y0 = matrix(ECI_y0,l)
    #print(ECI_y0_int_abPhiabphi)
    #print(ECI_y0_int_phiabphi)
    #print(ECI_y0)
    
    #                            /+oo
    # ECI_m[x] =  Y.new(x)$sd    |    ( (a+b*z) * Phi(a+b*z) + phi(a+b*z) ) * exp(-z^2/2) dz
    #                            /(m-Y(x)$mean)/Y(x)$sd
    #
    # a :=              ( m - Y(x)$mean ) / sqrt( Y_X(x)$sd2 - lambda_{X+x0}(x)' * [cov(x0,X)]^-1 * lambda_{X+x0}(x) )
    # b := - Y_X(x0)$sd * lambda_{X+x0}(x)' / sqrt( Y_X(x)$sd2 - lambda_{X+x0}(x)' * [cov(x0,X)]^-1 * lambda_{X+x0}(x) )
    
    a_m = ( m - my ) / sdy.new   # = matrix(m-my,l,l0)/sdy.new
    b_m = - lambdax.new * matrix(sdy0,l,l0,T)/ sdy.new   # non-vec: ( - t(lambdax.new) * sdy0 ) / sdy.new
    ECI_m_int_abPhiabphi = int_abPhiabphi(lower=m_lim,upper=Inf,a=array(a_m),b=array(b_m)) 
    ECI_m_int_phiabphi   =   int_phiabphi(lower=m_lim,upper=Inf,a=array(a_m),b=array(b_m))
    ECI_m = ECI_m_int_abPhiabphi + ECI_m_int_phiabphi
    ECI_m = matrix(ECI_m,l)
    #print(ECI_m_int_abPhiabphi)
    #print(ECI_m_int_phiabphi)
    #print(ECI_m)
    
    #print(sdy.new)
    
    sdy.new[sdynew.eq.0] = 0
    
    if (length(x.eq.x0)>0) {
      warning("Some x match x0");
      #print(length(ECI_y0))
      eci = insert.at.row(sdy.new * (ECI_y0 + ECI_m), x.eq.x0, matrix(0,length(x.eq.x0),l0)) / (2*pi)
      #print(length(eci))
    } else  eci = sdy.new * (ECI_y0 + ECI_m) / (2*pi)
  } else {
    warning("All x are matching an x0")
    eci = matrix(EI(x,model),nrow(x),l0)
  }
  x0.eq.X = row.equals(x0,model@X)
  if(length(x0.eq.X)>0) {#eci0 <<- eci
    warning(paste(collapse=";",length(EI(x,model))))
    warning(nrow(x))
    warning(paste(collapse=",",x0.eq.X))
    warning(paste(collapse=" ",dim(matrix(EI(x,model),nrow = nrow(x),ncol = length(x0.eq.X)))))
    warning(paste(collapse=",",dim(eci[(1:nrow(eci))[-x.eq.x0],x0.eq.X])))
    
    eci[(1:nrow(eci))[-x.eq.x0],x0.eq.X] =  matrix(EI(x,model),nrow = nrow(x),ncol = length(x0.eq.X))
  }
  #if(length(x0.eq.X)>0) eci1 <<- eci
  
  eci
  
  #list(eci=eci,a_y0=a_y0,b_y0=b_y0,ECI_y0=ECI_y0*sdy.new,a_m=a_m,b_m=b_m,ECI_m=ECI_m*sdy.new,sdy.new=sdy.new,
  #     ECI_Phi=matrix(ECI_m_int_abPhiabphi,l)*sdy.new+matrix(ECI_y0_int_abPhiabphi,l)*sdy.new,
  #     ECI_phi=matrix(ECI_y0_int_phiabphi,l)*sdy.new+matrix(ECI_m_int_phiabphi,l)*sdy.new)
}

# #debug
# set.seed(1);
# X=matrix(runif(10),ncol=1); 
# y=1-sin(2.8*X); 
# kmi <- km(design=X,response=y); 
# DiceView::sectionview.km(kmi,xlim=c(.52,.58),ylim=c(0,.02))
# abline(v=X); 
# abline(h=min(y))
# 
# ECI(x=0.54,x0=.55,kmi)
# ECI(x=0.56,x0=.55,kmi)
# ECI(x=0.549,x0=.55,kmi)
# 
# x=seq(f=.52,t=.581,l=101)
# eci.list=ECI(x = matrix(x,ncol=1),x0 = 0.55,model = kmi)
# plot(x,eci.list$eci/eci.list$sdy.new,type='l',ylim=c(-100,100))
# abline(v=.55)
# abline(v=kmi@X)
# abline(h=0)
# lines(x,eci.list$ECI_y0*eci.list$sdy.new*100000,col='red')
# lines(x,eci.list$a_y0,col='red',lty=2)
# lines(x,eci.list$b_y0,col='red',lty=3)
# lines(x,eci.list$ECI_m*eci.list$sdy.new*10000,col='blue')
# lines(x,eci.list$a_m,col='blue',lty=2)
# lines(x,eci.list$b_m,col='blue',lty=3)




#' @test sobol_inside(10,2,0,1)
#' @test sobol_inside(10,1,0,1)
require(randtoolbox)
sobol_inside <- function(n,dim,lower=0,upper=1,...) {
  s = sobol(n=n,dim=dim,...)
  #s * (upper-lower) - lower
  matrix(s * (upper-lower) - lower,ncol=dim)
}


#' @param x0
#' @param model 
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); IECI(x0=.5,model=kmi)
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); IECI(x0=matrix(c(.5,.6)),model=kmi)
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)IECI(x0=x,model=kmi,lower=c(0,0),upper=c(1,1)),dim=1,npoints=100); xx=seq(f=0,t=1,l=100); lines(xx,EI(xx,model=kmi),col='red'); abline(v=X)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); IECI(x0=c(.5,.5),model=kmi,lower=c(0,0),upper=c(1,1))
#' @test set.seed(12);X=matrix(runif(20),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)IECI(x0=x,model=kmi,lower=c(0,0),upper=c(1,1)),dim=2,npoints=20); DiceView::contourview.fun(function(x)EI(x=x,model=kmi),dim=2,npoints=20,add=TRUE,col='red'); points(X);
#' @test set.seed(1);X=matrix(runif(200),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)IECI(x0=x,model=kmi,lower=c(0,0),upper=c(1,1)),dim=2,npoints=30); DiceView::contourview.fun(function(x)EI(x=x,model=kmi),dim=2,npoints=30,add=TRUE,col='red'); points(X);
#' 
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); x=matrix(seq(f=0,t=1,l=101),ncol=1);plot(x,IECI(x0=x,model=kmi,integration.points=100))
#'
#' @perf set.seed(1);X=matrix(runif(100),ncol=2); y=rowMeans(sin(pi*X));    kmi <- km(design=X,response=y); system.time(for (i in 1:1000)IECI(x0=runif(2),model=kmi))
#' @perf set.seed(1);X=matrix(runif(100),ncol=2); y=rowMeans(sin(pi*X));    kmi <- km(design=X,response=y); integration.points = sobol_inside(100,2);precalc.data = precomputeUpdateData(kmi,integration.points);system.time(for (i in 1:1000)IECI(x0=runif(2),model=kmi,integration.points=integration.points,precalc.data=precalc.data))
#' @perf set.seed(1);X=matrix(runif(100),ncol=2); y=rowMeans(sin(pi*X));    kmi <- km(design=X,response=y); integration.points = sobol_inside(100,2);precalc.data = precomputeUpdateData(kmi,integration.points);Yx = predict_nobias_km(object = kmi, newdata = integration.points,  cov.compute = FALSE, se.compute = TRUE,low.memory=TRUE); system.time(for (i in 1:1000)IECI(x0=runif(2),model=kmi,integration.points=integration.points,precalc.data=precalc.data,Yx=Yx))
#' 
#' @perf set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); Rprof();for (i in 1:1000)IECI(x0=runif(1),model=kmi);summaryRprof()
#' @perf set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); Rprof(); integration.points = sobol_inside(10,1);precalc.data = precomputeUpdateData(kmi,integration.points);Yx = predict_nobias_km(object = kmi, newdata = integration.points,  cov.compute = FALSE, se.compute = TRUE,low.memory=TRUE); for (i in 1:1000)IECI(x0=runif(1),model=kmi,integration.points=integration.points,precalc.data=precalc.data,Yx=Yx);summaryRprof()
#' @test set.seed(1);X=matrix(runif(10),ncol=1); fun=function(x){x=x*1.2;1-1/2*(sin(12*x)/(1+x)+2*cos(7*x)*x^5+0.7)};y=fun(X);  kmi <- km(design=X,response=y); sectionview.km(kmi,ylim=c(0,1)); abline(v=X); xx=seq(f=0,t=1,l=1000); n=function(x){(x-min(x))/(max(x)-min(x))}; lines(xx,(fun(xx)),type='l');lines(xx,n(IECI(x0 = xx,model = kmi,lower=0,upper=1)),type='l',col='blue');lines(xx,n(EI(xx,model=kmi)),col='red'); 
IECI <- function(x0, model, lower=0,upper=1,plot.integration=NULL,...) {  
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
  #if (is.null(precalc.data)) precalc.data = precomputeUpdateData(model,integration.points)
  Yx = opts[['Yx']]
  #if (is.null(Yx)) Yx = predict_nobias_km(object = model, newdata = integration.points,  cov.compute = FALSE, se.compute = TRUE,low.memory=TRUE) 
  
  new_IECI = colSums(ECI(x=integration.points,x0=x0,model=model,precalc.data=precalc.data,Yx=Yx) * integration.weights)
  return(new_IECI)
}


#' @test set.seed(1);sample.resampling(function(x)matrix(1,nrow=nrow(x)),lower=0,upper=1,10)
#' @test set.seed(1);sample.resampling(function(x)matrix(1,nrow=nrow(x)),lower=c(0,0),upper=c(1,1),10)
#' @test set.seed(1); s=sample.resampling(function(x)x,lower=0,upper=1,10000); mean(s$sample*1/s$sample.density)
sample.resampling <- function(density.fun,lower,upper,N,...) {
  d=length(lower)
  s = sobol_inside(n=10*N,dim=d,lower,upper,...)
  probs=density.fun(s)
  ix = sample(1:nrow(s),size=N,replace=TRUE,prob=probs/sum(probs))
  points = matrix(s[ix,],ncol=d)
  points.fun=probs[ix]
  
  return(list(sample=points,sample.density=points.fun))
}

#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); constrIECI(x0=.5,model=kmi,constraint.fun=function(x)(x<0.49),integration.points=100)
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X);    kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)constrIECI(x0=x,model=kmi,constraint.fun=function(x)(x<0.5),integration.points=100),dim=1,npoints=100); DiceView::sectionview.fun(function(x)IECI(x0=x,model=kmi,integration.points=100),dim=1,npoints=100,add=TRUE,col='green'); xx=seq(f=0,t=1,l=100); lines(xx,EI(xx,model=kmi),col='red'); abline(v=X)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); IECI(x0=c(.5,.5),model=kmi); constrIECI(x0=c(.5,.5),model=kmi,constraint.fun=NULL)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); IECI(x0=c(.5,.5),model=kmi); constrIECI(x0=c(.5,.5),model=kmi,constraint.fun=function(x)sum(x)<1)
#' @test set.seed(12);X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)constrIECI(x0=x,model=kmi,constraint.fun=function(x)sum(x)>1),dim=2,npoints=20); DiceView::contourview.fun(function(x)IECI(x0=x,model=kmi),dim=2,npoints=20,add=TRUE,col='green'); DiceView::contourview.fun(function(x)EI(x=x,model=kmi),dim=2,npoints=20,add=TRUE,col='red'); points(X);
#' 
#' @test set.seed(1);X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); x=matrix(seq(f=0,t=1,l=101),ncol=1);plot(x,IECI(x0=x, model=kmi,integration.points=100),type='l',col='red');lines(x,constrIECI(x0=x,constraint.fun=function(x)(x-.5), model=kmi,integration.points=100),col='blue')
constrIECI <- function(x0, model, constraint.fun, lower=0, upper=1, ...) {
  #cat("      constrIECI.o1\n")
  opts=list(...)  
  integration.points=opts[['integration.points']]
  integration.weights=opts[['integration.weights']]
  
  if (is.null(constraint.fun)) return(IECI(x0,model, lower, upper, integration.points=integration.points, integration.weights=integration.weights, ...))

  if (is.null(integration.weights)) integration.weights=1
  
  if (!is.matrix(x0)) x0=matrix(x0,ncol=model@d)
  
  d=ncol(x0)
  if (is.null(integration.points) || !is.matrix(integration.points)){
    if (is.null(integration.points)) {
      n=10^d
    } else if (length(integration.points)==1) {
      n=integration.points
    } else n=nrow(integration.points)
    integration.points = sobol_inside(n,d,lower,upper)
  }
  
  constr_not_checked = which(!check.constr(integration.points,constraint.fun))
  #print(constr_not_checked)
  if (length(constr_not_checked)>0) integration.points = as.matrix(integration.points[-constr_not_checked,])
  
  if(length(integration.weights)==1) 
    integration.weights = rep(integration.weights,nrow(integration.points))
  else 
    if (length(constr_not_checked)>0)  integration.weights = integration.weights[-constr_not_checked]
  if (sum(integration.weights) != 1) 
    integration.weights = integration.weights / sum(integration.weights)
  
  #plot(integration.points,col=rgb(0,0,0,1-integration.weights),pch=20 )
  #print(integration.points)
  
  precalc.data = opts[['precalc.data']]
  if (!is.null(precalc.data) && length(constr_not_checked)>0) precalc.data=list(Kinv.c.olddata=precalc.data$Kinv.c.olddata[,-constr_not_checked],Kinv.F=precalc.data$Kinv.F,first.member=precalc.data$first.member[-constr_not_checked,])
  
  Yx = opts[['Yx']]
  if (!is.null(Yx) && length(constr_not_checked)>0) Yx=list(mean=Yx$mean[-constr_not_checked],sd=Yx$sd[-constr_not_checked])
  
  new_IEI = colSums(ECI(x=integration.points,x0=x0,model=model,precalc.data=precalc.data,Yx=Yx) * integration.weights)
  #cat("      /constrIECI.o1\n")
  return(new_IEI)
}

