#' @test plot(phi,xlim=c(-3,3))
phi=function(t) exp(-t^2/2)

#' / x
#' |     exp(-t^2/2) dt
#' / -oo
#' 
#' @ref Phi.approx = function(t,lower=-10) integrate(f=phi,lower = lower,upper=t)$value
#' @test Phi(1) - Phi.approx(1)
#' @test plot(Phi,xlim=c(-3,3))
Phi=function(t) sqrt(2*pi) * pnorm(t)

#'              / x
#' 2/sqrt(pi) * |   exp(-t^2) dx
#'              / 0
#'
#' @ref erf.approx = function(x) 2*(integrate(f=function(t)exp(-t^2), lower=0 , upper=x))$value/sqrt(pi)
#' @test erf(1) - erf.approx(1)
#' @test plot(erf,xlim=c(-3,3))
erf = function(x) 2*pnorm(x*sqrt(2))-1
# erf(x) = Phi(x*sqrt(2))*sqrt(2/pi)-1
# Phi(x) = (erf(x/sqrt(2))+1)*sqrt(pi/2)

#' /x    /y
#' |     |      1/(2*pi*sd_x*sd_y/sqrt(1-rho^2)) * exp(-1/(2*(1-rho^2)) * ( ((u-mean_x)/sd_x)^2 + ((v-mean_y)/sd_y)^2 - 2*rho * (u-mean_x)/sd_x * (v-mean_y)/sd_y) ) du dv
#' /-oo  /-oo
#'
#' @test x = seq(f=-6,t=6,l=50); z = matrix(NaN, 50,50); for (i in 1:50) for (j in 1:50) z[i,j] = pbnorm(x[i],0,1,x[j],0,1,0); persp(x,x,z)
library(pbivnorm)
#library(mnormt)
#' @test pbnorm(0,y=0)
#' @test pbnorm(x=-Inf,0,1,y=0,-5,1,-0.5)
#' @test pbnorm(x=Inf,0,1,y=0,5,1,-0.5)
#' @test pbnorm(x=100,mean_x=0,sd_x=1,y=0,mean_y=-10000,sd_y=1,rho=-0.95)
#' @test pbnorm(x=0,0,1,y=c(0,1),5,1,-0.5)
pbnorm <- function(x,mean_x=0,sd_x=1,y,mean_y=0,sd_y=1,rho=0) {
  l=100
  #stopifnot(length(x)==length(y))
  x.norm=(x-mean_x)/sd_x
  y.norm=(y-mean_y)/sd_y
  x.norm[which(x.norm>l)]=l
  y.norm[which(y.norm>l)]=l
  x.norm[which(x.norm<(-l))]=-l
  y.norm[which(y.norm<(-l))]=-l
  if (length(x.norm)==1 && length(y.norm)>1) x.norm=rep(x.norm,length(y.norm))
  if (length(y.norm)==1 && length(x.norm)>1) y.norm=rep(y.norm,length(x.norm))
  pbivnorm(x.norm,y.norm,rho)
  
  #varcov=diag(2)
  #varcov[1,2]=rho
  #varcov[2,1]=rho
  #pmnorm(cbind(x.norm,y.norm),varcov = varcov)
}

#todo : pbnorm.dx(x1,x2,mean_x=0,sd_x=1,y,mean_y=0,sd_y=1,rho=0)

#' @test y=5;a=10;b=10;integrate(f=Vectorize(function(u)integrate(f=Vectorize(function(t)exp(-(t+a+b*u)^2/2)),lower=-100,upper=0)$value*exp(-u^2/2)),lower=-100,upper=y)
#' @test 2*pi*pbnorm(y,0,1,0,-a,sqrt(1+b^2),-b/sqrt(1+b^2))


################ int(phi(a+b*u)*exp(-u^2/2), u = -infinity .. x) #####################

#' /u 
#' |  exp(-(a+b*x)^2/2) dt exp(-x^2/2) dx
#' /l 
#'   
#' @test plot(Vectorize(int_phiabphi.approx),xlim=c(-5,5))
int_phiabphi.approx <- function(upper,lower=-10,a=0,b=1) {
  integrate(f=function(u) phi(a+b*u)*exp(-u^2/2), lower=lower,upper=upper,abs.tol = 1E-15)$value 
}

#' @test int_phiabphi(3,-10,0,20000000000) - int_phiabphi.approx(3,-10,0,2000000000)
#' @test int_phiabphi(3,-10,10,20000000000) - int_phiabphi.approx(3,-10,10,2000000000)
#' @test int_phiabphi(3,-10,2,0) - int_phiabphi.approx(3,-10,2,0)
#' @test int_phiabphi(3,-10,2,200000) - int_phiabphi.approx(3,-10,2,200000)
#' @test int_phiabphi(3,-10,20,20) - int_phiabphi.approx(3,-10,20,20)
#' 
#' @test int_phiabphi(3,-10,2,-2) - int_phiabphi.approx(3,-10,2,-2)
#' @test int_phiabphi(3,-10,-2,-2) - int_phiabphi.approx(3,-10,-2,-2)
#' @test int_phiabphi(3,-10,-2,2) - int_phiabphi.approx(3,-10,-2,2)
#' 
#' @test int_phiabphi(c(3,4),c(-10,-11),c(-2,-3),c(2,3)) - c(int_phiabphi(3,-10,-2,2) , int_phiabphi(4,-11,-3,3))
#' 
#' @test int_phiabphi(3,-Inf,-2,2) 
#' @test int_phiabphi(Inf,0,-2,2) 
#' @test plot(int_phiabphi,xlim=c(-5,5))
int_phiabphi <- function(upper,lower=-10,a=0,b=1) {
  -phi(sqrt(a^2/(1+b^2)))/sqrt(1+b^2)*(Phi(((1+b^2)*lower+a*b)/sqrt(1+b^2)) - Phi(((1+b^2)*upper+a*b)/sqrt(1+b^2))) 
}


################ int(Phi(a+b*u)*exp(-u^2/2), u = -infinity .. x) #####################


#' /u /x
#' |  | exp(-t^2/2) dt exp(-x^2/2) dx
#' /l /-oo 
#'    
# = (2*pi) / 8 *((1 + erf(upper / sqrt(2))) ^ 2 - (1 + erf(lower / sqrt(2))) ^ 2)
#' @test int_Phiphi(3,0) == integrate(function(x)Phi(x)*phi(x),lower=0,upper=3)$value
#' @test int_Phiphi(+Inf,0)
#' @test int_Phiphi(0,-Inf)
#' @test plot(int_Phiphi,xlim=c(-5,5))
int_Phiphi <- function(upper,lower=-10) 0.5*(Phi(upper)^2 - Phi(lower)^2)

#' /u /a+b*x
#' |  | exp(-t^2/2) dt exp(-x^2/2) dx
#' /l /-oo 
#'   
#' @test int_Phiabphi.approx(3,-10,0,1) - int_Phiphi(3)
#' @test plot(Vectorize(int_Phiabphi.approx),xlim=c(-5,5))
int_Phiabphi.approx <- function(upper,lower=-10,a=0,b=1) {
  integrate(f=function(u) Phi(a+b*u)*exp(-u^2/2), lower=lower,upper=upper,abs.tol = 1E-15)$value 
}

#' @test int_Phiabphi(3,-10,0,200000000) - int_Phiabphi.approx(3,-10,0,200000000)
#' @test int_Phiabphi(3,-10,2,0) - int_Phiabphi.approx(3,-10,2,0)
#' @test int_Phiabphi(3,-10,2,2) - int_Phiabphi.approx(3,-10,2,2)
#' @test int_Phiabphi(3,-10,20,20) - int_Phiabphi.approx(3,-10,20,20)
#' @test int_Phiabphi(3,-10,2,-2) - int_Phiabphi.approx(3,-10,2,-2)
#' @test int_Phiabphi(3,-10,-2,-2) - int_Phiabphi.approx(3,-10,-2,-2)
#' @test int_Phiabphi(3,-10,-2,2) - int_Phiabphi.approx(3,-10,-2,2)
#' @test int_Phiabphi(0,-Inf,3,2)
#' @test int_Phiabphi(Inf,0,-2,2)
#' @test int_Phiabphi(upper=c(10,10),lower=c(-10,-10),a=c(-2,-2),b=c(2,0))
#' @test int_Phiabphi(upper=c(10,10,12),lower=c(-10,-10,-11),a=c(-2,-2,-3),b=c(2,0,1))
#' @test int_Phiabphi(upper=c(10,10,12),lower=-10,a=-2,b=2)
#' @test plot(int_Phiabphi,xlim=c(-5,5))
int_Phiabphi <- function(upper,lower=-10,a=0,b=1) {
  #if(b==0) return(Phi(a)*(Phi(upper)-Phi(lower)))
  #(pbnorm(upper,0,1,0,-a/abs(b),sqrt(1+b^2)/abs(b),-abs(b)/sqrt(1+b^2)) - pbnorm(lower,0,1,0,-a/abs(b),sqrt(1+b^2)/abs(b),-abs(b)/sqrt(1+b^2)) )* (2*pi)
  
  b.is.0 = unique(c(which(b==0),which(is.na(b)),which(!is.numeric(b))))
  if (length(b.is.0)>0) {
    b[b.is.0] = 1E-5
    a.orig=a
    if (length(a)>1) a[b.is.0] = 0
    upper.orig=upper
    if (length(upper)>1) upper[b.is.0] = 0
    lower.orig=lower
    if (length(lower)>1) lower[b.is.0] = 0
  }
  
  sqrt_1pb2 = sqrt(1+b^2)
  bosqrt_1pb2 = b/sqrt_1pb2
  
  if (any(abs(bosqrt_1pb2)>1)) stop("rho = ",-bosqrt_1pb2)
  
  check.length.equals.or.1(lower,upper)
  check.length.equals.or.1(upper,a,b)
  check.length.equals.or.1(lower,a,b)
  
  p=(pbnorm(upper,0,1,0,-a,sqrt_1pb2,-bosqrt_1pb2) - pbnorm(lower,0,1,0,-a,sqrt_1pb2,-bosqrt_1pb2) )* (2*pi)
  
  if (length(b.is.0)>0) 
    p[b.is.0] = Phi(a.orig[b.is.0])*(Phi(upper.orig[b.is.0])-Phi(lower.orig[b.is.0]))

  return(p)
}

check.length.equals.or.1 <- function(...) {
  all=list(...)
  l=1
  for(n in names(all)) {
    ln = length(all[[n]])
    if (l>1) stopifnot(ln != l)
    l=ln
  }
}

#' @test insert.at.1(a=c(1,2,3),pos=2,toinsert=-1)
#' @test insert.at.1(a=c(1,2,3),pos=1,toinsert=-1)
#' @test insert.at.1(a=c(1,2,3),pos=3,toinsert=-1)
insert.at.1 <- function(a, pos, toinsert){
  stopifnot(length(toinsert)==1)
  stopifnot(length(pos)==1)
  result=array(a,length(a))
  result[1:(pos-1)]=a[1:(pos-1)]
  result[pos]=toinsert
  result[(pos+1):(length(a)+1)]=a[pos:length(a)]
#   result <- vector("list",2*length(pos)+1)
#   result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos)))
#   result[c(FALSE,TRUE)] <- toinsert
#   unlist(result)
result
}

#' @test insert.at(a=c(1,2,3),pos=c(2,3),toinsert=c(-1,-2))
#' @test insert.at(a=runif(100),pos=1:10,toinsert=-runif(10))
insert.at <- function(a, pos, toinsert){
  stopifnot(length(toinsert)==length(pos))
  ix=sort(pos,index.return=T)$ix
  for (i in 1:length(toinsert))
    a = insert.at.1(a,pos[ix[i]],toinsert[ix[i]])
  a
}

#' @test insert.at.1col(m=matrix(c(1,2,3,4,5,6),2,3),col.pos=2,col.toinsert=c(-1,-1))
#' @test insert.at.1col(m=matrix(runif(101*7),7,101),col.pos=2,col.toinsert=-runif(7))
insert.at.1col <- function(m, col.pos, col.toinsert){
  stopifnot(all(col.pos<=ncol(m)))
  stopifnot(length(col.toinsert)==nrow(m))
  
  t(matrix(insert.at(m,(col.pos-1)*nrow(m)+1:nrow(m),col.toinsert),ncol(m)+1,nrow(m),byrow=T))
}

# insert.at.col <- function(m, cols.pos, cols.toinsert){
#   stopifnot(length(cols.pos)==length(cols.toinsert)/nrow(m))
#   for (i in 1:length(cols.pos))
#     m = insert.at.1col(a,cols.pos[i],cols.toinsert[,i])
#   m
# }

#' @test insert.at.1row(m=matrix(c(1,2,3,4,5,6),3,2),row.pos=2,row.toinsert=c(-1,-1))
insert.at.1row <- function(m, row.pos, row.toinsert){
  t(insert.at.1col(t(m),row.pos,t(row.toinsert)))
}

#' @test insert.at.row(m=matrix(c(1,2,3,4,5,6),3,2),rows.pos=c(2,3),rows.toinsert=matrix(runif(4),2,2))
#' @test insert.at.row(m=matrix(c(1,2,3,4,5,6),3,2),rows.pos=c(3,2),rows.toinsert=matrix(runif(4),2,2))
insert.at.row <- function(m, rows.pos, rows.toinsert){
  stopifnot(length(rows.pos)==length(rows.toinsert)/ncol(m))
  ix=sort(rows.pos,index.return=T)$ix
  for (i in 1:length(rows.pos))
    m = insert.at.1row(m,rows.pos[ix[i]],rows.toinsert[ix[i],])
  m
}



##################################################################

#' /u     /a+b*x
#' |  b*x*| exp(-t^2/2) dt exp(-x^2/2) dx
#' /l     /-oo 
#'   
#' @test int_bPhiabphi.approx(3,-10,0,1)
#' @test plot(Vectorize(int_bPhiabphi.approx),xlim=c(-5,5))
int_bPhiabphi.approx <- function(upper,lower=-10,a=1,b=1) {
  integrate(f=function(u) (b*u)*Phi(a+b*u)*exp(-u^2/2), lower=lower,upper=upper,abs.tol = 1E-15)$value
}

#' @test int_bPhiabphi(3,-10,0,0) - int_bPhiabphi.approx(3,-10,0,0)
#' @test int_bPhiabphi(3,-10,0,1) - int_bPhiabphi.approx(3,-10,0,1)
#' @test int_bPhiabphi(3,-10,0,-20000) - int_bPhiabphi.approx(3,-10,0,-20000)
#' @test int_bPhiabphi(3,-10,0,.2) - int_bPhiabphi.approx(3,-10,0,.2)
#' @test int_bPhiabphi(3,-10,0,-1) - int_bPhiabphi.approx(3,-10,0,-1)
#' @test int_bPhiabphi(3,-10,0,-2) - int_bPhiabphi.approx(3,-10,0,-2)
#' @test int_bPhiabphi(3,-10,0,-.2) - int_bPhiabphi.approx(3,-10,0,-.2)
#' @test int_bPhiabphi(3,-10,2,2) - int_bPhiabphi.approx(3,-10,2,2)
#' @test int_bPhiabphi(3,-10,3,2) - int_bPhiabphi.approx(3,-10,3,2)
#' @test int_bPhiabphi(3,-10,.3,2) - int_bPhiabphi.approx(3,-10,.3,2)
#' @test int_bPhiabphi(3,-Inf,.3,2)
#' @test int_bPhiabphi(Inf,0,.3,2)
int_bPhiabphi <- function(upper,lower=-10,a=0,b=1) {
  b*(Phi(a+b*lower)*phi(lower) - Phi(a+b*upper)*phi(upper)) + b^2*int_phiabphi(upper,lower,a,b)
}


################ int((a+b*u)*Phi(a+b*u)*exp(-u^2/2), u = -infinity .. x) #####################

#' /u         /a+b*x
#' |  (a+b*x)*| exp(-t^2/2) dt exp(-x^2/2) dx
#' /l         /-oo 
#'  
#' @test int_abPhiabphi.approx(3,-10,0,1)
#' @test plot(Vectorize(int_abPhiabphi),xlim=c(-5,5))
int_abPhiabphi.approx <- function(upper,lower=-10,a=1,b=1) {
  integrate(f=function(u) (a+b*u)*Phi(a+b*u)*exp(-u^2/2), lower=lower,upper=upper,abs.tol = 1E-15)$value
}

#' @test int_abPhiabphi(3,-10,0,0) - int_abPhiabphi.approx(3,-10,0,0)
#' @test int_abPhiabphi(3,-10,0,1) - int_abPhiabphi.approx(3,-10,0,1)
#' @test int_abPhiabphi(3,-10,0,10000) - int_abPhiabphi.approx(3,-10,0,10000)
#' 
#' @test int_abPhiabphi(3,-10,0,.2) - int_abPhiabphi.approx(3,-10,0,.2)
#' @test int_abPhiabphi(3,-10,0,-1) - int_abPhiabphi.approx(3,-10,0,-1)
#' @test int_abPhiabphi(3,-10,0,-2) - int_abPhiabphi.approx(3,-10,0,-2)
#' @test int_abPhiabphi(3,-10,0,-.2) - int_abPhiabphi.approx(3,-10,0,-.2)
#' 
#' @test int_abPhiabphi(3,-10,2,2) - int_abPhiabphi.approx(3,-10,2,2)
#' @test int_abPhiabphi(3,-10,3,2) - int_abPhiabphi.approx(3,-10,3,2)
#' @test int_abPhiabphi(3,-10,.3,2) - int_abPhiabphi.approx(3,-10,.3,2)
#' 
#' @test int_abPhiabphi(3,-Inf,3,2)
#' @test int_abPhiabphi(Inf,0,3,2)
#' 
#' @test int_abPhiabphi(c(3,4),c(-10,-11),c(.3,.4),c(2,3)) - c(int_abPhiabphi(3,-10,.3,2) , int_abPhiabphi(4,-11,.4,3))
#' 
#' @test plot(function(x)int_abPhiabphi(x,a=0,b=10),xlim=c(-15,15)); plot(Vectorize(function(x)int_abPhiabphi.approx(upper=x,a=0,b=10)),add=T,col='blue',xlim=c(-15,15))
int_abPhiabphi <- function(upper,lower=-10,a=0,b=1) {
  #abs(b)*(Phi(a+b*lower)*phi(lower) - Phi(a+b*upper)*phi(upper)) + b^2*int_phiabphi(upper,lower,a,b) +  a*int_Phiabphi(upper,lower,a,b)
  b*(Phi(a+b*lower)*phi(lower) - Phi(a+b*upper)*phi(upper) + b*int_phiabphi(upper,lower,a,b)) + a*int_Phiabphi(upper,lower,a,b)
}

