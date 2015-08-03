# Integrated Expected Conditional Improvement

These R files contains an analytical implementation of the criterion "IECI" proposed for Optimization under unknown constraints (2010) by [Robert Gramacy and Herbert Lee](http://arxiv.org/pdf/1004.4027v2.pdf). 

This implementation is based on [DiceKriging](https://cran.r-project.org/web/packages/DiceKriging/index.html) R package.

## R usage

Objective function


```r
fun=function(x){-(1-1/2*(sin(12*x)/(1+x)+2*cos(7*x)*x^5+0.7))}

plot(fun,type='l',ylim=c(-1.2,0.2),ylab="y")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png) 

```r
X=data.frame(X=matrix(c(.0,.33,.737,1),ncol=1))
y=fun(X)
```

Expected Improvement (EI)


```r
par(mfrow=c(2,1))

library(DiceKriging)
set.seed(123)
kmi <- km(design=X,response=y,control=list(trace=FALSE),optim.method = 'gen')
```

```
## Warning in (function (fn, nvars, max = FALSE, pop.size = 1000,
## max.generations = 100, : NAs introduits lors de la conversion automatique
```

```r
library(DiceView)
sectionview.km(kmi,col_surf='blue',col_points='blue',xlim=c(0,1),ylim=c(-1.1,0),title = "",Xname = "x",yname="Y(x)")
plot(fun,type='l',add=T)
abline(h=min(y),col='blue',lty=2)
text(labels="m",x=0.05,y=min(y)+0.05,col='blue')
legend(bg = rgb(1,1,1,.4),0.3,0,legend = c("(unknown) objective function","conditionnal gaussian process","observations"),col=c('black','blue','blue'),lty = c(1,2,0),pch = c(-1,-1,20))

source("EI.R")
.xx=seq(f=0,t=1,l=501)
plot(.xx,EI(.xx,kmi),col='blue',type='l',xlab="x",ylab="EI(x)")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 

Expected Conditional Improvement (ECI) by Monte Carlo estimation (! slow !)


```r
par(mfrow=c(2,1))

xn=0.4

sectionview.km(kmi,col_surf='blue',col_points='blue',xlim=c(xn-0.2,xn+0.2),ylim=c(-1.1,-.7),title = "",Xname = "x")
plot(fun,type='l',add=T)
abline(h=min(y),col='blue',lty=2)
text(labels="m",x=0.35,y=min(y)-0.05,col='blue')


abline(v=xn,col='red')
text(labels="xn",x=xn+0.01,y=-0.75,col='red')

Yn=predict.km(kmi,newdata=xn,type="UK",checkNames = F)
.yy = seq(f=Yn$lower,t=Yn$upper,l=11)
#lines(x=dnorm(.xxx,yn$mean,yn$sd)/100+xn,y=.xxx,col='red')
points(x=rep(xn,length(.yy)),y=.yy,col=rgb(1,0,0,0.2+0.8*dnorm(.yy,Yn$mean,Yn$sd)*sqrt(2*pi)*Yn$sd),pch=20)
points(x=xn,y=.yy[8],col=rgb(0,1,0),pch=1)
points(x=xn,y=.yy[4],col=rgb(0,1,0),pch=1)

plot(.xx,EI(.xx,kmi),col='blue',type='l',xlab="x",ylab="ECI(xn,x)",xlim=c(xn-0.2,xn+0.2),ylim=c(0,0.05))
text(labels="EI(x)",x=.4,y=0.03,col='blue')

for (yn in .yy) {
  kmii <- km(design=rbind(X,xn),response=rbind(y,yn),control = list(trace=FALSE))
  lines(.xx,5*EI(.xx,kmii),col=rgb(1,0,0,0.2+0.8*dnorm(yn,Yn$mean,Yn$sd)*sqrt(2*pi)*Yn$sd),lty=2)
}
yn = .yy[8]
kmii <- km(design=rbind(X,xn),response=rbind(y,yn),control = list(trace=FALSE))
lines(.xx,5*EI(.xx,kmii),col=rgb(0,1,0),lty=2)
yn = .yy[4]
kmii <- km(design=rbind(X,xn),response=rbind(y,yn),control = list(trace=FALSE))
lines(.xx,5*EI(.xx,kmii),col=rgb(0,1,0),lty=6)

source("IECI_mc.R",echo = FALSE)
lines(.xx,Vectorize(function(.x)5*ECI_mc_vec.o2(.x,xn,kmi))(.xx),col='red',type='l')
text(labels="ECI(xn,x) = E[ EI(x) | yn~Y(xn)] (x50)",x=.3,y=0.045,col='red')
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 


Expected Conditional Improvement (ECI) by exact calculation (! fast !)


```r
par(mfrow=c(2,1))

sectionview.km(kmi,col_surf='blue',col_points='blue',xlim=c(xn-0.2,xn+0.2),ylim=c(-1.1,-.7),title = "",Xname = "x")
plot(fun,type='l',add=T)
abline(h=min(y),col='blue',lty=2)
text(labels="m",x=0.35,y=min(y)-0.05,col='blue')


abline(v=xn,col='red')
text(labels="xn",x=xn+0.01,y=-0.75,col='red')

Yn=predict.km(kmi,newdata=xn,type="UK",checkNames = F)
.yy = seq(f=Yn$lower,t=Yn$upper,l=11)
#lines(x=dnorm(.xxx,yn$mean,yn$sd)/100+xn,y=.xxx,col='red')

plot(.xx,EI(.xx,kmi),col='blue',type='l',xlab="x",ylab="ECI(xn,x)",xlim=c(xn-0.2,xn+0.2),ylim=c(0,0.05))
text(labels="EI(x)",x=.4,y=0.03,col='blue')

source("IECI.R",echo = FALSE)
lines(.xx,5*ECI(matrix(.xx,ncol=1),xn,kmi),col='red',type='l')
```

```
## Warning in ECI(matrix(.xx, ncol = 1), xn, kmi): skipping x[]: 201
```

```
## Warning in ECI(matrix(.xx, ncol = 1), xn, kmi): Some x match x0
```

```r
text(labels="ECI(xn,x) = E[ EI(x) | yn~Y(xn)] (x50)",x=.3,y=0.045,col='red')
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png) 


And finally, Integrated Expected Conditional Improvement (IECI) calculated at the point xn: 

```
> IECI_mc(xn,kmi)
[1] 0.002014966
> IECI(xn,kmi)
[1] 0.002115793
```


Overall criterions values (EI and IECI):


```r
par(mfrow=c(1,1))

sectionview.km(kmi,ylim=c(-1,0))
abline(v=X$X,lty=2)

n=function(x){(x-min(x))/(max(x)-min(x))-1}

xx=seq(f=0,t=1,l=1000)
lines(xx,(fun(xx)),type='l')

lines(xx,n(IECI(x0 = xx,model = kmi,lower=0,upper=1)),type='l',col='blue')
```

```
## Warning in ECI(x = integration.points, x0 = x0, model = model, precalc.data
## = precalc.data, : 100
```

```
## Warning in ECI(x = integration.points, x0 = x0, model = model, precalc.data
## = precalc.data, : 100
```

```
## Warning in ECI(x = integration.points, x0 = x0, model = model, precalc.data
## = precalc.data, : 1,1000
```

```
## Warning in ECI(x = integration.points, x0 = x0, model = model, precalc.data
## = precalc.data, : 100 2
```

```
## Warning in ECI(x = integration.points, x0 = x0, model = model, precalc.data
## = precalc.data, : 0,2
```

```r
text(labels="IECI(x)",x=0.6,y=0.9-1,col='blue')

lines(xx,n(EI(xx,model=kmi)),col='red')
text(labels="EI(x)",x=0.6,y=0.1-1,col='red')
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png) 



## Analytical development of IECI:


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_01.png" height="30">


### Computable form of ECI

ECI is viewed in two main parts : updated realizations where the new conditional point 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_37.png" height="30">
 is below current minimum 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_38.png" height="30">
, and updated realizations where the new conditional point 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_39.png" height="30">
 is over current minimum 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_40.png" height="30">
 (which then does not change this current minimum for EI calculation):


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_02.png" height="30">


* 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_41.png" height="30">
 part: 

<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_03.png" height="30">
 

<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_04.png" height="30">

Change variable: 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_42.png" height="30">

  -  
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_43.png" height="30">

  -  
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_44.png" height="30">

  -  
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_45.png" height="30">

  
So,

<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_05.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_06.png" height="30">

Considering that ([@KrigingUpdate]),

<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_07.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_08.png" height="30">


Finally, 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_09.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_10.png" height="30">

Given,

<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_11.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_12.png" height="30">





* 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_46.png" height="30">
 part: 

<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_13.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_14.png" height="30">

Change variable: 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_47.png" height="30">

  -  
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_48.png" height="30">

  -  
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_49.png" height="30">

  -  
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_50.png" height="30">


So,

<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_15.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_16.png" height="30">


Finally, 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_17.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_18.png" height="30">

Given,

<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_19.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_20.png" height="30">




These integrals need the following expressions to become tractable:

* 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_51.png" height="30">

* 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_52.png" height="30">


We have (see below):


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_21.png" height="30">


And


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_22.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_23.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_24.png" height="30">


It should be noticed that this exact calculation of ECI is vectorizable, which means that synchronized computations may be performed in an efficient manner.


### Detailed calculation of integrals in 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_53.png" height="30">
 and 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_54.png" height="30">


#### 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_55.png" height="30">



<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_25.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_26.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_27.png" height="30">



#### 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_56.png" height="30">



Changing variable: 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_57.png" height="30">
,

<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_28.png" height="30">


<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_29.png" height="30">


   * Then, 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_30.png" height="30">

    
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_31.png" height="30">

   * And, 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_32.png" height="30">

  We use the bi-normal cumulative density approximation [@Genz1992]: 
  
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_33.png" height="30">

  Changing variable 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_58.png" height="30">
, 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_34.png" height="30">

and thus, 

<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_35.png" height="30">
 
<img src="https://rawgit.com/IRSN/IECI/master//figure/eq_no_36.png" height="30">




