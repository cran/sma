##############################################################
#
# File: Rsingle.R
#
# Created by B. M. Bolstad
# Created on Nov 20,2000 
#
# last modified by bmb
# last modified on Nov 15, 2001
#
# History:
# Nov 15, 2001 Added plot.single.slides function
# Nov 20, 2000 Initial version, created by combining
#              func.churchill.s, func.newton.s
# Jan 18, 2000 Make fixes to complete integration into sma
#
#
###############################################################


#############################################
#
# file :	func.churchill.s
# aim  :	implement methods in
#		Sapir and Churchill
#
# Created:	bmb:	7 June 2000
# Last mod:     bmb:	13 Nov 2000
#
#
# History:	
#	
# 7 June 2000 - initial versions
# 13 Nov 2000 - modifications to allow
#		integration with sma
# 15 Nov 2000 - Continued modification
#
##############################################

##############################################
#
# Internal functions - normalisation functions
# are defunct in sma framework
#
#############################################

	
##############################################
# function to do churchills normalisation.
# input is green and red values
#
# g - green
# r - red
# wts - weights to get robust regression
#
##############################################
	
chur.norm.func <- function(g,r,wts=NA){
# g is green
# r is red

if (is.na(wts)){
    wts <- rep(1,length(g))
}

tmp.fit <- lm(log(r) ~ log(g),weights=wts)
tmp.param <- tmp.fit$coef

#return orthogonal residual
list(param=tmp.param,resid=cos(atan(tmp.param[2]))*resid(tmp.fit))

# misread equation from summary
#tmp.xint <- -tmp.param[1]/tmp.param[2]
#tmp.adj <- log(g) - tmp.xint
#tmp.opp <- fitted(tmp.fit)
#tmp.diff - resid(tmp.fit)*cos(atan(tmp.opp/tmp.adj))

}

##############################################
#
# Wrapper routine to allow input of A and M 
# rather than x,y or g,r to norm function
# basically just converts to x,y then calls
# standard routine
#
# A - inputs
# M - inputs
# weights - to do robust regression
#
###############################################

chur.wrapper.func <- function(A,M,wts=NA){	
	x <- 2^(A - M/2)	
	y <- 2^(A + M/2)
	chur.norm.func(x,y)
}

###############################################
#
# Routine to take input M (log(R/G)) and return
# orthogonal residuals used in Sapir and
# Churchill
#
#
###############################################
	
chur.M.to.e.func <- function(M){
	M/sqrt(2)
}


	
################################################
#
# Function to work if within boundaries basically
# so we can deal with uniform distribution
#
################################################
	
in.func <- function(x,a,b){
(x>=a) & (x<= b)
}


#################################################
#
# Fit Mixture model using EM algorithm
#
# e - orthogonal residuals
# theta - starting parameter estimates
# maxits - maximum iterations for the EM algorithm
# a - lower bounds of uniform distribution
# b - upper bounds of uniform distribution
#
# returns final theta, final posterior probabilities of being different
#	
##################################################
	
chur.em.func <- function(e,theta=c(0.5,1),maxits=50,a=0,b=1)
{
  p <- theta[1]
  s2 <- var(e)
  notdone <- TRUE
  iter <- 1
  while (notdone)
    {
      # E-step
      z <- p*(1/(b-a))*in.func(e,a,b)/(p*in.func(e,a,b)*(1/(b-a))+ (1-p)*(1/sqrt(2*pi*s2)*exp(-e^2/(2*s2))))
      
      # M-step
      p <- sum(z)/length(e)
      s2 <- sum((1-z)*e^2)/(length(e)-sum(z))
#     print(cbind(p,s2))
      iter <- iter +1
      notdone <- (iter < maxits)
    }

  theta[1] <- p
  theta[2] <- s2
  #print(theta)
  list(theta=theta,pp=z)
}

############################################
#
# Churchill mixture model pdf
#
# x - data
# theta - fitted parameters
# a - lower bound for uniform distribution
# b - upper bound for uniform distribution
#
############################################

chur.pdf <- function(x,theta,a=0,b=1)
{
	p <- theta[1]
	s2 <- theta[2]
	(1-p)*(1/sqrt(2*pi*s2)*exp(-x^2/(2*s2))) + p * in.func(x,a,b)*(1/(b-a))
}




#################################################
#
# Function to calculate the values (of M) above 
# which all points have higher posterior probability
# than specified level
#
#
#################################################

givelim <-function(pp,p,s,b,a){
  A <- p*1/(b-a)
  B <- (1-p)*1/sqrt(2*pi*s)
  sqrt(2)*sqrt(-2*log(-A*(pp-1)/(pp*B)))*sqrt(s)
}



###############################################
#
# Only the following functions should be exposed
# to the outside world
#
###############################################
	

#############################################
#
#  Function to perform Churchill Sapir on
#  a set of slides when given 
#
#
#
#############################################

stat.ChurSap <- function(RG,layout,pp=0.95,norm="p", pout=TRUE, image.id=1, ...)
{
  MA <-stat.ma(RG, layout, norm, pout=FALSE)
  stat.ChurSap.ma(MA,pp,pout,image.id,...)
	
}

###############################################
#
# Given we already have M,A (normalised as desired)
# values perform Churchill-Sapir on set of slides
#
################################################
	
stat.ChurSap.ma <- function(MA,pp,pout=TRUE, image.id=1,...)
{
  tmpM <- MA$M[,image.id]
  orthogres <- chur.M.to.e.func(tmpM)
  ind <- !is.na(orthogres)
  pptmp <- rep(NA,length(orthogres))
  orthogres <- orthogres[!is.na(orthogres)]
  A <- min(orthogres)
  B <- max(orthogres)
  #print(orthogres)
  #print(A)
  #print(B)
  em<-chur.em.func(orthogres,a=A,b=B)
#  limits <-givelim(pp,em$theta[1],em$theta[2],B,A)
  if (pout==TRUE){
    plot(MA$A[,image.id],MA$M[,image.id],cex=0.6,xlab="A",ylab="M")
    limits <-givelim(pp,em$theta[1],em$theta[2],B,A)
    abline(limits,0)
    abline(-limits,0)
  } else {
    limits <-givelim(pp,em$theta[1],em$theta[2],B,A)
    pptmp[ind] <- em$pp	
    list(limits=limits,theta=em$theta,pp=pptmp)
  }
}
	


#################################
#
# File:	 func.newton.s
# Aim:	 Implementation of Newtons method
#
# Created ???? June 2000
#
# Last modified:	 Nov 2000
#
# History: 15 Nov 2000 Initial modifications to allow 
#		       integration with sma
#	   18 Nov 2000 Adding deriv information 
#
##################################
	
###########################################################################
# Reading 
# Modify from s.mn2a
###########################################################################

#matt.newton.func <- function(name)
#{
### This based on Newton's normalization.
#  data <- read.table(name,header=T,sep=",") 
#  nspot <- nrow(data)

# Background adjustment (very simple)
#  x <- data[,2]  ## Ch1 green  
#  y <- data[,3]  ## Ch2 red
# Normalization
#  totcy3 <- sum.na( x[x>0] )
#  totcy5 <- sum.na( y[y>0] )
  
#  x <- x/totcy3
#  y <- y/totcy5
  
  
# Rescale to help with underflow problem 10^5 (does not affect shape params)
#  x <- x*100000
#  y <- y*100000
  
#  xx <- x[!is.na(x)]
#  yy <- y[!is.na(y)]

#  list(x=xx, y=yy)
#}

###########################################################################
# newton.func.rotate <- 
# Our modify newton's function.
###########################################################################
###########################################################################
# A) newton.func.rotate <- 
###########################################################################

newton.plot.rotate <- function(A, M, theta,...)
{
## Input A and M
  ind <- is.na(A) | is.na(M) | is.infinite(A) | is.infinite(M)
  A <- A[!ind]
  M <- M[!ind]
  plot(A, M,xlab="A", ylab="M", type="n",...)
  vec1 <- seq(range(A)[1], range(A)[2], length=150)
  vec2 <- seq(range(M)[1], range(M)[2], length=150)
##  theta <- fits[1,]
  logbf <- lod2(A,M,theta=theta)
  bf <- outer(vec1,vec2,"lod2",theta=theta)
  bar <- contour(vec1,vec2,bf,levels=c(0,1,2), save=TRUE, plotit=TRUE, add=TRUE,
		 labex=0, lwd=2 )
  points( A[logbf >=0], M[logbf>=0], cex=.6 , col=2,...)
  points( A[logbf < 0], M[logbf< 0], cex=.6 , col=3,...)
##  box()

}



###########################################################################
# chen.func.rotate <-
###########################################################################
chen.func.rotate <- function(A, M, err=0.01, ...)
{
## reading in A 
## reading in M
  xx <- 2^(A - M/2)
  yy <- 2^(A + M/2)
  tt <- xx/yy
  chat <- sqrt( mean( (tt-1)^2/(1+tt^2) ) )
  tmp01 <- chen.poly(chat, err=err)
 
  abline(h =  log(tmp01[1],2), lty=2, lwd=1.5, ...)
  abline(h =  log(tmp01[2],2), lty=2, lwd=1.5, ...)
}

###########################################################################
# chen.poly (Chen's method)
###########################################################################
chen.poly <- function(cv,err=.01)
        {
        # part of table 2 from Chen et al
        bar <- rbind( c(.979, -2.706, 2.911, -2.805 ),
                      c(.989, 3.082, -2.83, 28.64),
                      c(.9968, -3.496,4.462, -5.002),
                      c( .9648,4.810,-15.161,78.349) )
        if( err==.05 ){
         coef <- bar[1,]
         tmp1 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
         coef <- bar[2,]
         tmp2 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
         }
        if( err==.01 )
         {
         coef <- bar[3,]
         tmp1 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
         coef <- bar[4,]
         tmp2 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
         }
        return( c(tmp1,tmp2) )
        }


#lod <- function(x,y,theta)
#        {
        # Log_(10) posterior odds
        # x = channel 1 intensity
        # y = channel 2 intensity
        # theta = (aa,a0,nu,pp)
#        aa <- theta[1]; a0 <- theta[2]; x0 <- theta[3]
#         y0 <- x0; z0 <- x0; pp <- theta[4]
#        tmp <- log( pp ) - log(1-pp) +
#                a0*( log(x0) + log(y0) - log(z0) ) +
#                (2*aa+a0)*log(x+y+z0) -
#                (aa+a0)*( log(x+x0) + log(y+y0) ) +
#                2*lgamma(aa+a0) - lgamma(a0) - lgamma(2*aa+a0)
#        return(tmp/2.3)
 #       }

	
###########################################################################
# Calculating lod based on A vs M plot
###########################################################################
lod2 <- function(A, M, theta)
{
# A = log(x*y)/2
# M = log(y/x)
# Log_(10) posterior odds
# x = channel 1 intensity (green)
# y = channel 2 intensity (red)
# theta = (aa,a0,nu,pp)

        x <- 2^(A - M/2)
        y <- 2^(A + M/2)
        aa <- theta[1]
        a0 <- theta[2]
        x0 <- theta[3]
        y0 <- theta[3]
        z0 <- theta[3]
        pp <- theta[4]
        tmp <- log(pp) - log(1 - pp) + a0 * (log(x0) + log(y0) - log(z0)) + (2 * 
                aa + a0) * log(x + y + z0) - (aa + a0) * (log(x + x0) + log(y + 
                y0)) + 2 * lgamma(aa + a0) - lgamma(a0) - lgamma(2 * aa + a0)
        return(tmp/2.3)
}


###########################################################################
# s.em (EM algorithm)
###########################################################################
nploglik <- function(theta,xx=xx,yy=yy,zz=zz)
        {
         # xx,yy are intensities in the two channels; zz=P(b!=c|xx,yy)
         # theta=(aa,a0,nu)
         # (I'll separately optimize pp=P(zz=1); hence npl.. for partial loglik
         aa <- theta[1]; a0<-theta[2]; x0<-theta[3]; y0<- theta[3];
         z0 <- theta[3]
         n <- length(xx)
         # 
         # Complete data loglikelihood
         ll <- (aa-1)*sum(log(xx)+log(yy))  +
                sum(zz)*( 2*(lgamma(aa+a0)-lgamma(aa)-lgamma(a0) ) ) +
                sum(zz)*a0*(log(x0)+log(y0)) +
                (n-sum(zz))*( lgamma(2*aa+a0)-2*lgamma(aa)-lgamma(a0) ) +
                (n-sum(zz))*a0*log(z0) -
                (aa+a0)*sum( zz*( log(x0+xx) + log(y0+yy) ) ) -
                (2*aa+a0)*sum( (1-zz)*( log(z0+xx+yy) ) )
        return(-ll)
        }


func.em <- function(A, M, theta=c(2,1.2,2.7,.4))
{
# EM  algorithm
# input starting values, A and M
# Beta hyperparameter for p
pprior <- 2
# starting value
notdone <- TRUE
iter <- 1
x <- 2^(A - M/2)
y <- 2^(A + M/2)
ind <- is.na(x) | is.na(y) | is.infinite(x) | is.infinite(y)
xx <- x[!ind] 
yy <- y[!ind]
xx <- xx 
yy <- yy

n <- length(xx)
while( notdone )
        {
         aa <- theta[1]; a0<-theta[2]; x0<-theta[3]; y0<- theta[3];
         z0 <- theta[3]; pp <- theta[4]
        # E-step 
        tmp <- log( pp ) - log(1-pp) +
                a0*( log(x0) + log(y0) - log(z0) ) +
                (2*aa+a0)*log(xx+yy+z0) -
                (aa+a0)*( log(xx+x0) + log(yy+y0) ) +
                2*lgamma(aa+a0) - lgamma(a0) - lgamma(2*aa+a0)
        zz <- 1/( 1 + exp(-tmp) )
        # M-step
        #fit <- nlminb( start=c(2,1.2,2.7), objective=nploglik, lower=c(1,0,0),  theta=theta, xx=xx, yy=yy, zz=zz )
	#fit <- optim( par=c(theta[1],theta[2],theta[3]), fn=nploglik,lower=c(1,10^-300,10^-300),gr =nploglikderiv, method="L-BFGS-B",control=list(trace=T),xx=xx, yy=yy, zz=zz )
         fit <- optim( par=c(theta[1],theta[2],theta[3]), fn=nploglik,lower=c(1,10^-300,10^-300),gr =nploglikderiv, method="L-BFGS-B",xx=xx, yy=yy, zz=zz )

         
        # Add a prior on pp
        theta <- c( fit$par,  ( pprior + sum( zz ) )/(2*pprior+n ) )
        #print(round(theta,4) )
        iter <- iter + 1
	 notdone <- (iter < 40)
       } 
theta
}

###########################################################################
# Old function ## without rotation (based on logG vs logR plot
###########################################################################
###########################################################################
# newton.func 
###########################################################################
#
#newton.plot.func <- function(xx, yy, theta,chen=T, chen.err = 0.01)
#{
#  plot( xx, yy, log="xy", pch=".",xlab="Cy3", ylab="Cy5", type="n")
#  vec <- log10(seq(range(c(xx, yy))[1], range(c(xx,yy ))[2], length=150))
##  theta <- fits[1,]
#  logbf <- lod(xx,yy,theta)
#  bf <- outer(10^(vec),10^(vec),"lod",theta=theta)
#  bar <- contour(vec,vec,bf,levels=c(0,1,2), save=T, plotit=T, add=T,
#		 labex=0, lwd=2 )
##  u <- bar$"0"$x
##  v <- bar$"0"$y
##  ind <- is.na(u)
##  u[ind]<- range(vec)
##  v[ind]<- range(vec)
  ## polygon(10^u,10^v,col=4)
#  points( xx[logbf >=0], yy[logbf>=0], cex=.6 , col=2)
#  points( xx[logbf < 0], yy[logbf< 0], cex=.6 , col=3)
##  box()
#  if(chen)
#    chen.func(xx, yy, err=chen.err)
#}


#####################################################################
# chen.func 
#
#
#####################################################################
#chen.func <- function(xx, yy, err=0.01)
#{
#  tt <- xx/yy
#  chat <- sqrt( mean( (tt-1)^2/(1+tt^2) ) )
#  tmp01 <- chen.poly(chat, err=err)
#  # Sandrine changed this: need log10 because plots are in log10 scale
#  abline( log(tmp01[1],10), 1, lty=2, lwd=1.5)
#  abline( log(tmp01[2],10), 1, lty=2, lwd=1.5)
#}


######################################################################
# wrapper function for SMA to perform Newtons method on data
#
#
######################################################################

stat.Newton <- function(RG,layout,norm="p",image.id=1,pout=TRUE){
	MA <-stat.ma(RG, layout, norm, pout=FALSE)
	stat.Newton.ma(MA,image.id,pout)

}


stat.Newton.ma <- function(MA,image.id=1,pout=TRUE){
	M <- MA$M[,image.id]
	A <- MA$A[,image.id]

        ind <- is.na(M) | is.na(A) 
	
	theta <- func.em(A[!ind],M[!ind])

        if (pout == TRUE){
          newton.plot.rotate(A,M,theta)
        } else {
          logodds <- rep(NA,length(M))
          logodds[!ind] <- lod2(A[!ind],M[!ind],theta)
          list(theta=theta,lod=logodds)
        }
      }


#####################################################
#
# nploglikderiv
#
# derivative of loglikelihood function hopefully
# fix broken optim (redundant, found proper scaling)
#
####################################################

	

nploglikderiv <- function(theta,xx=xx,yy=yy,zz=zz){
         # xx,yy are intensities in the two channels; zz=P(b!=c|xx,yy)
         # theta=(aa,a0,nu)
         # (I'll separately optimize pp=P(zz=1); hence npl.. for partial loglik
         aa <- theta[1]; a0<-theta[2]; x0<-theta[3]; y0<- theta[3];
         z0 <- theta[3]
	 nu <- theta[3]
         n <- length(xx)
	daa <- sum(log(xx)+log(yy)) + sum(zz)*(2*digamma(aa+a0)-2*digamma(aa)) +(n-sum(zz))*(2*digamma(2*aa+a0) - 2*digamma(aa)) - sum(zz*(log(x0+xx) + log(y0+yy))) - 2*sum((1-zz)*(log(z0+xx+yy)))
	da0 <- sum(zz)*(2*digamma(aa+a0)-2*digamma(a0)) + sum(zz)*(log(x0)+log(y0)) + (n-sum(zz))*(digamma(2*aa+a0)-digamma(a0)) + (n-sum(zz))*log(z0) - sum(zz*(log(x0+xx)+log(y0+yy))) - sum((1-zz)*log(z0+xx+yy))
	dnu <- 2*sum(zz)*a0/nu + (n-sum(zz))*a0/nu - (aa+a0)*sum(zz*(1/(nu+xx) + 1/(nu+yy))) - (2*aa+a0)*sum((1-zz)/(nu+xx+yy))
	
	c(-daa,-da0,-dnu)
}




chen.plot.rotate <- function(A, M,pout=TRUE){
  ind <- is.na(A) | is.na(M) | is.infinite(A) | is.infinite(M)
  A <- A[!ind]
  M <- M[!ind]
  if (pout==TRUE){
    plot(A, M,cex=0.6,xlab="A", ylab="M")
    chen.func.rotate(A, M, err=0.01, col=4)
    chen.func.rotate(A, M, err=0.05, col=5)
  } else {
      xx <- 2^(A - M/2)
      yy <- 2^(A + M/2)
      tt <- xx/yy
      chat <- sqrt( mean( (tt-1)^2/(1+tt^2) ) )
      tmp01 <- chen.poly(chat, err=0.01)
      tmp05 <- chen.poly(chat, err=0.05)
      list(lower01=log(tmp01[1],2),upper01=log(tmp01[2],2),lower05=log(tmp05[1],2),upper05=log(tmp05[2],2))
  } 
}



stat.Chen <- function(RG,layout,norm="p",image.id=1,pout=TRUE){
  MA <-stat.ma(RG, layout, norm, pout=FALSE)
  stat.Chen.ma(MA,image.id,pout)
}


stat.Chen.ma <- function(MA,image.id,pout=TRUE){
  chen.plot.rotate(MA$A[,image.id],MA$M[,image.id],pout) 
}



plot.single.slide <- function(x,layout,norm="p",image.id=1,...){
  #RG <- x
  MA <- stat.ma(x, layout, norm, pout = FALSE)
  Newton <- stat.Newton.ma(MA,image.id, pout=FALSE)
  ChurSap <- stat.ChurSap.ma(MA,pp=0.95,pout=FALSE,image.id)
  ChurSap2 <- stat.ChurSap.ma(MA,pp=0.99,pout=FALSE,image.id)
  Chen <- stat.Chen.ma(MA,image.id,pout=FALSE)
  newton.plot.rotate(MA$A[,image.id],MA$M[,image.id],Newton$theta,...)
  abline(h=Chen$lower01,col=4,lty=2,lwd=2)
  abline(h=Chen$upper01,col=4,lty=2,lwd=2)
  abline(h=Chen$lower05,col=5,lty=2,lwd=2)
  abline(h=Chen$upper05,col=5,lty=2,lwd=2)
  abline(h=ChurSap$limit,col=6,lty=3,lwd=2)
  abline(h=-ChurSap$limit,col=6,lty=3,lwd=2)
  abline(h=ChurSap2$limit,col=9,lty=3,lwd=2)
  abline(h=-ChurSap2$limit,col=9,lty=3,lwd=2)
}


