###########################################################################
# Statistics for Microarray Analysis
# Bayesian Method
#
# Date : November 21, 2000
#
# History:
#
# Authors: Ingrid Lönnstedt  and Yee Hwa (Jean) Yang.
# First written by Ingrid Lönnstedt and modify by Jean.
##########################################################################
#This method calculate a lodscore (lods) for each gene in an
#experiment, using the normalized M-values (output from stat.ma), the
#number of slides (nb), and the number of replicates for each gene
#within each slide (nw). If there are j replicates within slides, the
#vectors of M-values for each slide should be on the form M11, ...,
#M1j, M21, ...M2j, ..., Mgj, where g is the number of genes. 

#Xprep 	is a list containing means, sums of squares etc needed for the lods
#	is produced by setup.bayesian (via stat.bayesian or plot.bayesian)
#	is calculated from X, nb and nw
#	depends only on the data, not on the prior parameters
#	Once calculated there is no need for X, nb or nw

#para	are the parameters needed for the calculation of lods
#	can be produced by stat.bayesian or plot.bayesian
#	are estimated from Xprep

#Contents in Xprep:

#nw=number of duplicates within slides (default is 1)
#nb=number of duplicates between slides
#Mbar=overall mean for a gene
#SSB=sum of squares between
#SSW=sum of squares within

#Contents in para:

#v, a =parameters in the prior for the variance.
#c=parameter in the prior for the mean.
#p=proportion of genes that are differentially expressed.
#k=ratio of variances between:within slides.


######################################################
# Main Program
######################################################

######################################################/**
# \name{stat.bayesian}
# \alias{stat.bayesian}
# 
# \title{Calculates an Odds Ratio of Each Gene in a Multi-slide microarray
# Experiment}
# \description{
# This function takes the normalized expression estimates from a
# multi-slide microarray experiment (M-values output
# by stat.ma) and returns an odds ratio for each gene: log( Pr(the gene is
# differentially expressed) / Pr(the gene is not differentially expressed)
# ). The parameter estimates of the Bayesian model used, as well as some
# data structures which are useful when presenting the lodscore
# graphically are also in the output.}
# 
# \usage{
# stat.bayesian(X=NULL, nb=NULL, nw=1, Xprep=NULL, para=list(p = 0.01,
# v = NULL, a = NULL, c = NULL, k = NULL)) 
# }
# 
# \arguments{
#   \item{X}{List containing matrix of (normalized) log expression ratios
#     \eqn{M = log_2 (R/G)} (E.g. output from stat.ma())}
#   \item{nb}{Number of slides containing spots for a gene (common for all
#   genes). Unnecessary argument if nw=1.}
#   \item{nw}{Number of spots for a gene within each slide (common to all
#     genes).Default is 1.}
#   \item{Xprep}{Some data structures useful in graphical
#     presentation. (This is calculated only if not already supplied as
#     input. See details!)} 
#   \item{para}{Estimates of the parameters used in the Bayesian
#     calculations. (These are calculated only if not
#     already supplied as input. See details!)}
# }
# \details{
# 
#   Xprep and para are optional input, but they are always in the
#   output. If Xprep is supplied as input, X, nb and nw are unnecessary
#   input. A subset of the parameters in para can be specified in the
#   input, allowing the function to estimate only the others.
# 
#   Xprep is a list containing
#   
#   \item{nb}{Number of slides containing spots for a gene (common to all
#   genes).}
#   \item{nw}{Number of spots for a gene within each slide 
#             (common to all genes).}
#   \item{Mbar}{Overalll means for each gene.}
#   \item{SSB}{Sum of squares between slides for each gene.}
#   \item{SSB}{Sum of squares within slides for each gene.}
# 
#   para is a list of parameters common to all genes containing
# 
#   \item{p}{Probability that a random gene is differentially
#     expressed. Default is 0.01.
#   \item{v,a}{Parameters in the prior for the variance. If the observed
#     variances for the genes are (V), then v and a is such that
#     (va/V)~Chi2(v).}
#   \item{c}{Parameter in the prior for the mean expression ratio.}
#   \item{k}{An approximate ratio of the variances between to within slides,
#     only needed when both nw>1 and nb>1}
# 
# }
# 
# \value{
#   
#   A list of
#  
#   \item{Xprep}{Some data structures useful in graphical
#     presentation. See details!} 
#   \item{para}{Estimates of the parameters used in the Bayesian
#     calculations. See details!}
#   \item{lods}{The log odds ratio for each gene.}
# }
# 
# 
# \author{
#     Ingrid Lönnstedt \email{ingrid@math.uu.se} \cr
#     Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
# }
# 
# \seealso{\code{\link{plot.bayesian}}}
# 
# \examples{
# data(MouseArray)
# ## mouse.setup <- init.grid() 
# ## mouse.data <- init.data() ## see \emph{init.data} 
# ## mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# ## mouse.bayesian<-stat.bayesian(X=mouse.lratio)
# ## plot(mouse.bayesian$Xprep$Mbar, mouse.bayesian$lods)
# 
# }
# 
# \keyword{microarray, bayesian, lodsratio}
#*/#####################################################

stat.bayesian <- function(X=NULL, nb=NULL, nw=1, Xprep=NULL, para=list(p=0.01, v = NULL, a=NULL, c = NULL,k=NULL))
  {
    ## Input X = List(A, M) (output from stat.ma) as well as nb (and nw if >1)
    ## Xprep and para are calculated and used for calculating lods.
    ## If Xprep is given in the function input, X, nb and nw are unnecessary.

    ## Setting up Data
    if(is.null(Xprep))
      {
    	Xprep<- setup.bayesian(X=X, nb=nb, nw=nw)
      }
    nb<-Xprep$nb
    nw<-Xprep$nw
    Mbar<-Xprep$Mbar
    SSW<-Xprep$SSW
    SSB<-Xprep$SSB

    ## Setting up parameters
     if(is.null(para$v) | is.null(para$a)) 
     {
	va <- va.func(SSB/(nb-1), vstart=list(v0=0.1, vn=nb+3, vstep=0.1), astart=list(a0=0.001, an=0.05, astep=0.001))
	para$v<-va$v
	para$a<-va$a
     }

     if(is.null(para$k)) 
     {
	if(is.null(SSW)) para$k<-0
	else para$k<-median(SSB/(nb-1)/SSW*nb*(nw-1))
     }

     if(is.null(para$c)) para$c<-c.min(para=para, Xprep=Xprep)

     lods<-lods.func(Xprep, para)
     list(Xprep=Xprep, lods=lods,para=para)
}

######################################################/**
# \name{plot.bayesian}
# \alias{plot.bayesian}
# 
# \title{Plots an Odds Ratio of Each Gene in a Multi-slide microarray
# Experiment}
# \description{
# This function takes the normalized expression estimates from a
# multi-slide microarray experiment (M-values output
# by stat.ma) and plots an odds ratio for each gene: log( Pr(the gene is
# differentially expressed) / Pr(the gene is not differentially expressed)
# ) vs the gene-specific average M-value. Alternatively, the output of
# stat.bayesian() can be given as input. The resulting plot is the same,
# but tedious calculations don't have to be done all over again.}
# 
# \usage{
# plot.bayesian(X=NULL, nb=NULL, nw=1, lods=NULL, Xprep=NULL, para=list(p
# = 0.01, v = NULL, a = NULL, c = NULL, k = NULL))
# }
# 
# \arguments{
#   \item{X}{List containing matrix of (normalized) log expression ratios
#     \eqn{M = log_2 (R/G)} (E.g. output from stat.ma())}
#   \item{nb}{Number of slides containing spots for a gene (common for all
#   genes). Unnecessary argument if nw=1.}
#   \item{nw}{Number of spots for a gene within each slide (common to all
#     genes).Default is 1.}
#   \item{lods}{The log odds ratio for each gene. Not necessary input. See
#     details!} 
#   \item{Xprep}{Some data structures useful in this graphical
#     presentation. Not necessary input. See details!} 
#   \item{para}{Estimates of the parameters used in the Bayesian
#     calculations. Not necessary input. See details!}
# }
# \details{
# 
#   lods, Xprep and para are optional input. If Xprep is supplied as
#   input, then X, nb and nw are unnecessary input. If lods is supplied, para
#   is unnecessary. A subset of the parameters
#   in para can be specified in the input, allowing the function to
#   estimate only the others. 
# 
#   Xprep is a list containing
#   
#   \item{nb}{Number of slides containing spots for a gene (common to all
#   genes).}
#   \item{nw}{Number of spots for a gene within each slide (common to
#   all genes).} 
#   \item{Mbar}{Overalll means for each gene.}
#   \item{SSB}{Sum of squares between slides for each gene.}
#   \item{SSB}{Sum of squares within slides for each gene.}
# 
#   para is a list of parameters common to all genes containing
# 
#   \item{p}{Probability that a random gene is differentially
#     expressed. Default is 0.01.
#   \item{v,a}{Parameters in the prior for the variance. If the observed
#     variances for the genes are (V), then v and a is such that
#     (va/V)~Chi2(v).}
#   \item{c}{Parameter in the prior for the mean expression ratio.}
#   \item{k}{An approximate ratio of the variances between to within slides,
#     only needed when both nw>1 and nb>1}
# 
# }
# 
# \author{
#     Ingrid Lönnstedt \email{ingrid@math.uu.se} \cr
#     Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
# }
# 
# \seealso{\code{\link{stat.bayesian}}}
# 
# \examples{
# data(MouseArray)
# ## mouse.setup <- init.grid() 
# ## mouse.data <- init.data() ## see \emph{init.data} 
# ## mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# #Alternative 1
# ## mouse.bayesian<-stat.bayesian(X=mouse.lratio)
# ## plot.bayesian(Xprep=mouse.bayesian$Xprep, lods=mouse.bayesian$lods)
# 
# #Alternative 2
# ## plot.bayesian(X=mouse.lratio)
# 
# #My changes
# ## my.para<-mouse.bayesian$para
# ## my.para$p<-0.05
# ## plot.bayesian(Xprep=mouse.bayesian$Xprep, para=my.para)
# 
# }
# 
# \keyword{microarray, bayesian, lodsratio}
#*/#####################################################
 

plot.bayesian <- function(X=NULL, nb=NULL, nw=1, lods=NULL, Xprep=NULL, para=list(p=0.01, v = NULL, a=NULL, c = NULL,k=NULL))

  {
    ## Input X = List(A, M) (output from stat.ma) as well as nb (and nw if >1)
    ## Xprep and para are calculated and used for plotting lods.
    ## If Xprep is given in the function input, X, nb and nw are unnecessary.
    
    ## Setting up Data
    if(is.null(Xprep))
      {
    	Xprep<- setup.bayesian(X=X, nb=nb, nw=nw)
      }
    nb<-Xprep$nb
    nw<-Xprep$nw
    Mbar<-Xprep$Mbar
    SSW<-Xprep$SSW
    SSB<-Xprep$SSB

    ## Setting up parameters
    if(is.null(lods))
    {
       if(is.null(para$v) | is.null(para$a)) 
	{
		va <- va.func(SSB/(nb-1), vstart=list(v0=0.1, vn=nb+3, vstep=0.1), astart=list(a0=0.001, an=0.05, astep=0.001))
		para$v<-va$v
		para$a<-va$a
       	}
       if(is.null(para$k)) 
	{
		if(is.null(SSW)) para$k<-0
		else para$k<-median(SSB/(nb-1)/SSW*nb*(nw-1))
	}
        if(is.null(para$c)) para$c<-c.min(para=para, Xprep=Xprep)
    
       lods<-lods.func(Xprep, para)
    }

    ## Plotting
    	plot(Mbar, lods, xlab="Mean", ylab="lodsratio", main="Lodsratio vs Mean", type="n")
    	text(Mbar, lods, cex=1)
}


    
######################################################
#Functions
######################################################

setup.bayesian <- function(X, nb=NULL, nw=1)
  {
    if (nw == 1)
    {
	nb<-ncol(X$M)
	SSW<-NULL
	SSB<-apply(X$M,1,var.na)*(nb-1)
	Mbar<-apply(X$M,1,mean.na)
    }

    if (nw > 1)
    {
	
	if (nb > 1)
	{
  	  Mtmp<-NULL
	  for (i in 1:nb)
	  {
	    for (j in 1:nw)
	    {
	      Mtmp<-cbind(Mtmp,X$M[seq(j,nrow(X$M),nw),i])
	    }
	  }
	  Mbar<-apply(Mtmp,1,mean.na)

	  Mslide<-NULL
	  SSW<-rep(0,nrow(X$M)/nw)
	  SSB<-rep(0,nrow(X$M)/nw)
	  for (i in 1:nb)
	  {
	    Mslide<-cbind(Mslide,apply(Mtmp[,((i-1)*nw+1):(i*nw)],1,mean.na))
	    SSW<-SSW+apply((Mtmp[,((i-1)*nw+1):(i*nw)]-Mslide[,i])**2,1,sum.na)	
	    SSB<-SSB+nw*((Mslide[,i]-Mbar)**2)
	  }
	}

	if (nb == 1)
	{

  	  Mtmp<-NULL

          for (j in 1:nw)
          {
	    Mtmp<-cbind(Mtmp,X$M[seq(j,length.na(X$M),nw)])
	  }

	  SSB<-apply(Mtmp,1,var.na)*(nw-1)
	  nb<-nw
	  nw<-1
	  SSW<-NULL
	  Mbar<-apply(Mtmp,1,mean.na)
	}
	
    }	
	list(Mbar=Mbar, SSB=SSB, SSW=SSW, nb=nb, nw=nw)
  }

#Lods calculates the logodds ratio.

lods.func<-function(Xprep=list(Mbar=Mbar, SSB=SSB, SSW=SSW, nb=nb, nw=nw), para=list(p=p, v=v, a=a, k=k)){
   Mbar <- Xprep$Mbar            #overall means
   nb<-Xprep$nb
   nw<-Xprep$nw
   SSB <- Xprep$SSB              #sums of squares between slides
   SSW <- Xprep$SSW              #sums of squares within slides  
   if(is.null(Xprep$SSW)) SSW <- rep(0, length.na(SSB))

   p <- para$p
   v <- para$v
   a <- para$a
   c <- para$c
   k <- para$k

   odds1<-p/(1-p)
   odds2<-c/((nb*nw) + c)
   odds3<-(v * a + SSB + nb * nw * (Mbar^2)+k*SSW)/(v*a+SSB+ nb * nw *(Mbar^2) + k * SSW-(nb * nw * Mbar)^2/(c + nb * nw))
   odds<-odds1*(odds2**(1/2))*(odds3**((nb*nw+v)/2-1))
   log(odds)
 }


######################################################

#va.func fits a chi^2 distribution to the variance (finds the parameters v (df) and a(scale)). It only checks the values of v and a in a range of values. The range of these might need to be changed.
#va2.func also fits a chi^2 distribution to the variance, but it does it by iteration (nlm function ~ Newton Raphson). va2.func tends to get stuck in local maximas if just used starting at the estimates from the method of moments (default).
#va.func is good, but to get a finer estimate (not really necessary!), one could use va2.func with starting values equal to the estimates from va.func. 


sq.func <- function(v, a, MSB,nclass=100, pout = F, ...)
  {
    index <- seq(1, length.na(MSB), round(length.na(MSB)/nclass))
    hst<-hist((log(v*a/MSB))[index], plot=FALSE, freq=F, nclass=10)
    theo<-dchisq(exp(hst$mids), df=v)
    sum.na((hst$density-theo)**2)
  }


va.func <- function(MSB, vstart=list(v0=0.1, vn=nb+3, vstep=0.1), astart=list(a0=0.0001, an=0.05, astep=0.0001), pout=F, aseq=NULL, vseq=NULL){

  if(is.null(vseq)) vseq<-seq(vstart$v0, vstart$vn, vstart$vstep)
  if(is.null(aseq)) aseq<-seq(astart$a0, astart$an, astart$astep)
  
  sq<-matrix(0, ncol=length.na(vseq), nrow=length.na(aseq))

  i<-0
  for (a in aseq){
    i<-i+1
    j<-0
    for (v in vseq){
      j<-j+1
      sq[i,j] <- sq.func(v, a, MSB, nclass=100)
    }
  }  
  
  minsq<-min(sq)
  for (i in 1:length.na(aseq)){
    for (j in 1:length.na(vseq)){
      if (sq[i,j]==minsq){
        minv<-j
        mina<-i
      }
    }
  }
  v<-vseq[minv]
  a<-aseq[mina]

  if(pout){
    hst<-hist(log(v*a/MSB),freq=F, plot=F)
    x<-exp(hst$mids)
    fx<-dchisq(x, v)
    hst<-hist(log(v*a/MSB),freq=F, ylim=c(0,max(fx,hst$density)))
    points(hst$mids,fx, type="p", pch=16)

  }

  list(v=v, a=a)
}

sq2.func <- function(va=c(v,a))
  { 
      sq.func(va[1], va[2], MSB=MSB, nclass=100)
  }

va2.func <- function(MSB, vstart=list(v0=NULL, vn=100, vstep=0.1), astart=list(a0=NULL, an=100, astep=0.001)){

  if(is.null(vstart$v0) | is.null(astart$a0)){
    Sbar<-mean.na(MSB)
    Sbar2<-mean.na(MSB)**2
    S2bar<-mean.na(MSB**2)
  
    if(is.null(vstart$v0)){
      ## estimate using MOM
      vstart$v0<-(4*S2bar-2*Sbar2)/(S2bar-Sbar2)
    }
    
    if(is.null(astart$a0)){
      ## estimate using MOM
      astart$a0<-S2bar*Sbar/(2*S2bar-Sbar2)
    }
  }
  tmp<-nlm(sq2.func, p=c(vstart$v0, astart$a0), stepmax=max(vstart$vn,astart$an), steptol=min(vstart$vstep, astart$astep))
  list(v=tmp$estimate[1], a=tmp$estimate[2])
}

######################################################
#Lf calculates the (-1)*likelihood of the data as a function of c.
#This function is not used by any main functions.

lf.func <-function(point){
c<-point

f1<-gamma((v+nw*nb)/2)/gamma(v/2)
f2<-pi
f3<-c/(c+nw*nb)
f4<-k
if (k==0){
	f4<-1
	SSW<-rep(0, length.na(SSB))
}
f5<-v*a
f6<-v*a+nw*nb*Mbar**2+SSB+k*SSW-(nw*nb*Mbar)**2/(c+nw*nb)
f7<-v*a+nw*nb*Mbar**2+SSB+k*SSW

fg1<-f1*(f2**(-nw*nb/2))*(f3**(1/2))*(f4**((nw-1)*nb/2))*(f5**(v/2-1))*(f6**(-((v+nw*nb)/2-1)))
fg2<-f1*(f2**(-nw*nb/2))*(f4**((nw-1)*nb/2))*(f5**(v/2-1))*(f7**(-((v+nw*nb)/2-1)))
fg<-log(p*fg1+(1-p)*fg2)
f<-sum.na(fg)
-f
}

######################################################
#dc.min and c.min are used to estimate c.
#If there is no c within (0,1], the result is NULL
#It is possible to change the range to for example (1,2] when calling c.min.

dc.min<-function(c, para, Xprep)
 {
  para$c<-c
  if (is.null(Xprep$SSW)) Xprep$SSW<-rep(0,length.na(Xprep$SSB))
  l<-stat.bayesian(Xprep=Xprep, para=para)$lods
  T<-(1:(nrow(X$M)/Xprep$nw))[rank(l)>(length.na(l)-round(length.na(l)*para$p))]

  #Posterior estimates of mean and variance
  tau<-(Xprep$nb*Xprep$nw*+para$v)/(para$v*para$a+Xprep$nb*Xprep$nw*Xprep$Mbar[T]^2+ Xprep$SSB[T]+ para$k*Xprep$SSW[T]-(Xprep$nb*Xprep$nw*Xprep$Mbar[T])^2/(Xprep$nb*Xprep$nw+c))
  mu2<-(Xprep$nb*Xprep$nw*Xprep$Mbar[T]/(Xprep$nb*Xprep$nw+c))^2-1/c/tau+1/tau/(Xprep$nb*Xprep$nw+c)

  sum.na(mu2*c*tau)/length.na(T)-1

}

c.min<-function(u=10^(-7),l=1, para, Xprep, ndigit=5)
 {
 if (dc.min(u,para=para,Xprep=Xprep) * dc.min(l,para=para,Xprep=Xprep) <0){
  if (abs(u-l)> 10^(-ndigit))
  {
    new<-0.5*(u+l)
	
    if (dc.min(u,para=para,Xprep=Xprep) * dc.min(new,para=para,Xprep=Xprep) <0 )	c.min(u=u, l=new, para=para, Xprep=Xprep)
    else  c.min(u=new, l=l, para=para, Xprep=Xprep)
  }	   
  else round(0.5*(u+l),digits=ndigit)
 }
 }
