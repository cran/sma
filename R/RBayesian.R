###########################################################################
# Statistics for Microarray Analysis
# Bayesian Method
#
# Date : October 2, 2001
#
# History:
#
# Authors: Ingrid Lönnstedt  and Yee Hwa (Jean) Yang.
# Written by Ingrid with help from Jean.
##########################################################################

#stat.bayesian calculates a lodscore (lods) for each gene in an experiment, using the normalized M-values (output from stat.ma), the number of slides (nb), and the number of replicates for each gene within each slide (nw). If there are j replicates within slides, the vectors of M-values for each slide should be on the form M11, ..., M1j, M21, ...M2j, ..., Mgj, where g is the number of genes.

#stat.bay.est calculates a lodscore (lods) for each gene in an experiment, using independent, sufficient statistics for the effect and its variance of each gene. 

#plot.bayesian plots the results of any of the above, highlighting genes which meet userdefined criteria.

######################################################
# Main Program
######################################################

stat.bayesian <- function(M=NULL, nb=NULL, nw=1, Xprep=NULL, para=list(p=0.01, v = NULL, a=NULL, c = NULL,k=NULL))
  {
    ## Input M (output from stat.ma) as well as nb (and nw if nw>1)
    ## Xprep and para are calculated and used for calculating lods.
    ## If Xprep is given in the function input, M, nb and nw are unnecessary.

    ## Setting up Data
    if(is.null(Xprep))
      {
    	Xprep<- setup.bayesian(M=M, nb=nb, nw=nw)
      }
    nb<-Xprep$nb
    nw<-Xprep$nw
    Mbar<-Xprep$Mbar
    SSW<-Xprep$SSW
    SSB<-Xprep$SSB

    ## Setting up parameters
     if(is.null(para$v) | is.null(para$a)) 
     {
	va <- va.func(Vest=SSB/(nb-1), k=nb)
	if (is.null(para$v)) para$v<-va$v
	if (is.null(para$a)) para$a<-va$a
     }

     if(is.null(para$k)) 
     {
	if(is.null(SSW)) para$k<-0  else para$k<-median(SSB/(nb-1)/SSW*nb*(nw-1))
     }

     if(is.null(para$c)) para$c<-c.func(Xprep=Xprep, para=para)

     lods<-lods.func(Xprep, para)
     list(Xprep=Xprep, lods=lods,para=para)
}


stat.bay.est <- function(M=NULL, Xprep=NULL, para=list(p=0.01, v = NULL, a=NULL, c = NULL))
  {
    ## Input Xprep or M (output from stat.ma).
    ## If M is given, stat.bay.est assumes the experiment consists of ncol(M) microarray slides all measuring the same effect (which will be stimated by Mbar)
    ## Para is calculated and used for calculating lods.
    
    ## Setting up Data
    if(is.null(Xprep))
      {
    	Xprep$Mbar<-apply(M,1,mean.na)   #Effect estimate
        Xprep$Vest<-apply(M,1,var.na)    #Variance estimate
        Xprep$k<-ncol(M)                 #Variance constant
        Xprep$f<-Xprep$k-1               #Degrees of freedom
      }
    Mbar<-Xprep$Mbar
    Vest<-Xprep$Vest
    k<-Xprep$k
    f<-Xprep$f

    ## Setting up parameters
     if(is.null(para$v) | is.null(para$a)) 
     {
	va <- va.func(Vest=Vest, k=k)
	if (is.null(para$v)) para$v<-va$v
	if (is.null(para$a)) para$a<-va$a
     }

     if(is.null(para$c)) para$c<-c.func.est(Xprep=Xprep, para=para)

     lods<-lods.func.est(Xprep, para)
     list(Xprep=Xprep, lods=lods,para=para)
}


plot.bayesian<-function(x=NULL,Mbar=x$Xprep$Mbar,lods=x$lods, type="t",spec=50, ch=NULL, col='black',...){
  #bay <- x
  index<-NULL
  if(type=="t") index<-(1:length(lods))[lods>sort(lods[!is.na(lods)])[length(lods[!is.na(lods)])-spec]]   
  if(type=="c") index<-(1:length(lods))[lods >= spec]
  if(type=="i") index<-spec                       
  plot(Mbar, lods, xlab="Effect estimate", ylab="Lodsratio", main="Lodsratio vs Effect estimate", type="n")
  if (!is.null(index)){
  points(Mbar[-index], lods[-index], pch='.')    
  if(is.null(ch))  text(Mbar[index],lods[index],labels=index,col=col)  else points(Mbar[index],lods[index],pch=ch,col=col)
}
  if(is.null(index)){points(Mbar, lods, pch='.')}
}
  
    
######################################################
#Functions
######################################################

setup.bayesian <- function(M, nb=NULL, nw=1)
  {
    if (nw == 1)
    {
	nb<-ncol(M)
	SSW<-NULL
	SSB<-apply(M,1,var.na)*(nb-1)
	Mbar<-apply(M,1,mean.na)
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
	      Mtmp<-cbind(Mtmp,M[seq(j,nrow(M),nw),i])
	    }
	  }
	  Mbar<-apply(Mtmp,1,mean.na)

	  Mslide<-NULL
	  SSW<-rep(0,nrow(M)/nw)
	  SSB<-rep(0,nrow(M)/nw)
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
	    Mtmp<-cbind(Mtmp,M[seq(j,length(M),nw)])
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

####################################################################
#Lods (or lods.func.est) calculates the logodds ratio.

lods.func<-function(Xprep=list(Mbar=Mbar, SSB=SSB, SSW=SSW, nb=nb, nw=nw), para=list(p=p, c=c, v=v, a=a, k=k)){
   Mbar <- Xprep$Mbar            #overall means
   nb<-Xprep$nb
   nw<-Xprep$nw
   SSB <- Xprep$SSB              #sums of squares between slides
   SSW <- Xprep$SSW              #sums of squares within slides  
   if(is.null(Xprep$SSW)) SSW <- rep(0, length(SSB))

   p <- para$p
   v <- para$v
   a <- para$a
   c <- para$c
   k <- para$k

   odds1<-p/(1-p)
   odds2<-1/(1+nb*nw*c)
   odds3<-a + (SSB + k*SSW)/nw/nb   
   odds4<-(odds3+Mbar^2)/(odds3+Mbar^2/(1+nb*nw*c))
   odds<-odds1*(odds2**(1/2))*(odds4**(v+nw*nb/2))
   log(odds)
 }

lods.func.est<-function(Xprep=list(Mbar=Mbar, Vest=Vest, k=k, f=f), para=list(p=p, c=c, v=v, a=a)){

   Mbar<-Xprep$Mbar
   Vest<-Xprep$Vest
   k<-Xprep$k
   f<-Xprep$f

   p <- para$p
   v <- para$v
   a <- para$a
   c <- para$c

   odds1<-p/(1-p)
   odds2<-1/(1+k*c)
   odds3<-a + f*Vest/k  
   odds4<-(odds3+Mbar^2)/(odds3+Mbar^2/(1+k*c))
   odds<-odds1*(odds2**(1/2))*(odds4**(v+f/2+1/2))
   log(odds)
 }

######################################################

#va.func estimates v and a by the method of moments, so that a*k/(2*sigma^2) ~Gamma(v,1), Vest are the genewise estimates of sigma^2 and k a constant such that the expected variance of the effect estimates are sigma^2/k. 

va.func<-function(Vest, k){
av.var<-mean.na(Vest)
var.var<-sum.na((Vest-mean.na(Vest))^2)/(length.na(Vest)-1)
vhat<-(2*var.var+av.var^2)/var.var
ahat<-av.var/k*2*(vhat-1)
list (v=vhat, a=ahat)
}


######################################################
#c.func (or c.func.est) uses ls.variance and sq.func to estimate c. It uses a least squares estimate so that for all Mbar-values, Mbar~N(0,simga^2/k) and for the top p proportion of the genes, the averages are ~N(0,c*sigma^2).

c.func<-function(Xprep, para){

#Estimate the variance sigma²
  var.est<-mean.na(Xprep$SSB/(Xprep$nb-1))
  start<-var.est/10/Xprep$nb
  end<-var.est*10/Xprep$nb
  
  sigma2<-ls.variance(X=Xprep$Mbar[!is.na(Xprep$Mbar)], var.start=start, var.stop=end, nclass=100)*Xprep$nb
  
  sigma2


#Estimate c*sigma²
  para$c<-1.5
  if (is.null(Xprep$SSW)) Xprep$SSW<-rep(0,length(Xprep$SSB))
  l<-stat.bayesian(Xprep=Xprep, para=para)$lods
  top.set<-(1:(length(Xprep$Mbar)/Xprep$nw))[!is.na(l)][rank(l[!is.na(l)])>(length.na(l)-round(length.na(l)*para$p))]

  var.est<-var.na(Xprep$Mbar[top.set])
  start<-var.est/10
  end<-var.est*10
  
  csigma2<-ls.variance(X=Xprep$Mbar[top.set], var.start=start, var.stop=end, zeros=FALSE)
  
  csigma2/sigma2
}


c.func.est<-function(Xprep, para){

#Estimate the variance sigma²
  var.est<-mean.na(Xprep$Vest)
  start<-var.est/10/Xprep$k
  end<-var.est*10/Xprep$k
  
  sigma2<-ls.variance(X=Xprep$Mbar[!is.na(Xprep$Mbar)], var.start=start, var.stop=end, nclass=100)*Xprep$k
  
  sigma2


#Estimate c*sigma²
  para$c<-1.5
  l<-stat.bay.est(Xprep=Xprep, para=para)$lods
  top.set<-(1:(length(Xprep$Mbar)))[!is.na(l)][rank(l[!is.na(l)])>(length.na(l)-round(length.na(l)*para$p))]

  var.est<-mean.na(Xprep$Vest[top.set])
  start<-var.est/10
  end<-var.est*10
  
  csigma2<-ls.variance(X=Xprep$Mbar[top.set], var.start=start, var.stop=end, zeros=FALSE)
  
  csigma2/sigma2
}


ls.variance<-function(X, var.start, var.stop, var.steps=100, nclass=NULL, zeros=TRUE) {

  var.seq<-seq(var.start, var.stop, (var.stop-var.start)/var.steps)
  sq<-rep(0,length(var.seq))
  i<-0
  for (v in var.seq){
    i<-i+1
      sq[i] <- sq.func(X=X, var=v, nclass=nclass, zeros=zeros)
  }  
  
  minsq<-min(sq)
  for (i in 1:length(var.seq)){
      if (sq[i]==minsq)  min.v<-i
  }
  var.seq[min.v]
}

sq.func<-function(X,var,nclass=NULL, zeros=TRUE){

    if (is.null(nclass)) index<-1:length(X) else index <- seq(1, length(X), round(length(X)/nclass))
    hst<-hist(X[index], plot=FALSE, freq=FALSE, nclass=20)

    if (zeros) sumindex<-(1:length(hst$counts)) else sumindex<-(1:length(hst$counts))[hst$density > 0]

    theo<-dnorm(hst$mids[sumindex], mean=0 ,sd=sqrt(var))
    sum.na((hst$density[sumindex]-theo)**2)
}












