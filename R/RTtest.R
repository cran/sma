###########################################################################
# Statistics for Microarray Analysis
# T-test
#
# Date : March 19, 2001
#
# History:
#    March 19, 2001: Some of the plot functions from Rarray.R.
#					
#
# Authors: Sandrine Dudoit and Yee Hwa (Jean) Yang.
##########################################################################

##########################################################################
# Test statistics for multiple slides only
##########################################################################

########################################################################/**
# \name{stat.t2}
# 
# \alias{stat.t2}
# \alias{t2stat.func}
# 
# \title{Two-sample t-statistics}
# 
# \description{
# Computes two-sample t-statistics for each gene in a multi-slide
# microarray experiment. 
# }
# 
# \usage{
# stat.t2(X, cl, x.ratio=FALSE, var.equal=TRUE, ...)
# }
# 
# \arguments{
#  \item{X}{if x.ratio=F, X is a list containing two components.  The
#  first component is a matrix of log intensity ratios
#  \eqn{M=\log_2 (R/G)} and the second component is the 
#   average log intensities \eqn{A = log_2 \sqrt{RG}}{A =
#   log_2(R*G)/2}, such as the output 
#   from \code{\link{stat.ma}}. If x.ratio=T, X is a matrix of log
#   expression ratios only.   The rows of X correspond to genes and 
#   columns correspond to different hybridizations, that is different
#   slides.}  
# 
#  \item{cl}{vector of class labels. Must consist of integers 1 and 2.}
# 
#  \item{x.ratio}{logical flag: if TRUE, the matrix X contains only
#  log intensity ratios, if FALSE, X is a list containing two
#  components.  The first component is a matrix of log expression
#  ratios and the second component contains average log
#  intensities A.}
# 
#  \item{var.equal}{logical flag: if TRUE, the variances of the class
#  1 and class 2 parent populations are assumed equal.} 
# 
#  \item{\dots}{other parameters used in \code{\link{t.test}}. }
# }
# 
# \value{
# List containing the following components
# 
#   \item{t}{the two-sample t-statistic for each gene;}
# 
#   \item{Num }{the numerator of the t-statistic for each gene, with
# names attribute "Num";}
# 
#   \item{Denominator}{the denominator of the t-statistic for each gene, with
# names attribute "Den";}
# 
#   \item{n1}{number of class 1 observations used to calculate the
#   t-statistic for each gene;}
# 
#   \item{n2}{number of class 2 observations used to calculate the
#     t-statistics for each gene;}   
# 
#     \item{Average A}{if x.ratio=F, the average across all
#     hybridizations of \eqn{A = log_2 \sqrt{RG}}{A = log_2(R*G)/2}, 
# with names attribute "A.bar", if x.ratio=T, NULL is returned.}
# }
# 
# \references{ D. Freedman, R. Pisani, and
# R. Purves. (1998). Statistics, 3rd ed. NewYork: W.W. Norton.} 
# 
#  
# \note{\code{\link{t2stat.func}} is called by \code{\link{stat.t2}}
# and is not typically used on its own.}         
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} 
# }
#   
# \seealso{\code{\link{t2stat.func}}, \code{\link{plot.t2}},
# \code{\link{plot.qq}}, \code{\link{t.test}}.} 
# 
# \examples{
# data(MouseArray)
# ## mouse.setup <- init.grid()
# ## mouse.data <- init.data() ## see \emph{init.data}
# ## mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# cl <- c(rep(1,3), rep(2,3))
# mouse.t2 <- stat.t2(mouse.lratio, cl)
# }
# 
# \keyword{T-test.}
#*/#########################################################################

stat.t2<-function(X, cl, x.ratio=FALSE, var.equal=TRUE,  ...)
{ 
  if(!x.ratio){
    n <- ncol(X)/2
    res<-t(apply(X$M,1,function(z) t2stat.func(z,cl,var.equal, ...)))
    A.bar<-apply(X$A,1,mean.na)
  }
  if(x.ratio){
    n <- ncol(X)
    res<-t(apply(X,1,function(z) t2stat.func(z,cl,var.equal,  ...)))
    A.bar<-NULL
  }
  list(t=res[,"t"],Num=res[,"Num"],Den=res[,"Den"],n1=res[,"n1"],n2=res[,"n2"], A.bar = A.bar)
}


##########################################################################
# Internal Function called by stat.t2
##########################################################################

t2stat.func<-function(x,cl,var.equal=TRUE, ...)
{
  x.ok<-x[!(is.na(x) | is.infinite(x))]
  cl.ok<-cl[!(is.na(x)| is.infinite(x))]
  
  x1<-x.ok[cl.ok==1]
  x2<-x.ok[cl.ok==2]
  n1<-length(x1)
  n2<-length(x2)

  if((n1>2) & (n2>2))
    {
      tmp<-t.test(x1, x2,  var.equal=var.equal, ...)
      tstat<--(tmp$stat)
      num<-tmp$est[2]-tmp$est[1]
      den<-num/tstat
      res<-c(tstat,num,den,n1,n2)
    }
  if((n1<=2) | (n2 <=2))
      res<-c(NA,NA,NA,n1,n2)
  names(res) <- c("t", "Num", "Den", "n1", "n2")
  res
}

############################################################################
#                              End of file
############################################################################
