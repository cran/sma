###########################################################################
# Statistics for Microarray Analysis
# Misc plots
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
# Exploratory plots for single slide
##########################################################################

########################################################################/**
# \name{plot.svb}
# 
# \alias{plot.svb}
# \alias{svb.func}
# 
# \title{Plot of Signal vs. Background}
# 
# \description{
# Produces a scatter plot of background corrected signal intensities
# and background intensities.
# }
# 
# \usage{
# plot.svb(x, channel="R", image.id=1, S.isbgcor=F, ...)
# }
# 
# \arguments{
#   \item{x}{a numeric list of signal and background intensities, can
#   be raw or background corrected data.}   
# 
#  \item{channel}{the specific channel to which the intensities to be
#    considered, correspond to, that is, either red or green. The default
#    channel is red.}
# 
#  \item{image.id}{integer value; the index of the slide which is considered}
# 
#  \item{S.isbgcor}{logical flag, equal to TRUE if the signal intensities in
#    x contain background corrected signal intensities instead of raw
#    signal intensities. By default this is set to FALSE.}
# 
#  \item{\dots}{graphical parameters may also be supplied as arguments
#  to the function (see \code{\link{par}}).} 
# 
# }
# 
# \value{a plot is created on the current graphics device.}
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Jessica Mar}
# }
#   \seealso{\code{\link{plot}}.} 
# 
# \examples{
# data(MouseArray)
# # mouse.setup <- init.grid()
# # mouse.data <- init.data() ## see \emph{init.data}
# 
# plot.svb(mouse.data, "green", 3) 
# ## thiscreates a plot of the signal versus background intensities 
# ## for the green channel, using data collected from the third slide. 
# }
# 
# \keyword{microarray, background}
# 
#*/######################################################################/**
plot.svb <- function(x, channel="R", image.id=1, S.isbgcor=FALSE, ...)
   {if(is.list(x))
     if ((channel=="R") | (channel=="r") | (channel=="red"))
     {svb.func(x$R[, image.id], x$Rb[, image.id], S.isbgcor=FALSE, ...)}
       
     if ((channel=="G") | (channel=="g") | (channel=="green"))
     {svb.func(x$G[, image.id], x$Gb[, image.id], S.isbgcor=FALSE, ...)}
   }

svb.func <- function(Signal, Bg, S.isbgcor = FALSE, ...)
{
  if(S.isbgcor){
    ind <- log.na(Signal, 2) < quantile(log.na(Signal, 2), 0.75,  na.rm=TRUE)
    plot(log.na(Bg,2)[ind], log.na(Signal,2)[ind], 
	 xlab="Background", ylab="Signal", ...)
  }
  if(!S.isbgcor){
    ind <- log.na(Signal-Bg, 2) < quantile(log.na(Signal-Bg, 2), 0.75,  na.rm=TRUE)
    plot(log.na(Bg,2)[ind], log.na(Signal-Bg,2)[ind], 
	 xlab="Background", ylab="Signal", ...)
  }
}

#######################################################
# plot.print.tip.lowess - a function to print individual
# lowess curves for each print tip super imposed onto an
# m vs a  plot.
#
# TO be added a linetype palette, functionality to add
# an index
#
#######################################################

plot.print.tip.lowess <- function (x, layout, norm = "n", image.id = 1, palette = rainbow(layout$ngrid.r*layout$ngrid.c), lty.palette = rep(1,layout$ngrid.r*layout$ngrid.c), ...) 
{
    tmp <- ma.func(R = x$R[, image.id], G = x$G[, image.id], 
        Rb = x$Rb[, image.id], Gb = x$Gb[, image.id], layout, 
        norm = norm, pout = FALSE, ...)
    plot(tmp$A, tmp$M, xlab = "A", ylab = "M", ...)
    npin <- layout$ngrid.r * layout$ngrid.c
    nspot <- layout$nspot.c * layout$nspot.r
    pin <- rep(1:npin, rep(nspot, npin))
    for (i in 1:npin) {
        index <- pin == i
        tM <- tmp$M[index]
        tA <- tmp$A[index]
        ind2 <- is.na(tM) | is.na(tA) | is.infinite(tM) | is.infinite(tA)
        smoothnum <- lowess(tA[!ind2], tM[!ind2])
        lines(approx(smoothnum), col = palette[i], lty = lty.palette[i], ...)
    }
}





##########################################################################
# Diagnostic plots for multiple slides only
##########################################################################

########################################################################/**
# \name{plot.qq}
# 
# \alias{plot.qq}
# 
# \title{ Histogram and Normal Quantile-Quantile plot}
# 
# \description{Produces a histogram and a normal Quantile-Quantile plot
# of the data. The points corresponding to genes with statistics
# less/greater than a user defined threshold are highlighted. The
# histogram and Q-Q plots are displayed on the same page. 
# }
# 
# \usage{
# plot.qq(x, name, low=-5, high=5)
# }
# 
# \arguments{
# 
#  \item{x}{a numeric vector containing the statistics whose histogram
#  and Q-Q plot will be produced. Missing values (NAs) are allowed.}
# 
#  \item{name}{title for the plots.}
# 
#  \item{low}{lower threshold: points with statistic < low are colored
#  in green.}
# 
#  \item{high}{upper threshold: points with statistic > high are
#  colored in red.} 
# } 
# 
# \references{Chambers, J. M., Cleveland, W. S., Kleiner, B. and
# Tukey, P. A. (1983). Graphical Methods for Data Analysis. Wadsworth,
# Belmont, California. 
#  
# Hoaglin, D. C., Mosteller, F. and Tukey, J.  W., editors
#  (1983). Understanding Robust and Exploratory Data Analysis. Wiley,
#  New York.         
# }
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} 
# }
# 
# \seealso{\code{\link{plot.spatial}}, \code{\link{plot.t2}},
# \code{\link{stat.t2}}, \code{\link{hist}}, \code{\link{qqnorm}}.} 
# 
# \examples{
# data(MouseArray)
# ## mouse.setup <- init.grid()
# ## mouse.data <- init.data() ## see \emph{init.data}
# ## mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# ## Calculation of t-statistics
# ## cl <- c(rep(1,3), rep(2,3))
# ## mouse.t2 <- stat.t2(mouse.lratio, cl)
# 
# ## Diagnostic plots
# plot.qq(mouse.t2$t, "Mouse")
# }     
# 
# \keyword{microarray, histogram, qqplot.}
# 
#*/######################################################################/**

plot.qq <- function(x,name, low=-5, high=5,...)
{
  par(mfrow=c(2,1))
  hist(x,xlab="t",nclass=100,main=paste(name," Histogram and quantile-quantile plot of t-statistics", sep=":"),col=9,cex=0.8)
  tmp<-qqnorm(x,plot=FALSE)
  plot(tmp,pch=".",xlab="Quantiles of standard normal",ylab="t")
  points(tmp$x[tmp$y<low],tmp$y[tmp$y<low],pch="*",col=6,cex=1.2)
  points(tmp$x[tmp$y>high],tmp$y[tmp$y>high],pch="*",col=2,cex=1.2)
}

 
#########################################################################/**
# \name{plot.qqline}
# 
# \alias{plot.qqline}
# 
# \title{Add Line Going Through the Quantiles of a Q-Q Plot}
# 
# \description{
# This function adds a line to a quantile-quantile plot which passes
# through user defined quantiles. This function is similar to, but
# more general than, \code{\link{qqline}} because the reference
# distribution need not be the standard normal distribution and the
# quantiles need not be the first and third quartiles. \cr 
# Graphical parameters may be given as arguments to \code{plot.qqline}. 
# }
# 
# \usage{
# plot.qqline(x, y, a=0.25, ...)
# }
# 
# \arguments{
# \item{x}{the reference (first) sample for the Q-Q plot, for a normal Q-Q
#   plot this would be the quantiles of a N(0,1) random sample.}
# \item{y}{the data.}
# }
# \item{a}{a number between 0 and 1. A line is drawn which connects the
#   \code{a} and \code{1-a} quantile points. The default line passes
#   through the first and third quantiles.}
#  \item{\dots}{graphical parameters may also be supplied as arguments
#  to the function (see \code{\link{par}}).} 
# }
# 
# \seealso{\code{\link{qqplot}}, \code{\link{qqnorm}}, \code{\link{qqline}}.  }
# 
# \examples{
# data(MouseArray)
# # mouse.setup <- init.grid()
# # mouse.data <- init.data() ## see \emph{init.data}
# # mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# ## Calculation of t-statistics
# ## cl <- c(rep(1,3), rep(2,3))
# ## mouse.t2 <- stat.t2(mouse.lratio, cl)
# 
# ## Diagnostic plots
# plot.qq(mouse.t2$t, "Mouse")
# 
# ## Using the QQline function
# q <- quantile(rnorm(1000))
# plot.qqline(q, mouse.t2$t)
# }    
# 
# \keyword{Q-Q plots, quartiles}
# 
#*/#########################################################################

plot.qqline<-function(x,y,a=0.25, ...)
{
    y <- quantile(y[!is.na(y)],c(a, 1-a))
    x <- quantile(x[!is.na(x)],c(a, 1-a))
    points(x,y,...)
    slope <- diff(y)/diff(x)
    int <- y[1]-slope*x[1]
    abline(int, slope, ...)
}                       

#########################################################################/**
# \name{plot.scale.box}
# \alias{plot.scale.box}
# 
# \title{Box plots for microarray}
# \description{
# Produce box-and-whisker plot(s) of the given (grouped) values.
# }
# \usage{
# plot.scale.box(x, layout, x.names=NULL, ...)
# }
# 
# \arguments{
#   \item{x}{a vector or a matrix.}
#   \item{layout}{ a list specifying the dimensions of the spot matrix
#    and the grid matrix.  This can be generated by calling
#    \code{\link{init.grid}}.}
#   \item{x.names}{group labels which will be printed under each boxplot.}
#   \item{\dots}{further arguments to the default boxplot method and graphical
#     parameters may also be passed as arguments, see \code{\link{par}}.}
# }
# \details{
#   If x is a vector, this function will produce n boxplots where n is
#   number of print-tips groups.   If x is a matrix, this function will
#   produce n boxplots where n is number of columns in the matrix.  
# }
# 
# \author{Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu}}
# 
# \seealso{\code{\link{boxplot}}, \code{\link{bxp}}}
# 
# \examples{
#      data(MouseArray)
#      # mouse.setup <- init.grid() 
#      # mouse.data <- init.data() ## see \emph{init.data} 
#      mouse.lratio <- stat.ma(mouse.data, mouse.setup)
#      ## Producing boxplots for different print-tips groups.
#      plot.scale.box(mouse.lratio$M[,1], mouse.setup)
#      ## Producing boxplots for different slides.
#      plot.scale.box(mouse.lratio$M)
# }
# \keyword{boxplots, microarray}
# 
#*/#########################################################################


plot.scale.box <- 
function(x, layout=NULL, x.names=NULL, ...)
  {
    n <- layout$nspot.r * layout$nspot.c * layout$ngrid.r * layout$ngrid.c
    nperpin <- layout$nspot.r * layout$nspot.c
    npin <- layout$ngrid.r * layout$ngrid.c
    
    if(is.vector(x)){
      if((length(x) != n) & (is.null(x.names)))
        {
          stop(" Error: Length of vector different from total number of spots and vector has no row.name.\n")
        }
      if ((length(as.vector(x)) != n) & (!is.null(x.names)))
        {
          y <- x; x <- rep(NA, n);
          x[as.integer(x.names)] <- y
      }
      xmat <- matrix(x, nrow = nperpin)
       vect <- TRUE
    }

    if(is.matrix(x))
      xmat <- x
    
    boxplot(data.frame(xmat), names=x.names, ...)
  }



########################################################################/**
# \name{plot.t2}
# 
# \alias{plot.t2}
# 
# \title{Diagnostic Plots for Two-Sample t-statistics}
# 
# \description{
# Plots of two-sample t-statistics, |t-numerator| and
# t-denominator against average A, and plot of |t-numerator| against
# t-denominator. For each spot on a given slide, \eqn{A = log_2
# \sqrt{RG}}{A = log_2(R*G)/2}, where (R,G) denotes the red and green
# fluorescence intensity pair. Points with t-statistics exceeding user
# defined thresholds are highlighted. 
# } 
# 
# \usage{
# plot.t2(X, main.title="T plots", low=-5, high=5)
# }
# 
# \arguments{
# 
#  \item{X}{output from the function \code{\link{stat.t2}}.}
# 
#  \item{main.title}{title for the plot.}
# 
#  \item{low}{lower threshold for t-statistic: points with t<low are
#  colored in green.} 
# 
#  \item{high}{upper threshold for t-statistic: points with t>high are
#  colored in red.} 
# }
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} 
# }
# \seealso{\code{\link{stat.t2}}, \code{\link{t2stat.func}},
# \code{\link{plot}}, \code{\link{t.test}}.} 
# 
# \examples{
# data(MouseArray)
# # mouse.setup <- init.grid()
# # mouse.data <- init.data() ## see \emph{init.data}
# # mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# ## Calculation of t-statistics
# ## cl <- c(rep(1,3), rep(2,3))
# ## mouse.t2 <- stat.t2(mouse.lratio, cl)
# 
# ## Diagnostic plots
# plot.t2(mouse.t2, "Mouse")
# }    
# 
# \keyword{microarray, ttest.}
# 
#*/######################################################################/**


plot.t2 <- function(x, main.title="T plots", low=-5,high=5,...)
{
  par(mfrow=c(2,2),oma=c(1,1,3,1))

  lowt<-x$t < low
  hight<-x$t > high
 
  # 1. t vs. avg. A
  plot(x$A.bar,x$t,xlab="average A",ylab="t",pch=".",
       main="t vs. average A",cex=0.8)
  if(sum.na(lowt)>0)
    points(x$A.bar[lowt],x$t[lowt],pch="*",col=6,cex=1.5)
  if(sum.na(hight)>0)
    points(x$A.bar[hight],x$t[hight],pch="*",col=2,cex=1.5)

   # 2. t_denom vs. avg. A
   plot(x$A.bar,x$Den,xlab="average A",ylab="t denominator",
	pch=".",main="t denominator vs. average A",cex=0.8)
  if(sum.na(lowt)>0)
    points(x$A.bar[lowt],x$Den[lowt],pch="*",col=6,cex=1.5)
  if(sum.na(hight)>0)
    points(x$A.bar[hight],x$Den[hight],pch="*",col=2,cex=1.5)

  # 3. |t_num| vs. avg. A
  plot(x$A.bar,abs(x$Num),xlab="average A",ylab="|t numerator|",
       pch=".",main="|t numerator| vs. average A",cex=0.8)
  if(sum.na(lowt)>0)
    points(x$A.bar[lowt],abs(x$Num[lowt]),pch="*",col=6,cex=1.5)
  if(sum.na(hight)>0)
    points(x$A.bar[hight],abs(x$Num[hight]),pch="*",col=2,cex=1.5)

  # 4. t_denom vs. |t_num|
  plot(abs(x$Num) , x$Den ,xlab="|t numerator|",ylab="t denominator",
       pch=".",main="t denominator vs. |t numerator|",cex=0.8)
  if(sum.na(lowt)>0)
    points(abs(x$Num[lowt]),x$Den[lowt],pch="*",col=6,cex=1.5)
  if(sum.na(hight)>0)
    points(abs(x$Num[hight]),x$Den[hight],pch="*",col=2,cex=1.5)

  par(mfrow=c(1,1))
  mtext(main.title, line=4,cex=1.5)
}

