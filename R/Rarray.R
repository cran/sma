###########################################################################
# Statistics for Microarray Analysis
# Exploratory analysis - Mainly preprocessing.
#
# Date : August 9, 2000
# Last update : May 17, 2001
#
# History:
#   May 17, 2001: Fix to norm.scale.func
#   March, 19: Splitting Rarray in to smaller files.  
#              Including Comments at the start of each function.
#   Nov, 20: Change the argument on plot.mva...it's not usable otherwise.
#            Bug fix ma.func
#   Nov, 13: Ben's Bug fix on stat.ma
#   Nov, 10: Change data structure from matrix to list of matrix.  
#   Sept, 28: Bug fix: ma.func
#   Feb 20, 2003 - bug fix to ma.func (As suggested by G. Smyth)
#   Apr 27, 2003 - fix bug in ma.func (when both R-Rb and G-Gb are negative should M give NA)
#
# Authors: Sandrine Dudoit and Yee Hwa (Jean) Yang.
##########################################################################


##########################################################################
#  stat.gnames
#  History:  
#     March 19, 2001:  remove infinite values from the ordering.
#
##########################################################################

#########################################################################/**
# \name{stat.gnames}
# 
# \alias{stat.gnames}
# 
# \title{Sort Genes According to the Value of a Statistic}
# 
# \description{
# Lists genes and corresponding statistics in decreasing order of the
# statistics. This function applies to any type of statistic, including
# log ratios, one and two-sample t-statistics, and F-statistics. Missing
# values are ignored, as in \code{\link{sort}(..., na.last=NA)}. 
# }
# 
# \usage{
# stat.gnames(x, gnames, crit=0.05)
# }
# 
# \arguments{
#  \item{x}{a numeric vector containing the statistics for each
#  gene. Missing values (NAs) are allowed. }
#  
# \item{gnames}{a character vector containing the gene names.}
# 
#  \item{crit}{specifies the number of genes to be returned. If crit <
#  1, the crit*100\% genes with the largest x values are listed. If crit
#  >= 1, the crit genes with the largest x values are listed. }
# }
# 
# \value{
# List containing the following components 
#   \item{gnames}{gene names sorted in decreasing order of the
#  statistics in x.}
#  \item{t}{statistics sorted in decreasing order.}
# }
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} }
# 
# \seealso{\code{\link{stat.t2}}, \code{\link{order}}, \code{\link{sort}}.}
# 
# \examples{
# ## Calculating log ratio and performing a t test.
# data(MouseArray)
# ## mouse.setup <- init.grid()
# ## mouse.data <- init.data() ## see \emph{init.data}
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# cl <- c(rep(1,3), rep(2,3))
# mouse.t2 <- stat.t2(mouse.lratio, cl)
# 
# ## Looking at gene names
# ## Finding the top 10 t-statistics
# stat.gnames(abs(mouse.t2$t), mouse.gnames, crit=10)
# 
# ## Finding the top 1% of t-statistics
# stat.gnames(abs(mouse.t2$t), mouse.gnames, crit=0.01)
# 
# ## Finding the 10 extreme M values in the first slide
# stat.gnames(abs(mouse.lratio$M[, 1]), mouse.gnames, crit=10)
# }
# 
# \keyword{microarray.}
#*/#########################################################################

stat.gnames<-function(x, gnames, crit=0.05)
{
    ind <- is.infinite(x)
    x <- x[!ind]
    if (crit < 1) {
        which <- rev(order.na(x, na.last = FALSE))[1:(round(length(x) * 
            crit))]
        if (sum(is.na(x)) > (length(x) - round(length(x) * crit))) 
            warning("NA exists under your selection criteria")
    }
    if (crit >= 1) {
        which <- rev(order.na(x, na.last = FALSE))[1:crit]
        if (sum(is.na(x)) > (length(x) - crit)) 
            warning("NA exists under your selection criteria")
    }
    if (is.matrix(gnames) | is.data.frame(gnames)) 
      {
	gnames <- gnames[!ind, ]
        res <- list(gnames = gnames[which, ], t = x[which])
      }
    if (is.vector(gnames)) 
      {
	gnames <- gnames[!ind]
        res <- list(gnames = gnames[which], t = x[which])
      }
    res
}


##########################################################################
# Calculation of log ratios and normalization
# History:
#    March, 19, 2001: Modify the function so that stat.ma works with data 
#                     with no background values.
##########################################################################

#########################################################################/**
# \name{stat.ma}
# 
# \alias{stat.ma}
# \alias{ma.func}
# \alias{norm.l.func}
# \alias{norm.pin.func}
# \alias{norm.scale.func}
# 
# \title{Calculation of log Intensity Ratios and Average log Intensities}
# 
# \description{
# Computes the log intensity ratio \eqn{M = log_2 (R/G)} and the mean log
# intensity \eqn{A = log_2 \sqrt{RG}}{A = log_2(R*G)/2}, where R and G
# represent the fluorescence
# intensities in the red and green channels, respectively. Logarithms base
# 2 are used instead of natural or decimal logarithms as intensities are
# typically integers between 1 and \eqn{2^{16}}. The log intensity
# ratios M are normalized using one of the five available methods. 
# }
# 
# \usage{
# stat.ma(RG, layout, norm="p", pout=FALSE, ...)
# }
# 
# \arguments{
#   \item{RG}{
#     a list with 4 elements, each represents a matrix with p rows for p
#     genes and n columns for n slides. \cr
#     The first element "R" contains the raw red intensities from slide
#     i=1,...,n .\cr
#     Similarly, the second element "G" contains the raw green
#     intensities. \cr
#     The third element "Rb"  contains the background red intensities and \cr
#     the fourth element "Gb" contains the  background green intensities.\cr
#     This list structure can be generated by the interactive function
#     \code{\link{init.data}}. }
#   
#   \item{layout}{a list specifying the dimensions of the spot matrix
#   and the grid  
#     matrix.  This can be generated by calling \code{\link{init.grid}}.}
# 
#   \item{norm}{Character string, one of "n", "m", "l", "p" or "s".  This
#     argument specifies the type of normalization method to be performed:
#     "n" no normalization between the 2 channels; "m"
#     \code{\link{median}} normalization, which sets the median of log
#     intensity ratios to zero; "l" global \code{\link{lowess}}
#     normalization; "p" print-tip group lowess normalization and "s"
#     scaled print-tip group lowess normalization. The default method is
#     set to print-tip normalization.}
#
#   \item{pout}{if TRUE, an M vs. A plot will be produced. Otherwise,
#   a matrix of log intensity ratios and average log intensities is
#   return.  By default pout is set to FALSE.  The option pout='TURE'
#   is not yet implemented.}
#
#   \item{\dots}{other parameters used in \code{\link{ma.func}}. }
# }
# 
# \value{
#   List containing the following components:
#   
#   \item{M}{Matrix of log expression ratios \eqn{M = log_2 (R/G)}}
#   \item{A}{Matrix of average log intensities \eqn{A = log_2
#       \sqrt{RG}}{A = log_2(R*G)/2}}
#   For the matrix in each of the components, rows correspond to genes
#   and columns correspond to different hybridizations, that is
#   different slides.  
# }
# 
# \references{S. Dudoit, Y. H. Yang, M. J. Callow, and T. P. Speed. Statistical
# methods for identifying differentially expressed genes in replicated
# cDNA microarray experiments (Statistics, UC Berkeley, Tech Report \# 578).  }
# 
# \note{\code{\link{ma.func}}, \code{\link{norm.l.func}},
#   \code{\link{norm.scale.func}} and \code{\link{norm.pin.func}} are called by \code{\link{stat.ma}} and are not typically used on their own.}
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Natalie Roberts, \email{nroberts@wehi.edu.au}
# }
# 
# \seealso{\code{\link{ma.func}}, \code{\link{norm.l.func}},
#   \code{\link{norm.pin.func}}, \code{\link{norm.scale.func}}, \code{\link{plot.mva}}, \code{\link{lowess}}.}
# 
# \examples{
# data(MouseArray)
# ## mouse.setup <- init.grid() 
# ## mouse.data <- init.data() ## see \emph{init.data} 
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# }
# 
# \keyword{microarray, log ratio.}
#*/#########################################################################

stat.ma <- function(RG, layout, norm="p", pout=FALSE, ...)
{
  n <- ncol(RG$R)
  res <- list(A=NULL, M=NULL)

  for(i in (1:n))
    {
##       RG <- apply(RG, 2, as.numeric)
      if(is.null(RG$R[,i])){
	stop(" Error: No data is given in RG$R \n")
      }
      if(is.null(RG$G[,i])){
	stop(" Error: No data is given in RG$G\n")
      }
      if(pout)
        stop("pout=TRUE is not implemented")
      tmp <-ma.func(R=RG$R[,i],G=RG$G[,i],Rb=RG$Rb[,i], Gb=RG$Gb[,i], layout=layout, norm=norm, pout=pout, ...)
      res$A<-cbind(res$A, tmp$A)
      res$M<-cbind(res$M, tmp$M)
    }
  res
}

##########################################################################
#  stat.norm.exp
#  History:  
#     March 19, 2001:  Original
#
##########################################################################

#########################################################################/**
# \name{stat.norm.exp}
# \alias{stat.norm.exp}
# \title{Normalization of log Intensity Ratios across slides / experiments.}
# 
# \description{
# Performs scale normalization across slides (experiments)}
# }
# 
# \usage{
# stat.norm.exp(X)
# }
# 
# \arguments{
#   \item{X}{X is a matrix of log intensity ratios \eqn{M=\log_2 (R/G)}
#   The rows of X correspond to genes and columns correspond to different 
#   hybridizations, that is different slides (experiments). 
# }
# 
# \value{
#   A matrix of normalized log intensity ratios across different slides. 
#   For the matrix in each of the components, rows correspond to genes
#   and columns correspond to different hybridizations, that is different 
#   slides.  This methods scale the matrix such that each column has the 
#   same median absolute deviation.
# }
# 
# \references{Y. H. Yang, S. Dudoit, P. Luu and T. P. Speed. 
#  Normalization for cDNA Microarray Data. (Statistics, UC Berkeley, 
#  Tech Report \# 589).  } 
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
# }
# 
# \seealso{\code{\link{ma.func}}, \code{\link{norm.l.func}},
#   \code{\link{norm.pin.func}}, \code{\link{norm.scale.func}}, \code{\link{plot.mva}}, \code{\link{lowess}}.}
# 
# \examples{
# data(MouseArray)
# ## mouse.setup <- init.grid() 
# ## mouse.data <- init.data() ## see \emph{init.data} 
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# mouse.norm.lratio <- stat.norm.exp(mouse.lratio$M)
# }
# 
# \keyword{microarray, log ratio, normalization.}
#*/#########################################################################

stat.norm.exp <- function(X)
  {
    n <- ncol(X)
    
    xmat.mad <- apply(X, 2, mad, na.rm=TRUE)
    
    denom <- (prod.na(xmat.mad))^(1/n)
    si <- xmat.mad / denom

    t(t(X) / si)
  }


##########################################################################
# Internal functions call by stat.ma
##########################################################################

##########################################################################
# ma.func
#
# March 20,  Remove pch="." from the code and set as an argument.
#            
##########################################################################                 
ma.func <- function (R, G, Rb, Gb, layout, norm = "p", pout, f = 0.3, extra.type="tci", crit1=0.025,crit2=crit1, nclass=10, labs=NULL, plot.type="n", col.ex=NULL, pch=pch, ...){
###
# extra.type ="t" for txt, ="p" for points, ="tci" for text ci, ="pci" for points ci ="lci" for lines confidence bands   
# crit is the size of pointwise confidence bands
# nclass (0 < nclass <1) is the proportion of points in each band i.e. smoothness of confbands
# plot.type ="n" plot normalised, ="r" raw data, ="b" both.
###                 
  if(is.null(Gb))
    cy3 <- G
  else
    cy3 <- G - Gb
  if(is.null(Rb))
    cy5 <- R
  else
    cy5 <- R - Rb

  A <- oA <- (log.na(cy3,2) +   log.na(cy5, 2))/2   #<- log.na(cy3 * cy5, 2)/2
  oM <- (log.na(cy5,2) - log.na(cy3,2))    #log.na(cy5/cy3, 2)
  if (norm == "n")
    M <- oM
  if (norm == "m")
    M <- oM - median(oM, na.rm = TRUE)
  if (norm == "l")
    M <- norm.l.func(oA, oM, f = f)
  if (norm == "p")
    M <- norm.pin.func(oA, oM, layout, f = f)
  if (norm =="s"){
    temp <- norm.pin.func(oA, oM, layout, f = f)
    M <- norm.scale.func(temp, layout)
  }
  if (pout){
    if(is.null(labs)) labs <- as.character(1:length(M))
    if(plot.type=="b") par(ask=TRUE) else par(ask=FALSE)
    if( ((plot.type == "b") | (plot.type == "r")) ){
##      par(mfrow = c(2, 1))
      plot(oA, oM, xlab = "A", ylab = "M", pch=pch, ...)
      if(extra.type== "t") text(oA,oM,labs,col=col.ex,...)
      if(extra.type== "p") points(oA,oM,col=col.ex,...)
      if(extra.type== "tci") plot.confband.text(oA,oM,crit1,crit2,nclass, labs,col=col.ex,...)
      if(extra.type== "pci") plot.confband.points(oA,oM,crit1,crit2,nclass,col=col.ex,...)
      if(extra.type== "lci") plot.confband.lines(oA,oM,crit1,crit2,nclass,col=col.ex,...)
      plot.smooth.line(oA, oM, f = 0.3)
    }
    if( ((plot.type == "b") | (plot.type =="n")) ){
      plot(A, M, xlab = "A", ylab = "Normalized M",pch=pch,...)
      if(extra.type== "t") txt(A,M,labs,col=col.ex,...)
      if(extra.type== "p") points(A,M,col=col.ex,...)
      if(extra.type== "tci") plot.confband.text(A,M,crit1,crit2,nclass,labs, col=col.ex,...)
      if(extra.type== "pci") plot.confband.points(A,M,crit1,crit2,nclass,col=col.ex,...)
      if(extra.type== "lci") plot.confband.lines(A,M,crit1,crit2,nclass,col=col.ex,...)
      plot.smooth.line(A, M, f = 0.3)
    }
    par(ask=FALSE)
##    par(mfrow = c(1, 1))
  }
  else{ list(M = M, A = A) }
}

norm.scale.func <- function(x, layout, x.names=NULL)
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

    if(is.matrix(x)){
      xmat <- x
      vect <- FALSE
                    }
    
    xmat.mad <- apply(xmat, 2, mad, na.rm=TRUE)
    sigma2 <- (1/nperpin) * exp((1/npin)*sum.na(log(xmat.mad)))
    si <- xmat.mad / (sigma2 * nperpin)
    xmat.s <- t(t(xmat) / si)

  if(vect)
      res <- as.vector(xmat.s)
    else
      res <- xmat.s
    res
  }

norm.l.func <- function(A, M, ...)
  {
    ind <- is.na(A) | is.na(M) | is.infinite(A) | is.infinite(M)
    smoothnum <- lowess(A[!ind], M[!ind], ...)
    lowesslratio <- M
    lowesslratio[!ind] <- M[!ind] - approx(smoothnum, xout = A[!ind])$y
    lowesslratio
  }

norm.pin.func <- function(A, M, layout, ...)
  {
    npin <- layout$ngrid.r * layout$ngrid.c
    pin <- c(0, rep(layout$nspot.r * layout$nspot.c, npin) * (1:npin))
    ind <- 1:length(M)
    lowessratio <- M
    for(j in 1:npin) {
      index <- ((pin[j] + 1) <= ind) & (ind <= pin[j + 1])
      tM <- M[index]
      tA <- A[index]
      ind2 <- is.na(tM) | is.na(tA) | is.infinite(tM) | is.infinite(tA)
      smoothnum <- lowess(tA[!ind2], tM[!ind2], ...)
      lowessratio[index][!ind2] <- tM[!ind2] - approx(smoothnum, xout = tA[!ind2])$y
    }
    lowessratio
  }


##########################################################################
#                                End of file
##########################################################################
