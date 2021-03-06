%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file was automatically generated by Rdoc on 2001-04-20 20:18:16.
% Do NOT modify this file, instead modify the source (sma/R/Rmva.R).
% Rdoc is copyright (C) 2001, Henrik Bengtsson, henrikb@braju.com.
% Download the Rdoc compiler at http://www.braju.com/R/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plot.confband.text}

\alias{plot.confband.text}

\title{Add Selected Text to an M vs A Plot}
  
\description{`text' draws the strings given in the vector `labs' at the
  coordinates given by `M' and `A'}

\usage{
plot.confband.text(x, M, crit1=0.025, crit2=crit1, nclass=10,
labs=NULL, output=FALSE, ...)  
}

\arguments{
 \item{x}{a vector giving the x-coordinates of the points in the scatter
          plot. In the microarray context, this could be a vector of
          average log intensities. ie A} 
 
 \item{M}{a vector giving the y-coordinates of the points in the scatter
   plot. In the microarray context, this could be a vector of log
   intensity ratios.} 
	
 \item{crit1}{The number of points to be highlighted on the M vs A plot.
   If crit1 < 1, the crit1*100\% spots with the smallest M values
   will be highlighted. If crit1 >= 1, the crit spots  with the
   smallest M values are highlighted.}  
 \item{crit2}{Similar to "crit1".   If crit2 < 1, the crit2*100\%
 spots with the largest M values will be highlighted. If crit2 >= 1,
 the crit2 spots with the largest M values are highlighted.}  
 \item{nclass}{A single number giving the approximate number of
 intensity depedent groups to consider.} 
 \item{labs}{ one or more character strings or expressions specifying the
   text to be written.  If this string is not specified, by
   default the index of the vector `M' will be used.}
 \item{output}{logical, defaulting to `FALSE'. If `TRUE' a vector
   containning the index to the vector `M' that are  
   highlighted.} 
 \item{\dots}{graphical parameters may also be supplied as arguments
 to the function (see \code{\link{par}}).}        
}

\note{An M vs A plot must be constructed \bold{prior} to the execution of this function.}

\examples{data(MouseArray)
## mouse.setup <- init.grid()
## mouse.data <- init.data()

plot.mva(mouse.data, mouse.setup) ## an M vs A plot

mouse.lratio <- stat.ma(mouse.data, mouse.setup)

plot.confband.text(mouse.lratio$A, mouse.lratio$M)
## 2.5\% of the spots with the largest and smallest M values are
## highlighted on the M vs A plot, and each spot is assigned the
## default label of its corresponding index value. 
}

\seealso{ \code{\link{plot.mva}}, \code{\link{stat.ma}},
  \code{\link{lines}}, \code{\link{matlines}},
  \code{\link{plot.confband.lines}}, \code{\link{plot.confband.points}} .
}

\keyword{aplot}
%\keyword{microarray, point-wise confidence band.}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
