\name{plot.bayesian}
\alias{plot.bayesian}

\title{Plots lodscore vs effect estimate (the output of stat.bayesian() or stat.bay.est()).}
\description{
  This function plots the results of stat.bayesian() or stat.bay.est()
  highlighting genes meeting userspecified criteria.
}
\usage{
plot.bayesian(x=NULL, Mbar=x$Xprep$Mbar, lods=x$lods, 
type="t", spec=50, ch=NULL, col="black",...)
}

\arguments{
  \item{x}{output from stat.bayesian() or stat.bay.est()}
  \item{Mbar}{Effect estimate for each gene.}
  \item{lods}{Posterior lodsratio for each gene.}
  \item{type}{Type of criteria for highlighting genes (see details).}
  \item{spec}{Specifies cutoff etc. depending on \code{type}.}
  \item{ch}{Character with which to selected highlight genes. (Default
    is the gene index.)}
  \item{col}{Colour of highlighted genes. (Default is black.)}
  \item{\dots}{Additional graphical parameters}
}

\details{
  If type='t', plot.bayesian() highlights the top spec genes.
  If type='c', plot.bayesian() highlights the genes with lods>=spec.
  If type='i', plot.bayesian() highlights the genes with indecies
  specified in spec.
}

\author{ Ingrid Lönnstedt \email{ingrid@math.uu.se} \cr}

\seealso{ \code{\link{stat.bayesian}},\code{\link{stat.bay.est}}}

\examples{
data(MouseArray)
## mouse.setup <- init.grid() 
## mouse.data <- init.data() ## see \emph{init.data} 
## mouse.lratio <- stat.ma(mouse.data, mouse.setup)
mouse.bayesian<-stat.bayesian(M=mouse.lratio$M)

plot.bayesian(mouse.bayesian)
plot.bayesian(mouse.bayesian,type='c',spec=0,col='blue')
plot.bayesian(mouse.bayesian,type='i',spec=c(1000,2000),col='green')
plot.bayesian(mouse.bayesian,type='i',spec=1:length(mouse.bayesian$lods))

}

\keyword{hplot}
%\keyword{ microarray, bayesian, lodsratio}





