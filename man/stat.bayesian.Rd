\name{stat.bayesian}
\alias{stat.bayesian}

\title{Calculates an Odds Ratio of Each Gene in a Multi-slide microarray
Experiment}
\description{
This function takes the normalized expression estimates from a
multi-slide microarray experiment (M-values output
by stat.ma) and returns an odds ratio for each gene: log( Pr(the gene is
differentially expressed) / Pr(the gene is not differentially expressed)
). The parameter estimates of the Bayesian model used, as well as some
data structures which are useful when presenting the lodscore
graphically are also in the output.}

\usage{
	stat.bayesian(M=NULL, nb=NULL, nw=1, Xprep=NULL, 
	para=list(p = 0.01, v = NULL, a = NULL, c = NULL, k = NULL))
}

\arguments{
  \item{M}{Matrix of (normalized) log expression ratios
    \eqn{M = log_2 (R/G)} (E.g. output from stat.ma())}
  \item{nb}{Number of slides containing spots for a gene (common for all
  genes). Unnecessary argument if nw=1.}
  \item{nw}{Number of spots for a gene within each slide (common to all
    genes).Default is 1.}
  \item{Xprep}{Some data structures useful in graphical
    presentation. (This is calculated only if not already supplied as
    input. See details!)} 
  \item{para}{Estimates of the parameters used in the Bayesian
    calculations. (These are calculated only if not
    already supplied as input. See details!)}
}

\details{

  Xprep and para are optional input, but they are always in the
  output. If Xprep is supplied as input, M, nb and nw are unnecessary
  input. A subset of the parameters in para can be specified in the
  input, allowing the function to estimate only the others.

  Xprep is a list containing
  
  \item{nb}{Number of slides containing spots for a gene (common to all
  genes).}
  \item{nw}{Number of spots for a gene within each slide (common to all genes).}
  \item{Mbar}{Overalll means for each gene.}
  \item{SSB}{Sum of squares between slides for each gene.}
  \item{SSB}{Sum of squares within slides for each gene.}

  para is a list of parameters common to all genes containing

  \item{p}{Probability that a random gene is differentially
    expressed. Default is 0.01.}
  \item{v,a}{Parameters in the prior for the variance. If the observed
    sums of squares between slides are MSB, then v and a is such that
    (nb*a/(2MSB)~Gamma(v,1))}
  \item{c}{Parameter in the prior for the mean expression ratio.}
  \item{k}{An approximate ratio of the variances between to within slides,
    only needed when both nw>1 and nb>1}

}

\value{
  
  A list of
 
  \item{Xprep}{Some data structures useful in graphical
    presentation. See details!} 
  \item{para}{Estimates of the parameters used in the Bayesian
    calculations. See details!}
  \item{lods}{The log odds ratio for each gene.}
}


\author{
    Ingrid Lönnstedt \email{ingrid@math.uu.se} \cr
    Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
}

\seealso{\code{\link{plot.bayesian}}}

\examples{
data(MouseArray)
## mouse.setup <- init.grid() 
## mouse.data <- init.data() ## see \emph{init.data} 
mouse.lratio <- stat.ma(mouse.data, mouse.setup)

mouse.bayesian<-stat.bayesian(M=mouse.lratio$M)
plot(mouse.bayesian$Xprep$Mbar, mouse.bayesian$lods)

}

\keyword{models}

%\keyword{microarray, bayesian, lodsratio}
