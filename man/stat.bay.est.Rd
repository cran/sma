\name{stat.bay.est}
\alias{stat.bay.est}

\title{Calculates an Odds Ratio for Each Gene in a Multi-slide Microarray Experiment. }
\description{This function takes independent, sufficient estimates of
  the effect and its variance for each gene in a multi-slide microarray experiment
and returns an odds ratio for each gene: log( Pr(the gene is
differentially expressed) / Pr(the gene is not differentially expressed)
). The parameter estimates of the Bayesian model used, as well as some
data structures which are useful when presenting the lodscore
graphically are also in the output.}

\usage{
stat.bay.est(M=NULL, Xprep=NULL, para=list(p = 0.01, v = NULL, a = NULL, c = NULL))
}

\arguments{
  \item{M}{Matrix of (normalized) log expression ratios
    \eqn{M = log_2 (R/G)} (E.g. output from stat.ma())}
  \item{Xprep}{A list containing the effect estimates and variances for
    the genes, as well as some constants needed for the odds ratio} 
  \item{para}{Estimates of the parameters used in the Bayesian
    calculations. (These are calculated only if not
    already supplied as input. See details!)}
}

%#para	are the hyperparameters needed for the calculation of lods
%#	can be produced by stat.bay.est
%#	are estimated from Xprep

%#Mbar=effect estimates for all genes (#genes x 1)
%#Vest=estimates of sigma^2 (effect variances) for all genes (#genes x 1)
%#k=constant s th Mbar~N(.,sigma^2/k) for all genes (1 x 1)
%#f=degrees of freedom for Vest (1 x 1)
\details{
  Xprep and para are optional input, but they are always in the
  output. If Xprep is supplied as input, M is unnecessary
  input. If Xprep is not supplied, stat.bay.est assumes the experiment
  consists of ncol(M) microarray slides all measuring the same effect
  (which will be stimated by Mbar). A subset of the parameters in para can be specified in the
  input, allowing the function to estimate only the others.

  Xprep is a list containing
  
  \item{Mbar}{effect estimates for all genes (\#genes x 1)}
  \item{Vest}{estimates of sigma\^2 (effect variances) for all genes (\#genes x 1)}
  \item{k}{constant so that Mbar~N(.,sigma\^2/k) for all genes (1 x 1)}
  \item{f}{degrees of freedom for Vest (1 x 1)}

  para is a list of parameters common to all genes containing

  \item{p}{Probability that a random gene is differentially
    expressed. Default is 0.01.}
  \item{v,a}{Parameters in the prior for the variance such that a*k/(2*sigma\^2) ~Gamma(v,1)}
  \item{c}{Parameter in the prior for the mean expression ratio.}
}

\value{
  A list of
 
  \item{Xprep}{Some data structures useful in graphical
    presentation. See details!} 
  \item{para}{Estimates of the parameters used in the Bayesian
    calculations. See details!}
  \item{lods}{The log odds ratio for each gene.}
}

\references{I. Lönnstedt and T. P. Speed.  Replicated Microarray Data.
  Statistical Sinica, Accepted, see http://www.stat.berkeley.edu/users/terry/zarray/Html/papersindex.html}

\author{ Ingrid Lönnstedt \email{ingrid@math.uu.se} \cr}

\seealso{\code{\link{stat.bayesian}},\code{\link{plot.bayesian}}}

\examples{
data(MouseArray)
## mouse.setup <- init.grid() 
## mouse.data <- init.data() ## see \emph{init.data} 
## mouse.lratio <- stat.ma(mouse.data, mouse.setup)

mouse.bayesian<-stat.bay.est(M=mouse.lratio$M)
plot(mouse.bayesian$Xprep$Mbar, mouse.bayesian$lods)

#alternatively

mouse.est<-apply(mouse.lratio$M,1,mean.na)
mouse.Vest<-apply(mouse.lratio$M,1,var.na)
n<-ncol(mouse.lratio$M)
k<-n
f<-n-1
mouse.Xprep<-list(Mbar=mouse.est,Vest=mouse.Vest,k=k,f=f)

mouse.bayest<-stat.bay.est(Xprep=mouse.Xprep)
plot(mouse.bayest$Xprep$Mbar, mouse.bayest$lods)

}


\keyword{models}
%\keyword{microarray, bayesian, lodsratio}








