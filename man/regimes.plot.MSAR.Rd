\name{regimes.plot.MSAR}
\alias{regimes.plot.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot MSAR time series with regimes 
}
\description{
 Plot MSAR time series with regimes materialized by gray boxes. }
\usage{
regimes.plot.MSAR(res, data, ex = 1, col.l = "red", nc = 1, 
ylim = NULL, xlab = "time", ylab = "series", d = NULL, dt = 1, lwd = 1, cex = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
list obtained from fit.MSAR fonction as result of MSAR fitting}
  \item{data}{
data to plot}
  \item{ex}{
number of sample}
  \item{nc}{component number (useful for multivariate time series)
}
  \item{col.l}{color of time series (default is red)
}

  \item{ylim}{
range for the plotted 'y' values, defaulting to the range of the finite values of 'y'}
  \item{xlab}{
a title for the x axis}
  \item{ylab}{
a title for the y axis}
  \item{d}{dimension to be plot (for multivariate cases). Default is 1. 
%%     ~~Describe \code{d} here~~
}
  \item{dt}{
time step (default=1)}
  \item{lwd}{
width of the line}
  \item{cex}{
symbols/text size}
}


\value{Returns a plot and the regimes time series. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(lynx)
T = length(lynx)
data = array(log(lynx),c(T,1,1))
theta.init = init.theta.MSAR(data,M=2,order=2,label="HH")
mod.lynx = fit.MSAR(data,theta.init)
regimes.plot.MSAR(mod.lynx,data,ylab="Captures number")

theta.init = init.theta.MSAR(data,M=2,order=2,label="NH",nh.transitions="logistic")
attributes(theta.init)
theta.init$A0 = mod.lynx$theta$A0
theta.init$A = mod.lynx$theta$A
theta.init$sigma = mod.lynx$theta$sigma
theta.init$prior = mod.lynx$theta$prior
theta.init$transmat = mod.lynx$theta$transmat
theta.init$par.trans = matrix(c(1,-1,-.2,.2),2,2)
Y = array(data[2:T,,],c(T-1,1,1))
Z = array(data[2:T,,],c(T-1,1,1))
mod.lynx = fit.MSAR(Y,theta.init,verbose=TRUE,MaxIter=20,covar.trans=Z)
regimes.plot.MSAR(mod.lynx,data,ylab="Captures number")

## Not run
# Fit Homogeneous MS-AR models - multivariate time series
#data(PibDetteDemoc)
#T = length(unique(PibDetteDemoc$year))-1
#N.samples = length(unique(PibDetteDemoc$country))
#PIB = matrix(PibDetteDemoc$PIB,N.samples,T+1)
#Dette = matrix(PibDetteDemoc$Dette,N.samples,T+1)
#Democratie = matrix(PibDetteDemoc$Democratie,N.samples,T+1)

#d = 2
#Y = array(0,c(T,N.samples,2))
#for (k in 1:N.samples) {
#   Y[,k,1] = diff(log(PIB[k,]))
#   Y[,k,2] = diff(log(Dette[k,]))
#}
#Democ = Democratie[,2:(T+1)] 
#theta.hh.1 = init.theta.MSAR(Y,M=4,order=1,label="HH")
#res.hh = fit.MSAR(Y,theta.hh.1,verbose=TRUE,MaxIter=200)
#par(mfrow=c(2,1))
#regimes.plot.MSAR(res.hh,Y,ex=30,ylab="GDP")
#regimes.plot.MSAR(res.hh,Y,ex=30,nc=2,ylab="Debt")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
\keyword{ MSAR }
\keyword{ latent regimes}% __ONLY ONE__ keyword per line
