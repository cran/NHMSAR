\name{cor.MSAR}
\alias{cor.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Empirical correlation functions comparison .}
\description{
Empirical correlation function of observed data and simulated data are plotted on the same figure. 
A fluctuation interval of simulations is added to help the comparison. 
}
\usage{
cor.MSAR(data,data.sim,lag=NULL,nc=1,alpha=.05,plot=FALSE,
   xlab="Time (days)",dt=1,ylab="Correlation",...)
}
\arguments{
  \item{data}{
observed (or reference) time series, array of dimension T*N.samples*d }
  \item{data.sim}{
  simulated time series, array of dimension T*N.sim*d.
N.sim have to be  K*N.samples with K large enough (for instance, K=100) }

  \item{lag}{
maximum lag at which to calculate the empirical auto-correlation function. 
Default floor(T/2) with T the length of each data sample. }
  \item{nc}{
number of component for which to calculate the empirical auto-correlation function. }
  \item{alpha}{
confidence level for computation of the fluctuation interval. Default= 0.05.}
  \item{plot}{if plot is TRUE plots are drawn (default is FALSE).}
  \item{xlab}{x axis label}
  	\item{dt}{default time step is equal to 1}
  \item{ylab}{y axis label}
    \item{...}{for optional plot arguments}

}
\details{The auto-correlation functions are computed from one or several independent realizations of the same length. 
}
\value{
	A list with the following elements:
	\item{C.data}{observed data acf}
	\item{C.sim}{simulated data acf}
	\item{CI.sim}{fluctuation interval for each lag}
	\item{lags}{abscissa for acfs}

}
\references{
Bessac, J., Ailliot, P., Monbet, V. (2015). Gaussian linear state-space model for wind fields in the North-East Atlantic. 
Environmetrics, 26(1), 29-38.}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}



\seealso{
cross.cor.MSAR, cor 
}
\examples{
## Not run
#data(Wind)
#T = dim(U)[1]
#N.samples = dim(U)[2]
#Y = array(U[,,1],c(T,N.samples,1))

#theta.init=init.theta.MSAR(Y,M=2,order=1,label="HH")
#res.hh = fit.MSAR(Y,theta.init,verbose=TRUE,MaxIter=10)
#Bsim = 2
#Ksim = Bsim*N.samples
#Y0 = array(Y[1,sample(1:dim(Y)[2],1,replace=T),],c(2,Ksim,1))
#Y.sim = simule.nh.MSAR(res.hh$theta,Y0 = Y0,T,N.samples = Ksim) 
#c = cor.MSAR(Y,Y.sim$Y)
#plot(c$lags/4,c$C.data,typ="l",xlab="Time (days)",ylab="ACF",xlim=c(0,8))
#abline(h=0,lty=3,col="gray")
#lines(c$lags/4,c$C.sim,col="red")
#lines(c$lags/4,c$CI.sim[1,],col="red",lty=2)
#lines(c$lags/4,c$CI.sim[2,],col="red",lty=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Auto-correlation function}
\keyword{ Validation }% __ONLY ONE__ keyword per line
