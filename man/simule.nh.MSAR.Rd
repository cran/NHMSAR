\name{simule.nh.MSAR}
\alias{simule.nh.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of (non) homogeneous Markov Stiwtching autoregressive models
}
\description{
simule.nh.MSAR simulates realisations of (non) homogeneous Markov Switching autoregressive models with Gaussian innovations}
\usage{
simule.nh.MSAR(theta, Y0, T, N.samples = 1, covar.emis = NULL, covar.trans = NULL,
link.ct = NULL,nc = 1,S0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
list of class MSAR including model parameters and a description of the model. See init.theta.MSAR for more details. 
}
  \item{Y0}{
  	Initial value. Array of dimension order*N.samples*d with order the AR order, N.samples the number of samples to be simulated and d the dimension of the considered data.
%%     ~~Describe \code{Y0} here~~
}
  \item{T}{Length of each realisation to be simulated
%%     ~~Describe \code{T} here~~
}
  \item{N.samples}{ number of samples to be simulated
%%     ~~Describe \code{N.samples} here~~
}
  \item{covar.emis}{ emission covariate or lag for non homogeneous models. Lag is used if the covariate is the lagged time series.
%%     ~~Describe \code{covar.emis} here~~
}
  \item{covar.trans}{transition covariate or lag for non homogeneous models. Lag is used if the covariate is the lagged time series.}
  \item{link.ct}{allows to specify a link function for non homogeneous transitions.  
  	}
  	  \item{nc}{allows to specify the components of the  vector to be considered as covariates  in the non homogeneous transitions   (default is the first component).
  	}
  	  \item{S0}{initial state of the Markov chain if not null
  	}
}

\value{List including
	\item{..$Y}{simulated observation time series}
	\item{..$S}{simulated Markov chain}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~}
\author{
Val\'erie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
fit.MSAR, init.theta.MSAR,valid_all}
\examples{
data(meteo.data)
data = array(meteo.data$temperature,c(31,41,1)) 
k = 40
plot(data[,k,1],typ="l",xlab=("time (days)"),ylab=("temperature (Celsius degrees)"))
T = dim(data)[1]
N.samples = dim(data)[2]
d = dim(data)[3]
# Fit Homogeneous MS-AR models
M = 2
order = 2
theta.init = init.theta.MSAR(data,M=M,order=order,label="HH") 
mod.hh = fit.MSAR(data,theta.init,verbose=TRUE,MaxIter=20)
# Simulation
yT = 31
Bsim = 1
Ksim = Bsim*N.samples
Y0 = array(data[1:2,sample(1:dim(data)[2],Ksim,replace=T),],c(2,Ksim,1))
Y.sim = simule.nh.MSAR(mod.hh$theta,Y0 = Y0,T,N.samples = Ksim) 
# Validation
# valid_all(data,Y.sim$Y,id=1,alpha=.05)

## Not run
#data(lynx)
#lyt <- log10(lynx)
#T = length(lynx)
#Y = array(lyt,c(T,1,1))
#theta = init.theta.MSAR(Y,M=2,order=2,label='NH',nh.transitions="logistic",ncov.trans=1)
#Z = array(lyt[1:(T-2)],c(T-2,1,1))
#res=fit.MSAR(lyt[3:T],theta,covar.trans=Z,verbose=TRUE)
#Y0 = lyt[1:2]
#Bsim = 20
#Y0 = array(data[1:2,sample(1:dim(data)[2],Bsim,replace=TRUE),],c(2,Bsim,1))
#Y.sim = simule.nh.MSAR(res$theta,Y0 = Y0,T,N.samples = Bsim,covar.trans=2)}

\keyword{ Simulation }
\keyword{ MSAR }% __ONLY ONE__ keyword per line
