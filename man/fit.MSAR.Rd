\name{fit.MSAR (NH-MSAR)}
\alias{fit.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit (non) homogeneous Markov switching autoregressive models}
\description{
Fit (non) homogeneous Markov switching autoregressive models by EM algorithm. Non homogeneity may be introduced at the intercept level or in the probability transitions. The link functions are defined in the initialisation step (running init.theta.MSAR.R). 
}
\usage{
fit.MSAR(data, theta, MaxIter = 100, eps = 1e-05, verbose = FALSE, 
   covar.emis = NULL, covar.trans = NULL, method = NULL, 
   constraints = FALSE, reduct=FALSE, K = NULL, d.y = NULL, 
   ARfix = FALSE,penalty=FALSE,sigma.diag=FALSE, sigma.equal=FALSE,
   lambda1=.1,lambda2=.1,a=3.7,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
array of univariate or multivariate series with dimension T*N.samples*d.
T: number of time steps of each sample, N.samples: number of realisations of the same stationary process, d: dimension.
}
  \item{theta}{
initial parameter obtained running function init.theta.MSAR.R; object of class MSAR.
}
  \item{MaxIter}{
maximum number of iteration for EM algorithm (default : 100)}
  \item{eps}{ Tolerance for likelihood.
%%     ~~Describe \code{eps} here~~
}
  \item{verbose}{
if verbose=TRUE, the value of log-likelihood is printed at each EM-algorithm's iteration}
  \item{covar.emis}{
array of univariate or multivariate series of covariate to take into account in the intercept of the autoregressive models. 
The link function is defined in the initialisation step (running init.theta.MSAR.R). }
  \item{covar.trans}{
array of univariate or multivariate series of covariate to take into account in the transition probabilities.
The link function is defined in the initialisation step (running init.theta.MSAR.R).}
  \item{method}{
permits to choice the optimization algorithm if numerical optimisation is required in M step. Default : "ucminf". Other choices : "L-BFGS-B", "BFGS"}
  \item{constraints}{
if constraints = TRUE constraints are added to theta in order that matrices A and sigma are diagonal by blocks. }
  \item{K}{
number of sites. For instance, if one considers wind at k locations, K=k. Or more generally number of independent groups of components.  }
  \item{d.y}{
dimension in each sites. For instance, if one considers only wind intensity than d.y = 1; but, if one considers cartesian components of wind, then d.y =2.}
\item{ARfix}{if TRUE the AR parameters are not estimated, they stay fixed at their initial value.}
\item{reduct}{if TRUE, autoregressive matrices and innovation covariance matrices are constrained to have the same pattern (zero and non zero coefficients) as the one of initial matrices.   }
  \item{sigma.diag}{If sigma.diag==TRUE the estimated covariance of the innovation will be diagonal 
  (default is FALSE) -  available only for HH models 
  }
  \item{sigma.equal}{If sigma.equal==TRUE the estimated covariance of the innovation will be the same in all regimes -  available only for models with homogeneous emission probabilities
  (default is FALSE)}
  \item{penalty}{choice of the penalty for the autoregressive matrices. Possible values are ridge (available for regression matrices only), lasso or SCAD (default). }
 \item{lambda1}{
penalization constant for the precision matrices. It may be a scalar or a vector of length M (with M the number of regimes). If it is equal to 0 no penalization is introduced for the precision matrices. }
  \item{lambda2}{
penalization constant for the autoregressive matrices (available only for MSAR of order 1). It may be a scalar or a vector of length M (with M the number of regimes). If lambda2 is a scalar, it is weighted by the prior in each regime.  }
\item{a}{fixed penalisation constant for SCAD penalty}
 \item{...}{other arguments}

}
\details{
The homogeneous MSAR model is labeled "HH" and it is written
 \deqn{
P(X_t|X_{t-1}=x_{t-1}) = Q_{x_{t-1},x_t}} with \eqn{X_t} the hidden univariate process defined on \eqn{\{1,\cdots,M \}}
 \deqn{
Y_t|X_t=x_t,y_{t-1},...,y_{t-p} = \alpha_0^{x_t}+\alpha_1^{x_t}y_{t-1}+...+\alpha_p^{x_t}y_{t-p}+\sigma \epsilon_t} with \eqn{Y_t} the observed process and \eqn{\epsilon} a Gaussian white noise. \eqn{Y_t} may be mutivariate. 

The model with non homogeneous emissions is labeled "HN" and it is written
 \deqn{P(X_t|X_{t-1}=x_{t-1}) = Q_{x_{t-1},x_t}} with \eqn{X_t} the hidden process
 \deqn{Y_t|X_t=x_t,y_{t-1},...,y_{t-p} = f(z_t,\theta_z^{x_t})+\alpha_1^{x_t}y_{t-1}+...+\alpha_p^{x_t}y_{t-p}+\sigma \epsilon_t} with \eqn{Y_t} the observed process, \eqn{\epsilon} a Gaussian white noise and \eqn{Z_t} a covariate. 

The model with non homogeneous transitions is labeled "NH" and it is written
\deqn{P(X_t|X_{t-1}=x_{t-1}) = q(z_t,\theta_{z_t})} with \eqn{X_t} the hidden process and \eqn{q} a link function which has a Gaussian shape by default. 
\deqn{Y_t|X_t=x_t,y_{t-1},...,y_{t-p} = \alpha_0^{x_t}+\alpha_1^{x_t}y_{t-1}+...+\alpha_p^{x_t}y_{t-p}+\sigma \epsilon_t} with \eqn{Y_t} the observed process, \eqn{\epsilon} a Gaussian white noise and \eqn{Z_t} a covariate.}
\value{
	For fit.MSAR and its methods a list of class "MSAR" with the following elements:
	
	
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list including:
\item{..$theta}{object of class MSAR containing the estimated values of the parameter and some descriptors of the fitted model. See init.theta.MSAR for a detailled description.   }
\item{..$ll_history}{log-likelihood for each iterations of the EM algorithm. }
\item{..$Iter}{number of iterations run before EM converged}
\item{..$Npar}{number of parameters in the model}
\item{..$BIC}{Bayes Information Criterion}
\item{..$smoothedprob}{smoothing probabilities \eqn{P(X_t|y_0,\cdots,y_T)}}


Penalized likelihood is considered if at least one of the lambdas parameters are non zero. When LASSO penalty is chosen, the LARS algorithm is used. When SCAD is chosen, a Newton-Raphson algorithm is run with a quadratic approximation of the penalized likelihood. 
	For the precision matrices penalization, the package glasso is used. 
Limit of this function: likelihood penalization only works for VAR(1) models
}
\references{
	Ailliot P., Monbet V., (2012), Markov-switching autoregressive models for wind time series. Environmental Modelling & Software, 30, pp 92-101. 
	Efron, B., Hastie, T., Johnstone, I., Tibshirani, R., et al. (2004). Least angle regression. The Annals of statistics, 32(2):407-499.

Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties. Journal of the American statistical Association, 96(456):1348-1360. 
	Hamilton J.D. (1989). A New Approach to the Economic Analysis of Nonstionary Time Series and the Business Cycle. Econometrica 57: 357-384.
}
\author{
Valerie Monbet, valerie.monbet at univ-rennes1.fr}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{init.theta.MSAR, regimes.plot.MSAR, simule.nh.ex.MSAR, depmixS4, MSBVAR
}
\examples{
# Fit Homogeneous MS-AR models - univariate time series
data(meteo.data)
data = array(meteo.data$temperature,c(31,41,1)) 
k = 40
T = dim(data)[1]
N.samples = dim(data)[2]
d = dim(data)[3]
M = 2
order = 2 # MSAR of order 2
theta.init = init.theta.MSAR(data,M=M,order=order,label="HH") 
mod.hh = fit.MSAR(data,theta.init,verbose=TRUE,MaxIter=20)
#regimes.plot.MSAR(mod.hh,data,ylab="temperatures") 
#Y0 = array(data[1:2,sample(1:dim(data)[2],1),],c(2,1,1))
#Y.sim = simule.nh.MSAR(mod.hh$theta,Y0 = Y0,T,N.samples = 1)

## Not run
# Fit Non Homogeneous MS-AR models - univariate time series
#data(lynx)
#T = length(lynx)
#data = array(log10(lynx),c(T,1,1))
#theta.init = init.theta.MSAR(data,M=2,order=2,label="HH")
#mod.lynx.hh = fit.MSAR(data,theta.init,verbose=TRUE,MaxIter=200)
#regimes.plot.MSAR(mod.lynx.hh,data,ylab="Captures number")

#theta.init = init.theta.MSAR(data,M=2,order=2,label="NH",nh.transitions="logistic")
attributes(theta.init)
#theta.init$A0 = mod.lynx.hh$theta$A0
#theta.init$A = mod.lynx.hh$theta$A
#theta.init$sigma = mod.lynx.hh$theta$sigma
#theta.init$transmat = mod.lynx.hh$theta$transmat
#theta.init$prior = mod.lynx.hh$theta$prior
#Y = array(data[2:T,,],c(T-1,1,1))
#Z = array(data[1:(T-1),,],c(T-1,1,1))
#mod.lynx = fit.MSAR(Y,theta.init,verbose=TRUE,MaxIter=200,covar.trans=Z)
#regimes.plot.MSAR(mod.lynx,Y),ylab="Captures number")

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
#theta.hh = init.theta.MSAR(Y,M=M,order=1,label="HH")
#res.hh = fit.MSAR(Y,theta.hh,verbose=TRUE,MaxIter=200)
#regime.hh = apply(res.hh$smoothedprob,c(1,2),which.max)

## Not run
# Fit Non Homogeneous (emission) MS-AR models - multivariate time series
#theta.hn = init.theta.MSAR(Y,M=M,order=1,label="HN",ncov.emis=1)
#theta.hn$A0 = res.hh$theta$A0
#theta.hn$A = res.hh$theta$A
#theta.hn$sigma = res.hh$theta$sigma
#theta.hn$transmat = res.hh$theta$transmat
#theta.hn$prior = res.hh$theta$prior
#Z = array(t(Democ[,2:T]),c(T,N.samples,1))
#res.hn = fit.MSAR(Y,theta.hn,verbose=TRUE,MaxIter=200,covar.emis=Z)

# Fit Non Homogeneous (transitions) MS-AR models - multivariate time series
#theta.nh = init.theta.MSAR(Y,M=M,order=1,label="NH",nh.transitions="gauss",ncov.trans=1)
#theta.nh$A0 = res.hh$theta$A0
#theta.nh$A = res.hh$theta$A
#theta.nh$sigma = res.hh$theta$sigma
#theta.nh$transmat = res.hh$theta$transmat
#theta.nh$prior = res.hh$theta$prior
#theta.nh$par.trans[1:2,1] = 10
#theta.nh$par.trans[3:4,1] = 0
#theta.nh$par.trans[,2] = 2
#Z = array(t(Democ[,2:T]),c(T,N.samples,1))
#res.nh = fit.MSAR(Y,theta.nh,verbose=TRUE,MaxIter=200,covar.trans=Z)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MSAR }
\keyword{ Model fitting }
\keyword{ Maximum likelihood }% __ONLY ONE__ keyword per line
