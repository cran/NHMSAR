\name{fit.MSAR.VM}
\alias{fit.MSAR.VM}
\title{
Fit von Mises (non) homogeneous Markov switching autoregressive models}

\description{
Fit von Mises (non) homogeneous Markov switching autoregressive models by EM algorithm. Non homogeneity may be introduced at the intercept level or in the probability transitions. The link functions are defined in the initialisation step (running init.theta.MSAR.VM.R). 
}
\usage{
fit.MSAR.VM(data, theta, 
           MaxIter = 100, eps = 1e-05, verbose = FALSE, 
           covar.emis = NULL, covar.trans = NULL, 
           method = NULL, constr = 0, ...)
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
  \item{constr}{
if constr = 1 constraints are added to theta  }
  \item{...}{other arguments}
}
\details{
The homogeneous MSAR model is labeled "HH" and it is written
 \deqn{
P(X_t|X_{t-1}=x_{t-1}) = Q_{x_{t-1},x_t}} with \eqn{X_t} the hidden univariate process defined on \eqn{\{1,\cdots,M \}}
 \deqn{
Y_t|X_t=x_t,y_{t-1},...,y_{t-p}} has a von Mises distribution with density \deqn{ p_2(y_t|x_t,y_{t-s}^{t-1}) = \frac {1}{b(x_t,y_{t-s}^{t-1})} \exp\left(\kappa_0^{(x_t)} \cos(y_t-\phi_0^{(x_t)})+
    \sum_{\ell=1}^s\kappa_\ell^{(x_t)} \cos(y_t-y_{t-\ell}-\phi_\ell^{(x)})\right)}
    which is equivalent to
\deqn{ p_2(y_t|x_t,y_{t-s}^{t-1}) =\frac{1}{b(x_t,y_{t-s}^{t-1})} \left|\exp\left([\gamma_0^{(x_t)} +
    \sum_{\ell=1}^s\gamma_\ell^{(x_t)} e^{iy_{t-\ell}}]e^{-iy_t}\right)\right|}
    \eqn{b(x_t,y_{t-s}^{t-1})} is a normalization constant. 
    
    Both the real and the complex formulation are implemented. In practice, the complex version is used if the initial \eqn{\kappa} is complex. 
    
    The model with non homogeneous transitions is labeled "NH" and it is written
\deqn{P(X_t|X_{t-1}=x_{t-1}) = q(z_t,\theta_{z_t})} with \eqn{X_t} the hidden process and \eqn{q} von Mises link function such that \deqn{p_1(x_t|x_{t-1},z_{t}) =\frac{  q_{x_{t-1},x_t}\left|\exp 
   \left(\tilde\lambda_{x_{t-1},x_t} e^{-iz_{t}} \right)\right|}
    {\sum_{x'=1}^M q_{x_{t-1},x'}\left|\exp 
   \left(\tilde\lambda_{x_{t-1},x'} e^{-iz_{t}} \right)\right|},
} 
with \eqn{\tilde\lambda_{x,x'}} a complex parameter (by taking \eqn{\tilde\lambda_{x,x'}=\lambda_{x,x'}
    e^{i\psi_{x,x'}}}). 
}
\value{
	For fit.MSAR and its methods a list of class "MSAR" with the following elements:
	
	
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list including:
\item{..$theta}{object of class MSAR containing the estimated values of the parameter and some descriptors of the fitted model. See init.theta.MSAR.VM for a detailled description.   }
\item{..$ll_history}{log-likelihood for each iterations of the EM algorithm. }
\item{..$Iter}{number of iterations run before EM converged}
\item{..$Npar}{number of parameters in the model}
\item{..$BIC}{Bayes Information Criterion}
\item{..$smoothedprob}{smoothing probabilities \eqn{P(X_t|y_0,\cdots,y_T)}}
}

\references{
Ailliot P., Bessac J., Monbet V., Pene F., (2014) Non-homogeneous hidden Markov-switching models for wind time series. JSPI.}
\author{
Valerie Monbet, valerie.monbet at univ-rennes1.fr}

\seealso{init.theta.MSAR.VM, regimes.plot.MSAR
}
\examples{
	## Not run
# data(WindDir)
# T = dim(WindDir)[1]
# N.samples = dim(WindDir)[2]
# Y = array(WindDir,c(T,N.samples,1))
# von Mises homogeneous MSAR
# M = 2
# order = 2
# theta.init = init.theta.MSAR.VM(Y,M=M,order=order,label="HH")
# res.hh = fit.MSAR.VM(Y,theta.init,MaxIter=3,verbose=TRUE,eps=1e-8)
## von Mises non homogeneous MSA
# theta.init = init.theta.MSAR.VM(Y,M=M,order=order,label="NH",ncov=1,nh.transitions="VM")
#theta.init$mu = res.hh$theta$mu
#theta.init$kappa = res.hh$theta$kappa
#theta.init$prior = res.hh$theta$prior
#theta.init$transmat = res.hh$theta$transmat
#theta.init$par.trans = matrix(c(res.hh[[M]][[order+1]]$theta$mu,.1*matrix(1,M,1)),2,2)
#Y.tmp = array(Y[2:T,,],c(T-1,N.samples,1))
#Z = array(Y[1:(T-1),,],c(T-1,N.samples,1))
#res.nh = fit.MSAR.VM(Y.tmp,theta.init,MaxIter=10,verbose=T,eps=1e-8,covar.trans=Z)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ von Mises MSAR }
