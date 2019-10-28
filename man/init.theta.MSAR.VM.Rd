\name{init.theta.MSAR.VM}
\alias{init.theta.MSAR.VM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initialisation function for von Mises MSAR model fitting}
\description{
Initialization before fitting von Mises (non) homogeneous Markov switching autoregressive models by EM algorithm. Non homogeneity may be introduce in the probability transitions. The link function is defined here.}
\usage{
init.theta.MSAR.VM(data, ..., M, order, 
                  regime_names = NULL, 
                  nh.emissions = NULL, nh.transitions = NULL, 
                  label = NULL, ncov.emis = 0, ncov.trans = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{
array of univariate or multivariate series with dimension T*N.samples*d with 
T: number of time steps of each sample, N.samples: number of realisations of the same stationary process, d: dimension}
 
  \item{M}{
number of regimes}
  \item{order}{
order of AR processes}
  \item{label}{
"HH" (default) for homogeneous MS AR model 
"NH" for non homogeneous transitions  }
  \item{regime_names}{
(optional) regime's names may be chosen}
  \item{nh.emissions}{
not available - under development.   }
  \item{nh.transitions}{
link function for non homogeneous transitions. Default: von Mises (see details).   }

  \item{ncov.emis}{
not available - under development. }
  \item{ncov.trans}{
number of covariates in NH model
}
 \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
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
return a list of class MSAR including 
	
\item{theta}{parameter}
\item{..$transmat}{transition matrix}
\item{..$prior}{prior probabilities}
\item{..$mu}{vector of intercepts}
\item{..$kappa}{matrix of 'AR' coefficients (not complex by default)}
\item{..$par.emis}{parameters of non homogeneous emissions (not used)} 
\item{..$par.trans}{parameters of non homogeneous transitions} 
\item{label}{model's label}
}

\references{
Ailliot P., Bessac J., Monbet V., Pene F., (2014) Non-homogeneous hidden Markov-switching models for wind time series. JSPI.}
\author{
Val\'erie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
fit.MSAR.VM}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ parameter initialization}

