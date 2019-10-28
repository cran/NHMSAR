\name{Mstep.nh.MSAR}
\alias{Mstep.nh.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M step of the EM algorithm.}
\description{
M step of the EM algorithm for fitting Markov switching auto-regressive models with non homogeneous transitions.}
\usage{
Mstep.nh.MSAR(data,theta,FB,covar=NULL,method=method,
ARfix=FALSE,reduct=FALSE,penalty=FALSE,sigma.diag=FALSE,sigma.equal=FALSE,
lambda1=lambda1,lambda2=lambda2,par = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
array of univariate or multivariate series with dimension T*N.samples*d.
T: number of time steps of each sample, N.samples: number of realisations of the same stationary process, d: dimension.
}
  \item{theta}{
model's parameter; object of class MSAR. See also init.theta.MSAR.
}
  \item{FB}{
Forward-Backward results, obtained by calling Estep.MSAR function}

  \item{covar}{
transitions covariates}
  \item{method}{
permits to choice the optimization algorithm. default is "ucminf", other possible choices are "BFGS" or "L-BFGS-B"}
\item{sigma.diag}{if TRUE the innovation covariance matrices are diagonal.}
  \item{sigma.equal}{If sigma.equal==TRUE the estimated covariance of the innovation will be the same in all regimes -  available only for models with homogeneous emission probabilities
  (default is FALSE)}
\item{reduct}{if TRUE, autoregressive matrices and innovation covariance matrices are constrained to have the same pattern (zero and non zero coefficients) as the one of initial matrices.   }
 \item{ARfix}{
if TRUE the AR parameters are not estimated, they stay fixed at their initial value. }
  \item{lambda1}{
penalization constant for the precision matrices. It may be a scalar or a vector of length M (with M the number of regimes). If it is equal to0 no penalization is introduced for the precision matrices. }
  \item{lambda2}{
penalization constant for the autoregressive matrices. It may be a scalar or a vector of length M (with M the number of regimes).  If it is equal to0 no penalization is introduced for the atoregression matrices.}
\item{penalty}{choice of the penalty for the autoregressive matrices. Possible values are ridge, lasso or SCAD (default). }
\item{par}{allows to give an initial value to the precision matrices. }
}

\value{
List containing
	\item{..$A0}{intercepts}
	\item{..$A}{AR coefficients}
	\item{..$sigma}{variance of innovation}
	\item{..$prior}{prior probabilities}
	\item{..$transmat}{transition matrix}
	\item{..$par.trans}{transitions parameters}

}

\references{
Ailliot P., Monbet V., (2012), Markov switching autoregressive models for wind time series. Environmental Modelling & Software, 30, pp 92-101.}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
fit.MSAR, init.theta.MSAR, Mstep.hh.MSAR}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{EM algorithm M step }

