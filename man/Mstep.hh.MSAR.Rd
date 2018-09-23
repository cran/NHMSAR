\name{Mstep.hh.MSAR}
\alias{Mstep.hh.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M step of the EM algorithm for fitting  homogeneous Markov switching auto-regressive models.}
\description{
M step of the EM algorithm for fitting   homogeneous Markov switching auto-regressive models, called in fit.MSAR.}
\usage{
Mstep.hh.MSAR(data, theta, FB,sigma.diag=FALSE,sigma.equal=FALSE)
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
Forward-Backward results, obtained by calling Estep.MSAR function
}
  \item{sigma.diag}{If sigma.diag==TRUE the estimated covariance of the innovation will be diagonal 
  (default is FALSE) -  available only for HH models 
  }
  \item{sigma.equal}{If sigma.equal==TRUE the estimated covariance of the innovation will be the same in all regimes -  available only for HH models
  (default is FALSE)
  }
}

\value{A list containing
	\item{A0}{intercepts}
	\item{A}{AR coefficients}
	\item{sigma}{variance of innovation}
	\item{prior}{prior probabilities}
	\item{transmat}{transition matrix}

}
\references{
Ailliot P., Monbet V., (2012), Markov switching autoregressive models for wind time series. Environmental Modelling & Software, 30, pp 92-101.}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}

\seealso{fit.MSAR, Estep.MSAR, Mstep.classif
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }
\keyword{ M step }