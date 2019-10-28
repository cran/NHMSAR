\name{Mstep.classif}
\alias{Mstep.classif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fit an AR model for each class of C}
\description{
fit an AR model for each class of C by maximum likelihood method.}
\usage{
Mstep.classif(data, C, order,sigma.diag=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
array of univariate or multivariate series with dimension T*N.samples*d.
T: number of time steps of each sample, N.samples: number of realisations of the same stationary process, d: dimension.
}
  \item{C}{
Class sequence}
  \item{order}{
order of AR models (all models will have the same order)}
  \item{sigma.diag}{
if TRUE the covariance matrices will be diagonal (default FALSE)}
}

\value{list containing
	\item{A0}{intercept}
	\item{A}{AR coefficients}
	\item{sigma}{variance of innovation}
	\item{LL}{log likelihood}
}
%\references{%% ~put references to the literature/web site here ~}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
fit.MSAR}
\examples{
	data(meteo.data)
data = array(meteo.data$temperature,c(31,41,1)) 
k = 40
T = dim(data)[1]
N.samples = dim(data)[2]
d = dim(data)[3]
order = 2
C = array(meteo.data>0,c(31,41,1))
res = Mstep.classif(data,C,order=order)
str(res)
}% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MSAR with a priori classification }
