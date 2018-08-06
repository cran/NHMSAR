\name{prediction.MSAR}
\alias{prediction.MSAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One step ahead predict for (non) homogeneous MSAR models}
\description{
computes one step ahead predict for (non) homogeneous MSAR models. A time series is given as input and a prediction is return for each time. These function is mainly usefull for cross-validation.  }
\usage{
prediction.MSAR(data, theta, covar.emis = NULL, covar.trans = NULL, ex = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
observed time series, array of dimension T*N.samples*d}
  \item{theta}{
object of class MSAR including the model's parameter}
  \item{covar.emis}{
covariate for emissions (if needed)}
  \item{covar.trans}{
covariate for transitions (if needed)}
  \item{ex}{
numbers of samples for which prediction has to be computed}
}

\value{
Returns a list with the following elements:
  \item{y.p}{the one step ahead prediction for each time of data time series}
  \item{var.p}{the associated variance}
  \item{pr}{the prediction probabilities for each regime }
}
%\references{}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Cond.prob.MSAR
}
\examples{
## Not run
#data(meteo.data)
#data = array(meteo.data$temperature,c(31,41,1)) 
#T = dim(data)[1]
#N.samples = dim(data)[2]
#d = dim(data)[3]
#M = 2
#theta.init = init.theta.MSAR(data,M=M,order=2,label="HH")
#res.hh.2 = fit.MSAR(data,theta.init,verbose=TRUE,MaxIter=200)
#y.p.2 = prediction.MSAR(data,res.hh.2$theta,ex=1:N.samples)
#RMSE.2 = mean((data-y.p.2$y.p)^2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MSAR }
\keyword{ prediction }
\keyword{ forecast }
\keyword{ cross-validation }
% __ONLY ONE__ keyword per line
