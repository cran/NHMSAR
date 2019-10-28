\name{Mstep.nh.MSAR.VM}
\alias{Mstep.nh.MSAR.VM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M step of the EM algorithm for von Mises MSAR models}
\description{
M step of the EM algorithm for fitting von Mises Markov switching auto-regressive models with non homogeneous transitions.}
\usage{
Mstep.nh.MSAR.VM(data, theta, FB, covar.trans = NULL, method = method, constr = 0)
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

  \item{covar.trans}{
transitions covariates}
  \item{method}{
permits to choice the optimization algorithm. default is "ucminf", other possible choices are "BFGS" or "L-BFGS-B"}
  \item{constr}{
if constr=1 contraints are added the the \eqn{kappa} parameters}
}

\value{
List containing
	\item{mu}{intercepts}
	\item{kappa}{von Mises AR coefficients}
	\item{prior}{prior probabilities}
	\item{transmat}{transition matrix}
		\item{..$par.trans}{transitions parameters}
}
\references{
Ailliot P., Bessac J., Monbet V., Pene F., (2014) Non-homogeneous hidden Markov-switching models for wind time series. JSPI.}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}
%\note{
%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
fit.MSAR.VM, init.theta.MSAR.VM, Mstep.hh.MSAR.VM}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, theta, FB, covar = covar.trans, method = method, 
    constr = 0) 
{
    order = attributes(theta)$order
    d = dim(data)[3]
    if (is.na(d) | is.null(d)) {
        d = 1
    }
    M = attributes(theta)$NbRegimes
    if (length(covar) == 1) {
        Lag = covar
        covar = array(data[(1):(T - Lag + 1), , ], c(T - Lag + 
            1, N.samples, d))
        data = array(data[Lag:T, , ], c(T - Lag + 1, N.samples, 
            d))
    }
    N.samples = dim(covar)[2]
    ncov.trans = dim(covar)[3]
    par.hh = Mstep.hh.MSAR.VM(data, theta, FB, constr)
    theta$transmat[which(theta$transmat < 1e-15)] = 1e-15
    theta$transmat = mk_stochastic(theta$transmat)
    trans = para_trans(theta$transmat)
    par.trans = theta$par.trans
    nh_transition = attributes(theta)$nh.transitions
    par.init = plie2(trans, par.trans)
    lxi = dim(FB$probSS)[3]
    if (order > 0) {
        deb = order + 1
    }
    else {
        deb = 1
    }
    resopt = ucminf(par.init, fn = loglik_nh_inp.VM, gr = NULL, 
        covar = array(covar[deb + (1:(lxi)), , ], c(lxi, N.samples, 
            ncov.trans)), xi = FB$probSS, nh_transition = nh_transition, 
        hessian = 0, control = list(trace = FALSE))
    res = deplie2(resopt$par)
    trans = res$trans
    par.trans = res$par
    transmat = para_trans_inv(trans)
    list(mu = par.hh$mu, kappa = par.hh$kappa, prior = par.hh$prior, 
        transmat = transmat, par.trans = par.trans)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mstep of EM algorithm }
