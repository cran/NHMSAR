\name{viterbi_path}
\alias{viterbi_path}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Viterbi path homogeneous MSAR models}
\description{
 Computes the Viterbi path associated to  an homogenenous HMM or MSAR model.
}
\usage{
viterbi_path(prior, transmat, obsmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior}{ prior probabilities
PRIOR(I) = Pr(X(1) = I)
}
  \item{transmat}{ transition matrice
TRANSMAT(I,J) = Pr(X(T+1)=J | X(T)=I)
}
  \item{obsmat}{emission probabilities
OBSMAT(I,t) = Pr(Y(t) | X(t)=I)}
}


\value{List including
  \item{..$gamma}{smoothing probabilities 
  	P(X(t)|Y(0),...,Y(T))}
  \item{..$xi }{two steps smoothing probabilities 
  	 	P(X(t),X(t+1)|Y(0),...,Y(T))}
  \item{..$loglik }{log likelihood}
  \item{..$M }{Number of regimes}
  \item{..$alpha}{intermediate component in the FB algorithm (forward)}
  \item{..$beta }{intermediate component in the FB algorithm (backward)}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~}
\author{
Valerie Monbet, valerie.monbet@univ-rennes1.fr}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{forwards_backwards.R, fit.MSAR, Estep.MSAR
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Viterbi }
