viterbi_path = function(prior, transmat, obsmat){
# VITERBI Find the most-probable (Viterbi) path through the HMM state trellis.
# path = viterbi(prior, transmat, obsmat)
#
# Inputs:
  # prior(i) = Pr(Q(1) = i)
  # transmat(i,j) = Pr(Q(t+1)=j | Q(t)=i)
  # obsmat(i,t) = Pr(y(t) | Q(t)=i) emission probabilities
  #
  # Outputs:
    # path(t) = q(t), where q1 ... qT is the argmax of the above expression.
  
  # delta(j,t) = prob. of the best sequence of length t-1 and then going to state j, and O(1:t)
  # psi(j,t) = the best predecessor state, given that we ended up in state j at t
  
  scaled = 1
  
  T = dim(obsmat)[2]
  Q = length(prior)
  
  delta = matrix(0,Q,T)
  psi = matrix(0,Q,T)
  path = matrix(0,1,T)
  scale = matrix(1,1,T)
  
  t=1;
  delta[,t] = prior * obsmat[,t]
  if (scaled){
    scale[t] = 1/sum(delta[,t]) # ou sum(delta[,t])
    delta[,t] = normalise(delta[,t])
  }
  psi[,t] = 0; # arbitrary value, since there is no predecessor to t=1
  for (t in 2:T){
    for (j in 1:Q){
      delta[j,t] = max(delta[,t-1] * transmat[,j])
      psi[j,t] = which.max(delta[,t-1] * transmat[,j])
      delta[j,t] = delta[j,t] * obsmat[j,t]
    }
    if (scaled){
      scale[t] = 1/sum(delta[,t])
      delta[,t] = normalise(delta[,t])
    }
  }

  path[T] = which.max(delta[,T])
  for (t in seq(T-1,1,-1)){
     path[t] = psi[path[t+1],t+1]
  }
  
  
  return(path)
}