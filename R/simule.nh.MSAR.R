simule.nh.MSAR <-
function(theta,Y0,T,N.samples = 1,covar.emis=NULL,covar.trans=NULL,link.ct = NULL,nc=1,S0 = NULL) {
	# If length(covar)==1, covar is built from Y with delay covar
#	if (!inherits(res, "MSAR")) 
#        stop("use only with \"MSAR\" objects")
  
if (!is.null(S0) & length(S0) != N.samples){stop("The length of S0 has to be equal to N.samples")}
M = attributes(theta)$NbRegimes
order = attributes(theta)$order
d <- attributes(theta)$NbComp
if (length(Y0[,1,1]) < order) {stop(paste("Length of Y0 should be equal to ",order,sep=""))}
label = attributes(theta)$label

L = 1
Y = array(0,c(T,N.samples,d))
S = matrix(0,T,N.samples)
Y[1:max(order,1),,] = Y0[1:max(order,1),,]
transition =  array(0,c(M,M,T,N.samples))

if (is.null(S0)){
for (ex in 1:N.samples) {
	S[1,ex] = which(rmultinom(1, size = 1, prob = theta$prior)==1)
}
} else {	  S[1,] = S0}

if (substr(label,1,1) == "N") {
	nh_transitions = attributes(theta)$nh.transitions
	if (length(covar.trans)==1) {
		L = covar.trans
	} else {
		#transition=nh_transitions(c(covar.trans),theta$par.trans,theta$transmat);
		for (ex in 1:N.samples) {transition[,,,ex]=nh_transitions(array(covar.trans[,ex,],c(T,1,dim(covar.trans)[3])),theta$par.trans,theta$transmat);}
	}
} else {
	transition = array(theta$transmat,c(M,M,T,N.samples))
}

if (max(order,L)>1) {
	for (ex in 1:N.samples) {
		for (t in 2:max(order,L)){ 
			S[t,ex] = which.max(rmultinom(1, size = 1, prob = theta$transmat[S[t-1,ex],]))
		}
	}
}

A0 = theta$A0
sigma = theta$sigma
f.emis = array(0,c(T,M,d))
if (substr(label,2,2) == "N") {
	par.emis = theta$par.emis
	nh_emissions = attributes(theta)$nh.emissions
	for (m in 1:M) {
		f.emis[,m,] = nh_emissions(covar.emis,as.matrix(par.emis[[m]])) # A voir si d>1
	}	
}

d.c = length(nc)

if(order>0){ A = theta$A}
	
if (d>1) {
		sq_sigma = list()
	for(i in 1:M){
		sq_sigma[[i]] = chol(sigma[[i]])
	}
} else
{    sq_sigma = numeric(M)
     for (i in 1:M){sq_sigma[[i]] = sqrt(sigma[[i]])}
	A = list()
	for (m in 1:M) {
		A[[m]] = list()
		for (o in 1:order) {A[[m]][[o]] = theta$A[m,o]}}
}
for (ex in 1:N.samples){

	for (t in max(c(2,order+1,L+1)):(T)){
		if (substr(label,1,1) == "N" & length(covar.trans)==1) {
			if (is.null(link.ct)) {z = Y[t-L,ex,nc,drop=FALSE]}
			else {
				z = link.ct(Y[t-L,ex,,drop=FALSE])
				#z = matrix(z,1,length(z))
			}
			transition[,,t,ex] = nh_transitions(z,theta$par.trans,theta$transmat)
		}
		S[t,ex] = which.max(rmultinom(1, size = 1, prob = transition[S[t-1,ex], ,t,ex]))
		if(order>0){
			for(o in 1:order){
				Y[t,ex,] = Y[t,ex,]+ A[[S[t,ex]]][[o]]%*%Y[t-o,ex,]
			}
		}
		Y[t,ex,] = Y[t,ex,] + A0[S[t,ex],] +  f.emis[t,S[t,ex],] + t(t(sq_sigma[[S[t,ex]]])%*%matrix(rnorm(d),d,1))
		
	}

}
return(list(S=S,Y=Y))

}
