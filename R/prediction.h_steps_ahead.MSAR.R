prediction.h_steps_ahead.MSAR <-
function(data,theta,horizon=1,covar.emis=NULL,covar.trans=NULL,ex=1){
	if (!inherits(theta, "MSAR")) 
        stop("use only with \"MSAR\" objects")
	label = attributes(theta)$label
	M = attributes(theta)$NbRegimes
	order = attributes(theta)$order
	T = dim(data)[1]
	N.samples = dim(data)[2]
	d = dim(data)[3]
	if (d>1) {A=theta$A
	} else {
		A = list()
		for (m in 1:M) {
			A[[m]] = list()
			for (o in 1:order) { 
				A[[m]][[o]] = theta$A[m,o]
			}
		}
	}
	y.p = array(0,c(horizon,length(ex),d))
	var.p = array(0,c(horizon,length(ex),d,d))
	Lag = 0
    if (substr(label,1,1) == "N" & length(covar.trans)==1) {
    	Lag = covar.trans+1
    	covar.trans = array(data[(1):(T-Lag+1),,],c(T-Lag+1,dim(data)[2],d))
    	data =  array(data[Lag:T,,],c(T-Lag+1,dim(data)[2],d)) 
    	T=T-1
    }	
	#browser()
	pr = array(0,c(T+horizon,M,N.samples))
	for (iex in 1:length(ex)) {
		g <- emisprob.MSAR(data[,ex[iex],],theta=theta,covar=covar.emis[,ex[iex],])
		if (substr(label,1,1) == "H") {
			FB = forwards_backwards(theta$prior, theta$transmat, g)
			transitions = array(theta$transmat,c(M,M,T))
			alpha = matrix(0,M,T+horizon)
      alpha[,(order+1):T] = FB$alpha
		} else {
    		ncov.trans = dim(covar.trans)[3 ]
    		if(is.null(ncov.trans) || is.na(ncov.trans)){ncov.trans=1}
			  #covar.trans=array(covar.trans,c(T-L,N.samples,ncov.trans))
    		ct = array(0,c(T,1,ncov.trans))
    		ct[1:dim(covar.trans)[1],,] = covar.trans[,ex[iex],]
			  transmat = theta$transmat
    		par.trans = theta$par.trans
    		nh_transition = attributes(theta)$nh.transitions
    		inp = ct[(order+1):T,1,] #see line 331
    		transitions = nh_transition(array(inp,c(length(inp)/ncov.trans,1,ncov.trans)),par.trans,transmat)
        FB = nhforwards_backwards(theta$prior, transitions, g)
        alpha = matrix(0,M,T+horizon)
        alpha[,(order+1):T] = FB$alpha
        ct = array(0,c(M,M,T))
        ct[,,(order+1):T] = transitions
        transitions = ct 
		}
		data_prev = data[(T-o+1):(T),ex[iex],] 
		for (t in (T):(T+horizon)) {
			if (t==T){pr[t,,iex] = t(matrix(alpha[,t-1]))%*%transitions[,,t] # P(S_t = s|y_{1:(t-1)})
			} else {pr[t,,iex] = transitions[,,t]%*%matrix(pr[t-1,,iex],M,1)}
			y.p[t-T+1,iex, ] = 0
			var.p[t-T+1,iex,,] =  0*theta$sigma[[1]]# same size as sigma, case d=1 ?
			for (m in 1:M) {
				var.p[t-T+1,iex,,] = var.p[t-T+1,iex,,]+pr[t,m,iex]*theta$sigma[[m]]
				A0.m = theta$A0[m,]
				if (substr(label,2,2)=="N") {
					A0.m = A0.m+attributes(theta)$nh.emissions(matrix(covar.emis[t,ex[iex],],1,dim(covar.emis)[3]),as.matrix(theta$par.emis[[m]]))
				} 
				y.hat = A0.m
				if (order>0) {
					for (o in 1:order) {
						y.hat = y.hat+A[[m]][[o]]%*%data_prev									}
				}
				y.p[t-T+1,iex,] = y.p[t-T+1,iex,]+pr[t,m,iex]*y.hat 
			}
			data_prev[1:(o-1),,] = data_prev[2:o,,]
			data_prev[o,,] = y.p[t-T+1,iex,]
		}	
	}
	return(list(y.p=y.p,var.p=var.p,pr=pr))
}
