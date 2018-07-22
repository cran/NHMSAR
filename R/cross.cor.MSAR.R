cross.cor.MSAR <-
function(data,X=NULL,nc1=1,nc2=2,lag=10,regime=0,CI = FALSE,Bsim=0,N.samples=1,add=FALSE,col=1,names = NULL,alpha=.1,ylab="Cross-Correlation",dt=1,ylim=c(-.1,1)) {
	T = dim(data)[1]
	N.samp = dim(data)[2]
	if (is.null(X)) {X = array(1,c(T,N.samp,1))}
	m1 = mean(data[,,nc1][X!=regime],na.rm=TRUE)
	sd1 = sd(data[,,nc1][X!=regime],na.rm=TRUE)
	m2 = mean(data[,,nc2][X!=regime],na.rm=TRUE)
	sd2 = sd(data[,,nc2][X!=regime],na.rm=TRUE)
	s.cc = matrix(0,1,length(-lag:lag))
	n.obs = matrix(0,1,length(-lag:lag))	
	T = dim(data)[1]
	if (CI==TRUE) {
		S.CC = matrix(0,Bsim+1,2*lag+1)
		N.OBS = matrix(0,Bsim+1,2*lag+1)	
	}
	cc.ex = matrix(0,dim(data)[2],2*lag+1)
	for (ex in 1:dim(data)[2]) {
			x1 = data[,ex,nc1]
			x2 = data[,ex,nc2]
			if (regime > 0) {
				x1[X[,ex,]==regime] = NA
			    x2[X[,ex,]==regime] = NA
			}
			cnt = 0 ; cc = NULL
			for (h in (lag+1):1) {
				cnt = cnt+1
				cc[cnt] = cor(x1[1:(T-h+1)],x2[h:T],use="pairwise.complete.obs")
				s.cc[cnt] = s.cc[cnt]+sum((x1[1:(T-h+1)]-m1)*(x2[h:T]-m2),na.rm=T)
				n.obs[cnt] = n.obs[cnt]+sum(!is.na((x1[1:(T-h+1)]-m1)*(x2[h:T]-m2)))
			}
			for (h in 2:(lag+1)) {
				cnt = cnt+1
				cc[cnt] = cor(x1[h:T],x2[1:(T-h+1)],use="pairwise.complete.obs")
				s.cc[cnt] = s.cc[cnt]+sum((x1[h:T]-m1)*(x2[1:(T-h+1)]-m2),na.rm=T)
				n.obs[cnt] = n.obs[cnt]+sum(!is.na((x1[h:T]-m1)*(x2[1:(T-h+1)]-m2)))
			}
			if (CI==TRUE) {if (length(which(((1:Bsim)*N.samples)==ex))>0) {
				S.CC[ex/N.samples+1,] = s.cc
				N.OBS[ex/N.samples+1,] = n.obs
			}}
			cc.ex[ex,] = cc 
	}
	if (add==FALSE) {
		plot((-lag:lag)*dt,s.cc/n.obs/sd1/sd2,typ="l",ylab=ylab,xlab="Lag (day)",ylim=ylim,col=col,lwd=2)
		grid()
		if (!is.null(names)){title(paste(names[nc1],", ",names[nc2],sep=""))}
	} else {
	  lines((-lag:lag)*dt,s.cc/n.obs/sd1/sd2,col=col)
	  #matlines((matrix((-lag:lag)*dt,2*lag+1,dim(cc.ex)[1])),t(cc.ex/sd1/sd2),col="gray",lty=1)
	}
    if (CI==TRUE) {
			IC = matrix(0,2*lag+1,2)
        	for (k in 1:(2*lag+1)) {
        		tmp.cc = diff(S.CC[1:(Bsim+1),k])
        		tmp.obs = diff(N.OBS[1:(Bsim+1),k])
        		IC[k,] = quantile(tmp.cc/tmp.obs,probs=c(alpha/2,1-alpha/2)) 
        	}
        	lines((-lag:lag)*dt,IC[,1]/sd1/sd2,lty=2,col=col)
        	lines((-lag:lag)*dt,IC[,2]/sd1/sd2,lty=2,col=col)
        } else {IC = NULL}
     list(ccf=s.cc/n.obs/sd1/sd2,lag=-lag:lag,CI=IC/sd1/sd2)
}
