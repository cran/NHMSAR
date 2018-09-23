valid_all.MSAR <-
function(data,simu,title="",id=1,alpha=.05,spaghetti=TRUE,mfrow=NULL,save=FALSE,output=FALSE,root.filename=" ",path=NULL,col="red",width=4,height=4){

if (is.null(mfrow)) {dev.new()} else {par(mfrow=mfrow)}
T = dim(data)[1]
N.samples = dim(data)[2]
N.sim = dim(simu)[2]
Bsim = N.sim/N.samples
qqp = qqplot(data[,,id],simu[,,id],pch=20,xlab='Observations',ylab='Simulations',cex=.6)
title(title)
abline(a=0,b=1)
q = matrix(0,Bsim,length(data[,,id]))
s.data = sort(data[,,id])
lens = length(s.data)
for (k in 1:Bsim) {
  tmp = 
	q[k,] = sort(simu[,((k-1)*N.samples+1):(k*N.samples),id])
}
if (spaghetti){
	for (k in 1:Bsim){lines(s.data,approx(1:length(data[,,1]),q[k,], n = lens)$y,col="gray")}
} else {
	lines(s.data,approx(1:length(data[,,1]),q[,k], n = lens)$y,col="gray")
	IC = matrix(0,length(data[,,id]),2)
	for (k in 1:length(data[,,id])) {
		IC[k,] = quantile(q[,k],probs=c(alpha/2,1-alpha/2))
		}
		lines(s.data,approx(1:length(data[,,1]),IC[,1], n = lens)$y,lty=2,lwd=1.5)
		lines(s.data,approx(1:length(data[,,1]),IC[,2], n = lens)$y,lty=2,lwd=1.5)
}

# Correlations
if (save) {
	filename = paste(path,"qqplot-",root.filename,".eps",sep="")
	dev.copy2eps(file=filename,width=width,height=height)
}
if (is.null(mfrow)) {dev.new()}
C = cor.MSAR(array(data[,,id],c(T,N.samples,1)),array(simu[,,id],c(T,N.sim,1)),lag=15)
plot(0:14,C$C.data,typ="l",ylab="Correlation",xlab="Time (days)",lwd=2)
title(title)
if (spaghetti){
	for (k in 1:Bsim){lines(0:14,C$C.sim[,k],col="gray")}
	lines(0:14,C$C.data,lwd=2)
} else {
	lines(0:14,C$C.sim,col=col,lwd=1.5)
	lines(0:14,C$CI.sim[1,],col=col,lty=3,lwd=1.5)
	lines(0:14,C$CI.sim[2,],col=col,lty=3,lwd=1.5)
}
if (save) {
	filename = paste(path,"Cor-",root.filename,".eps",sep="")
	dev.copy2eps(file=filename,width=width,height=height)
}

# Upcrossings
if (is.null(mfrow)) {dev.new()}
u = seq(min(data[,,id],na.rm=TRUE),max(data[,,id],na.rm=TRUE),length.out=20)
gr.d = ENu_graph(data[,,id],u)
k=1
u = seq(min(simu[,,id]),max(simu[,((k-1)*N.samples+1):(k*N.samples),id]),length.out=50)
gr = ENu_graph(simu[,((k-1)*N.samples+1):(k*N.samples),id],u,add=TRUE,col="gray",CI = FALSE,N.s.data=dim(data)[2])
gr.d = ENu_graph(data[,,id],u,ylim=c(range(gr.d$Nu,gr$Nu)))	
if (spaghetti){
	for (k in 1:Bsim){
		u = seq(min(simu[,,id]),max(simu[,((k-1)*N.samples+1):(k*N.samples),id]),length.out=50)
		gr = ENu_graph(simu[,((k-1)*N.samples+1):(k*N.samples),id],u,add=TRUE,col="gray",CI = FALSE,N.s.data=dim(data)[2])
	}
	lines(c(gr.d$F,1),c(gr.d$Nu,0),lwd=3)
} else {		
	u = seq(min(simu[,,id]),max(simu[,,id]),length.out=50)
	gr = ENu_graph(simu[,,id],u,add=TRUE,col=col,CI = TRUE,N.s.data=dim(data)[2])
	if (max(gr$Nu)>max(gr.d$Nu)) {
		plot(c(gr$F,1),c(gr$CI[,2],0),lty=3,col=col,typ="l",xlab="P(Y<u)",ylab="Intensity of upcrossings")
		lines(c(gr$F,1),c(gr$CI[,1],0),lty=3,col=col)
		lines(c(gr$F,1),c(gr$Nu,0),typ="l",col=col,xlab="P(Y<u)",ylab="Intensity of upcrossings")
		lines(c(gr.d$F,1),c(gr.d$Nu,0),lwd=3/2)
}
}
abline(v=.5,lty=2)
title(title)
if (save) {
	filename = paste(path,"ENu-",root.filename,".eps",sep="")
	dev.copy2eps(file=filename,width=width,height=height)
}
if (is.null(mfrow)) {dev.new()}
if (spaghetti){
	MDO = MeanDurOver(array(data[,,id],c(T,N.samples,1)),array(simu[,((k-1)*N.samples+1):(k*N.samples),id],c(T,N.samples,1)),u,plot=FALSE)
	plot(MDO$F,MDO$mdo.data,typ="l",ylim=range(c(MDO$mdo.data,MDO$mdo.sim)),log="y",lwd=3,xlab = "P(Y<u)",ylab="Mean sojourn duration over u (log scale)")
	lines(MDO$F.sim,MDO$mdo.sim,col="gray")
	for (k in 2:Bsim){
		MDO = MeanDurOver(array(data[,,id],c(T,N.samples,1)),array(simu[,((k-1)*N.samples+1):(k*N.samples),id],c(T,N.samples,1)),u,plot=FALSE)
		lines(MDO$F.sim,MDO$mdo.sim,col="gray")
	}
	lines(MDO$F,MDO$mdo.data,lwd=3)
} else { 
	MDO = MeanDurOver(array(data[,,id],c(T,N.samples,1)),array(simu[,,id],c(T,N.sim,1)),u,col=col)
}
title(title)
if (save) {
	filename = paste(path,"MeanDurOver-",root.filename,".eps",sep="")
	dev.copy2eps(file=filename,width=width,height=height)
}
if (is.null(mfrow)) {dev.new()}
if (spaghetti){
	MDU = MeanDurUnder(array(data[,,id],c(T,N.samples,1)),array(simu[,((k-1)*N.samples+1):(k*N.samples),id],c(T,N.samples,1)),u,plot=FALSE)
	plot(MDU$F,MDU$mdu.data,typ="l",ylim=range(c(MDU$mdu.data,MDU$mdu.sim)),log="y",lwd=3,xlab = "P(Y<u)",ylab="Mean sojourn duration over u (log scale)")
	lines(MDU$F.sim,MDU$mdu.sim,col="gray")
	for (k in 2:Bsim){
		MDU = MeanDurUnder(array(data[,,id],c(T,N.samples,1)),array(simu[,((k-1)*N.samples+1):(k*N.samples),id],c(T,N.samples,1)),u,plot=FALSE)
		lines(MDU$F.sim,MDU$mdu.sim,col="gray")
	}
	lines(MDU$F,MDU$mdu.data,lwd=3)
} else { 
	MDU = MeanDurUnder(array(data[,,id],c(T,N.samples,1)),array(simu[,,id],c(T,N.sim,1)),u,col=col) 
}
title(title)
if (save) {
	filename = paste(path,"MeanDurUnder-",root.filename,".eps",sep="")
	dev.copy2eps(file=filename,width=width,height=height)
}
if (output) {return(list(qqp=qqp,C=C,ENu.data=gr.d,ENu.simu=gr,MDO=MDO,MDU=MDU ))}
}
