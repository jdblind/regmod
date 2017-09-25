# --- work out coverage ---
ab = hab(xt)
hs = seq(ab[1],ab[2],len=(Jh+1))
phis = seq(phia,phib,len=(Jphi+1))
# work out sectors
n = length(z)
hk = numeric(Jh)
ds = matrix(0,nr=Jh,nc=n)
for(k in 1:Jh){
	hk[k] = sum(hs[k:(k+1)])/2 # center points
	ds[k,] =((xt[,3]-hk[k])^2) # distances to center points
}
bins = apply(ds,2,which.min)	
phik = numeric(Jphi)
dps = matrix(0,nr=Jphi,nc=n)
for(k in 1:Jphi){
	phik[k] = sum(phis[k:(k+1)])/2 # center points
	dps[k,] = ((rphs[,2]-phik[k])^2)
}
pbins = apply(dps,2,which.min)
# sample sizes per (phi,h)
tsizes = numeric(Jh)
sns = matrix(0,nc=Jphi,nr=Jh)
for(k in 1:Jh){
	ks = which(bins==k)
	tsizes[k] = length(ks)
	for(j in 1:Jphi){
		js = which(pbins[ks]==j)
		sns[k,j] = length(js)
	}
}
# --------------
nb = length(eroi$z.midpoints)
n = nrow(xt)
nv = round(max(1,n/nb)) 
zb = eroi$z.pasbins
zz = xt[,3]
tsizes2 = numeric(Jh)
for(k in 1:Jh){
	kinds = which((zz<=zb[k+1])&(zz>zb[k]))
	tsizes2[k] = length(kinds)
	
}
# --------------
summary(tsizes-tsizes2)

if(0){
	if(doquartz){quartz()}
	par(mfrow=c(2,2),pch=20)
	# sampling coverage
	plot(rphs[,2:3],xlab="phi",ylab="h",main=paste("Coverage",ptid))
	#
	plot(tsizes,t="h",lwd=7,main="transverse slices sizes",
		col='navy',xlab="slice",ylab="sample size") 
	#
	image(t(sns),col=gray(c(0:255)/255),axes=F,
		xlab="sector",ylab="h-elevation",main="sector sizes")
	axis(1,at=c(0,1),labels=c(1,Jphi))
	axis(2,at=c(0,1),labels=c(1,Jh))
	#
	plot(tsizes2,t="h",lwd=7,main="transverse slices sizes (2)",
		col='navy',xlab="slice",ylab="sample size") 

	if(doquartz){quartz()}
	par(mfrow=c(min(ceiling(Jh/5),5),5))
	rr = c(0,max(sns))
	for(k in 1:Jh){
		plot(sns[k,],t="h",lwd=3,col='navy',ylim=rr,
			xlab="sectors",ylab="sample size",main=paste("(",ptid,") slice ",k,sep=""))
		abline(h=seq(5,rr[2],by=5),col=8,lwd=.5,lty=1)
		points(sns[k,],t="h",lwd=6,col='navy')
	}
}