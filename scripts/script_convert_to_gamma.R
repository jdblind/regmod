# Reconstruct estimated profiles from fitg() 
v = eroi$efit[,1] 
zuy=eroi$zuy
zz=zuy[,1]   # h-range values
uu=zuy[,2]   # voxel phases
yy=zuy[,3]   # uptakes
nz=Jh
nn=length(yy) 
nv=round(max(1,nn/nz)) 
# zb=sort(zz)[c(1,c(1:(nz-1))*nv,nn)] 			# discretized z-grid
# zb[1]=zb[1]-max(.1e-9,.001*(zb[2]-zb[1]))   # so as to include boundaries
zb=pasbins(zz,nz,nv,nn)
ghat = zhat = uhat = uout = vtst = vout = numeric(nn)
ghat2 = zhat2 = numeric(nn)
zza = xp[,4]
for(k in 1:nz) { # yu is NOT SCALED HERE
	kinds = ((zz<=zb[k+1])&(zz>zb[k]))
	yu=yy[kinds]
	ghat[kinds] = yu/eroi$a[k]
	vk=v[kinds]
	ghat2[kinds] = exp(-vk^2/2)
}
# from Gaussian to Gamma profile
g2g = gauss.to.gam(v,alpha,doplot=FALSE,gquant=.992)
fac = g2g$fac

if(0){
	# checks...
	pdf("out/conversion_process.pdf")
	par(mfrow=c(1,2),font=2,font.lab=2,font.axis=2)
	# reproduce fig3.pdf
	yn = ghat
	wy=eroi$efit[,2]
	wy=wy/max(wy) 
	wcut=sort(wy)[.4*length(wy)]
	if(wcut==max(wy)){ 
		wcut=sort(wy)[.2*length(wy)]
	}
	plot(v[wy>wcut],yn[wy>wcut],pch=".",main="Profile",
			axes=F,xlab="[<- Core]     Phase    [Boundary ->]",ylab="Uptake")
	axis(1); axis(2)
	i.s = which(wy>wcut)
	v.s = v[i.s]; y.s = yn[i.s]
	i.s = which(!is.na(y.s))
	y.s = y.s[i.s]; v.s = v.s[i.s]
	lines(unismooth2(v.s,y.s),lwd=c(2),col=c(3))
	#####
	i.s = which(!is.na(ghat))
	ylimi = range(c(ghat[i.s]/g2g$fac,g2g$gam.fit))
	plot(v,ghat/g2g$fac,ylim=ylimi,main="Gaussian to Gamma fit",
		xlab="phases",ylab="uptake",axes=F)
	is = order(v)
	points(v[is],ghat2[is]/g2g$fac,col="navy",t='l',lwd=6)
	abline(v=v[which.max(ghat2)],col='navy')
	par(new=TRUE)
	is = order(g2g$gv)
	plot(g2g$gv[is],g2g$gam.fit[is],t='l',lwd=6,
		col='green',axes=F,ylim=ylimi,xlab="",ylab="")
	axis(1)
	points(g2g$gv,ghat/g2g$fac,col='gray',pch=20,cex=.6)
	par(new=TRUE)
	plot(v[is],ghat2[is]/g2g$fac,col="navy",t='l',lwd=6,xlab="",ylab="",
		ylim=range(c(ghat[i.s]/g2g$fac,g2g$gam.fit)),axes=F)
	par(new=TRUE)
	plot(g2g$gv[is],g2g$gam.fit[is],t='l',lwd=6,
		col='green',ylim=ylimi,xlab="",ylab="",axes=F)
	axis(2)
	abline(h=0,lwd=1)
	legend("topright",col=c("navy","green"),lwd=c(6,6),
		legend=c("Gauss","Gamma"))
	abline(v=g2g$gv[which.max(g2g$gam.fit)],col='green')
	#
	range(g2g$gv)
	v[which.max(ghat)]
	g2g$gv[which.max(g2g$gam.fit)]
	dev.off()
}
