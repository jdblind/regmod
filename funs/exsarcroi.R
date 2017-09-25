source("funs/S_rfuns.R")

extract.roi <- function(xr,pcut=.25,alpha=5,nb=26,nres=25,diro=NULL,dopas=TRUE){
#
# Arguments
# 	xr: 	input tumor file, formatted as an Amide output .tsv file.
# 			For example, xr=read.table("tumor.tsv")
#	pcut: 	set minimum fraction of values in ROI to discard
# 	alpha: 	used for setting boundary surface
#	diro:	allows to change the output directory (default output dir is 
#			./output_extractroi)
#	dopas:	if set to FALSE, deactivates PAS projection, and ROIT.txt = ROIA.txt
#			(ie xp=xa within code below)
# Outputs	
#	All outputs are created in ./output_extractroi
#	fig1.jpg:	 PAS-transformed data
#	fig2.jpg:	 Initial and projected volumes
#	fig22.jpg:	 Summary of uptake data (from msg())
#	fig3.jpg:	 Uptake profile synthesis from Gaussian fit (from fitg())
#	outzm.txt:	 tumor spine data
#	outsg.txt:	 tumor shape data
#	ROIA.txt: 	 input ROI data, possibly downsized (random sampling)
#	ROIT.txt: 	 preliminary projected ROI data (transform using arotx())
#	ROIN.txt: 	 extracted (projected) ROI data (output from roid())
#	NEWPARS.txt: output tumor characterization parameters, as follows:
# 				 newpars = c(
#				 - nlcore:		core nonlinearity
#				 - phcore/tmy:	core phase
#				 - phmn/tmy:	mean phase
#				 - disp/tmy:	spread/dispersion (boundary evaluation?)
#				 - mhet:		model heterogeneity
#				 - aphcore/tmy:	core phase weighted amplitude
#				 - aphmn/tmy:	mean phase weightded amplitude
#				 - aphmx/tmy:	75th-percentile phase boundary weighted amplitude
#				 )
#				 where tmy=max(.1e-20,mean(y,trim=.15)), and y=xr[,1]
# Note 
#	Volume sizes are identical for ROIA.txt and ROIT.txt.
# Value
#	Returns a list of uptake model parameters as follows:
#		- a: 	(core) amplitude (as a function of h only)
#		- b: 
#		- tau: 
#		- mu1: 
#		- mu2: 
#		- ...: 
#

	old.dir = getwd()
	if(length(diro)){	
		dir.create(diro,showWarnings=F)
		setwd(diro)
	} else {
		dir.create("./output_extractroi",showWarnings=F)
		setwd("./output_extractroi")
	}

	# Pre-process input VOI
	id=c(1:nrow(xr))
	ir=id
	if(nrow(xr)>100000) {  
		ir=sample(id,size=100000) 
	}
	xr = xr[ir,]
	x = xr[,3:5]
	# x=cbind(xr$V3,xr$V4,xr$V5) 
	y=xr[,1]
	n=length(y)
	yl=sort(y)[.1*n]
	ym=sort(y[y>yl])[.5*n]
	nn=min(1000,.1*n) 
	my=median(sort(y)[c((n-nn):n)])
	thres=max(.1*my,sort(y[y>yl])[pcut*n])
	id=c(1:n)
	ir=id
	if(n>85000){  
		ir=sample(id,size=85000) 
	}
	xvals=cbind(x,y)[ir,]
	nrr=nrow(xvals)
	write(file="thres.txt",thres)								# uptake threshold
	write(file="ROIA.txt",t(cbind(rep(1,nrr),xvals)),ncol=5)	# ROI DATA
	d3=max(diff(as.real(names(table(x[,3])))))
	d2=max(diff(as.real(names(table(x[,2])))))
	d1=max(diff(as.real(names(table(x[,1])))))
	dd=min(d1,d2,d3)*1.2
	
	# PAS Transform Data
	if(dopas){
		out.arot=arotx(x,y,thres) 
	} else { # make up dummy out.arot
		arotc = c(0,0,0)
		arots = c(1,1,1)
		arotGxi = diag(1,3)
		arotS = diag(1,3)
		out.arot=list(amide=cbind(arotc,t(arotGxi),arots),c=arotc,s=arots,Gxi=arotGxi,S=arotS)
	}
	out=out.arot$amide
	mu=out[,1]
	gam=t(out[,2:4]) 
	rrval=max(out[,5])/min(out[,5])
	xt = cbind(x[,1]-mu[1],x[,2]-mu[2],x[,3]-mu[3])%*%gam
	write(file="ROIT.txt",t(cbind(rep(1,nrow(xt)),y,xt)),ncol=5) # RAW PAS ROI DATA
	if(nb<=0){ ## added this 18/03/2014 as nb=Jh in regularized model
		nb = min(72,round((max(xt[,3])-min(xt[,3]))/dd)) 
	}
	# rrval=max(rrval/1.5,2)
		
	# Compute Spine, Volume, spine, sg scaling and Boundary Surface
	jpeg("fig1.jpg",quality=95)
	par(mfcol=c(2,2),pty="m")
		out=msg(xt,y,thres,alpha,dd,nb,nres)
	dev.off()
	
	bxt = bpts(out$z,out$m,out$th,out$sg)
	nxt = npts(out$z,out$m,out$th,out$sg) 	  # Surface and Normals - PAS domain
	bxo = cbind(bxt[,1],bxt[,2],bxt[,3])%*%t(gam) 
	bxo = cbind(bxo[,1],bxo[,2],1.4*bxo[,3])  # Surface             - Original domain
	nxo = cbind(nxt[,1],nxt[,2],nxt[,3])%*%t(gam) 
	nxo = cbind(nxo[,1],nxo[,2],1.4*nxo[,3])  # Normals             - Original domain
	
	write(file="outzm.txt",t(cbind(out$z,out$m)),ncol=3)			      # Tumor Spine
	write(file="outsg.txt",t(rbind(out$th,out$sg)),ncol=(1+nrow(out$sg))) # Tumor Shape
	
	# Obtain the extracted ROI data and write to a file 
	# [easily extended for boundary data selection]
	o=roid(xt,y,out$z,out$m,out$sg)
	o.roid=o
	nr=nrow(o$roi)
	zuy=o$zuy  
	write(file="ROIN.txt",t(cbind(rep(1,nr),o$roi)),ncol=5)			# ROI DATA
	
	# Visualize Boundary Surface
	jpeg("fig2.jpg",quality=95)
		par(mfrow=c(2,3))
		xv=c(0,1,0);rrx=1.1
		viewxv2(xv,bxt,bxo,nxo,"x3","x1","Original",rrx)
		xv=c(1,0,0);rrx=1.1
		viewxv2(xv,bxt,bxo,nxo,"x3","x2"," ",rrx)
		xv=c(0,0,1);rrx=1.1
		viewxv2(xv,bxt,bxo,nxo,"x2","x1"," ",rrx) 
		#
		xv=c(0,1,0);rrx=1.1
		viewxv(xv,cbind(bxt[,1],bxt[,2],bxt[,3]),cbind(nxt[,1],nxt[,2],nxt[,3]),
				"x3'","x1'","Transformed (PAS)",rrx) 
		xv=c(1,0,0);rrx=1.1
		viewxv(xv,cbind(bxt[,1],bxt[,2],bxt[,3]),cbind(nxt[,1],nxt[,2],nxt[,3]),
				"x3'","x2'"," ",rrx)
		xv=c(0,0,1);rrx=1.1
		viewxv(xv,cbind(bxt[,1],bxt[,2],bxt[,3]),cbind(nxt[,1],nxt[,2],nxt[,3]),
				"x2'","x1'"," ",rrx)
	dev.off()
	
	jpeg("fig22.jpg",quality=95) # NB: out = output of msg()
		par(mfrow=c(2,2),pty="m")
		hist(y,xlab="Uptake",main=" ",axes=F)
		axis(1)
		abline(v=thres,lwd=c(4),col=c(8)) 
		par(new=T)
		qqnorm(y,axes=F,datax=T,xlab=" ",ylab=" ",main=" ")
		#
		plot(out$m[,1],out$z,main=" ",axes=F,xlab="Center x',y' (mm)",
				ylab="z' (mm)",lwd=c(4),type="l",xlim=3*range(out$m),col=c(8))
		lines(out$m[,2],out$z,lwd=c(4),col=c(1))
		axis(1)
		axis(2)
		#
		plot(out$vol*(d1*d2*d3)/10^3,out$z,main=" ",axes=F,xlab="Volume",
				ylab="z' (mm)",lwd=c(3),type="l",col=c(1))
		axis(1)
		axis(2)
		#
		sg=out$sg
		sgn=t(out$sg) 
		# Surface Sheet
		ss=sgn
		xrm=(out$th-pi)/(pi)
		for(k in c(1:ncol(sgn))) {
			xrk = (sqrt(out$vol[k]/max(out$vol)))*xrm #(out$th-pi)/(pi)
			ss[,k]=approx(xrk,sgn[,k],xo=xrm,rule=1)$y 
		}
		image(xrm,1.*out$z,ss,axes=F,main=" ",xlim=c(-1,1)*1.5,ylim=1.2*range(1.*c(bxt)),
				xlab="Topography",ylab="",col=grey(seq(0,1,length=256))) # Meridian at x'=0 
	dev.off()
	
	# Synthesis of zuy data using Gaussian Model - Determine Phase Pattern of the Tumor
	nz=nrow(sg)
	nu=nres
	e=fitg(zuy,nz,nu)
	nz=e$nz
		
	jpeg("fig3.jpg",quality=95)
		par(mfrow=c(2,2),pty="m")
		#
		# (a)
		imn=e$ims*0
		imn=rbind(rep(0,nz),imn)
		rad=apply(sg,1,mean) 
		nz=length(rad)
		ff=(max(e$pars[,1])-min(e$pars[,1]))/max(rad)
		ff=2.
		for(k in 1:nz) { 
			sfac=rad[k]/max(rad)
			uk=c(1:(nu+1))*sfac
		 	imn[,k]=approx(uk,e$ims[,k],xout=c(0:(nu+1))*ff,rule=2,yright=0)$y 
		}
		image(c(0:(nu+1))*ff,e$pars[,1],-imn,col=grey((0:255)/256),
				ylab="z (mm)",xlab=" ",axes=F,main="Radial Uptake")
		lines((nu+1)*rad/max(rad),e$pars[,1],lty=c(3))
		axis(2)
		axis(1,at=c(1,(nu)),labels=c("C","B"),cex=.2)
		#
		# (b)
		v=e$fit[,1]
		vm=max(v)
		vm=1
		v=v/vm 
		yn=e$fit[,3] 
		wy=e$fit[,2]
		wy=wy/max(wy) 
		wcut=sort(wy)[.4*length(wy)]
		if(wcut==max(wy)){ 
			wcut=sort(wy)[.2*length(wy)]
		}
		plot(v[wy>wcut],yn[wy>wcut],pch=".",main="Profile",
				axes=F,xlab="[<- Core]     Phase    [Boundary ->]",ylab="Uptake")
		axis(1)
		lines(unismooth2(v[wy>wcut],yn[wy>wcut]),lwd=c(2),col=c(3))
		#lines(unismooth(v[wy>wcut],yn[wy>wcut],wy[wy>wcut]),lwd=c(2),col=c(8))
		#
		# (c)
		o=smooth.spline(e$pars[,1],e$pars[,8]/vm)
		plot(e$pars[,4]/vm,e$pars[,1],pch="*",type="n",	
				main="Phase",axes=F,xlab=" ",ylab="z (mm)",
				xlim=range(v[wy>wcut],o$y))
		lines(o$y,o$x,lwd=c(2),col=c(2))
		o=smooth.spline(e$pars[,1],e$pars[,4]/vm)
		lines(o$y,o$x,lwd=c(2))
		points(e$pars[,5]/vm,e$pars[,1],pch=".",type="n",col=c(4))
		points(e$pars[,6]/vm,e$pars[,1],pch=".",type="n",col=c(4))
		ol=smooth.spline(e$pars[,1],e$pars[,5]/vm)
		ou=smooth.spline(e$pars[,1],e$pars[,6]/vm)
		lines(ol$y,ol$x,pch="=",col=c(4))
		lines(ou$y,ou$x,pch="=",col=c(4))
		axis(1); axis(2)
		#
		# (d)
		plot(100-100*e$pars[,3],e$pars[,1],xlim=c(0,100),
				main="Heterogeneity",axes=F,xlab="1-R^2",ylab="z (mm)")
		axis(1); axis(2) 
	dev.off()
	
	# NEW Tumor Characterization parameters
	# Core Non-linearity (mm units) 
	lmo=lm(out$m[,1]~out$z,weights=e$pars[,10]) 
	fit1=lmo$fitted
	wrss1=sum(e$pars[,10]*(out$m[,1]-fit1)^2)/sum(e$pars[,10]*out$m[,1]^2)
	lmo=lm(out$m[,2]~out$z,weights=e$pars[,10]) 
	fit2=lmo$fitted
	wrss2=sum(e$pars[,10]*(out$m[,2]-fit2)^2)/sum(e$pars[,10])
	nlcore=sqrt(wrss1+wrss2)
	# Core Phase 
	phcore=sum(e$pars[,10]*e$pars[,7]*e$pars[,8])/sum(e$pars[,10])
	# Mean Phase 
	phmn=sum(e$pars[,10]*e$pars[,7]*e$pars[,4])/sum(e$pars[,10])
	# Spread [Boundary Evaluation] [dispersion]
	disp=sum(e$pars[,10]*e$pars[,7]*e$pars[,9])/sum(e$pars[,10])	
	# Model Heterogeneity
	mhet=sum(e$pars[,10]*(1-e$pars[,3])*100)/sum(e$pars[,10])  
	# e$pars[,3] is a model R-squared by pas slice.
	# Core Phase weighted Amplitude 
	aphcore=sum(e$pars[,10]*e$pars[,7]*exp(-.5*e$pars[,8]^2) )/sum(e$pars[,10])
	# Mean Phase weighted Amplitude 
	aphmn=sum(e$pars[,10]*e$pars[,7]*exp(-.5*e$pars[,4]^2)   )/sum(e$pars[,10])
	# 75 Percentile Phase [Boundary Evaluation] weighted Amplitude 
	aphmx = sum(e$pars[,10]*e$pars[,7]*exp(-.5*e$pars[,6]^2))/sum(e$pars[,10])	
	tmy=max(.1e-20,mean(y,trim=.15))	
	newpars = c(nlcore,phcore/tmy,phmn/tmy,disp/tmy,mhet,aphcore/tmy,aphmn/tmy,aphmx/tmy)
	write(file="NEWPARS.txt",newpars,ncol=8)
	setwd(old.dir)

	# output model parameters:
	#	- phi and h determine the product-basis for the evaluation splines 
	#	  (of lengths Jh and Jphi)
	#	- (c,s,xi,mu) are structural parameters
	# 	- (a,tau,b) are uptake profile parameters
	# Note: a is a(h), must be adapted outside to a(h,phi)
	# Jh = nrow(out$sg)
	Jh = out$nb
	Jphi = ncol(out$sg)
	lpars = list(msg.out = out,
				z.pasbins=out$z.pasbins,
				z.midpoints=out$z.midpoints,
				phi = out$th,
				h = out$z,
				Jh = Jh,
				Jphi = Jphi,
				Gxi = out.arot$Gxi,
 				c = out.arot$c,
				s = out.arot$s,
				xi = euler.angles(out.arot$Gxi,method=0),
				mu1 = out$m[,1],  # out = output of msg()
				mu2 = out$m[,2],
				sg = out$sg,
				tau.old = e$pars[,8]*e$pars[,9],
				a.old = e$pars[,7],
				o.ccore = out$ccore,
				### new outputs...
				efit = e$fit,
				epars = e$pars,
				roid = o.roid$roi,
				nz = e$nz,
				nk = e$nk,
				zuy = o.roid$zuy,
				th = o.roid$th,
				r = o.roid$r,
				u = o.roid$u,
				icut = o.roid$icut,
				a = e$ak,
				b = e$bk, 
				tau = e$tauk,
				muk = e$muk,
				sdk = e$sdk,
				vox.phase = e$vox.phase,
				phase = e$phase,
				bscale = e$pars[,9],
			 	# add-ons:
			 	gsg.ee=out$gsg.ee,
			 	gsg.sg0=out$gsg.sg0,
			 	S=out$S,
			 	lams=out$lams,
				ccore.mk = out$ccore[,3:4],
				ccore.sk = out$ccore[,5:8],
				gsg.th = out$th,
				gsg.muk = out$m,
				gsg.sg = out$sg,
				siv=out$siv,
				res=out$res,
				u.from.gsg=out$u,
				uk=e$uk,
				re=out$re,
				uu=out$uu # actually = v[siv]
				)
	return(lpars)
}
