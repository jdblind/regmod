prune.xv <- function(bgam,ngam){
	hs = unique(bgam[,3])
	LH = length(hs)
	xl = yl = numeric(LH)
	bgam0 = cbind(bgam, numeric(LH), numeric(LH), numeric(LH))
	for(hi in 1:LH){
		slki = which(bgam[,3]==hs[hi])
		slk = bgam[slki,]
		xl[hi] = length(unique(slk[,1]))
		yl[hi] = length(unique(slk[,2]))
		bgam0[slki,4] = hi
		bgam0[slki,5] = length(unique(slk[,1]))
		bgam0[slki,6] = length(unique(slk[,2]))
	}
	bgamo = bgam
	ngamo = ngam
	iout = which((xl==1)|(yl==1))
	ilist = NULL
	iof = NULL
	if(length(iout)>1){ # prune all but edge slices:
		if(!((length(iout)==2)&(diff(iout)[1]==1))){ 
			# ie do the following unless we have a case of single-point edges
			dxl = diff(xl)
			ik = 1
			while((xl[ik]==1)&(dxl[ik]==0)){
				iof = c(iof, which(bgam0[,4]==ik))
				ilist = c(ilist, ik)
				ik = ik+1
			}
			ik = length(xl)
			while((xl[ik]==1)&(dxl[(ik-1)]==0)){
				iof = c(iof, which(bgam0[,4]==ik))
				ilist = c(ilist, ik)
				ik = ik-1
			}
			if(!is.null(ilist)){
				iout = sort(iof)
				bgamo = bgam[-iout,]
				ngamo = ngam[-iout,]
			}
		}
	}
	return(list(bgam=bgamo,ngam=ngamo,ilist=ilist,iof=iof))
}

# Find corresponding boundary radius for that (h,phi) grid point:
f <- function(r,ap,bp,tp){
	x = tp+r/bp
	return((ap*((x^(alpha-1)*exp(-x)))-cutoff))
}
bnd.r <- function(th,h,phi,L=2000){
# Boundary radii given theta and (h,phi):
	x1 = ecbs(h,ab[1],ab[2],Jh)
	x2 = etpb12(c(phi,h),phia,phib,Jphi,ab[1],ab[2],Jh) 
	x3 = epcbs(phi,phia,phib,Jphi)
	ap = sum(th[a.inds]*x2)
	bp = sum(th[b.inds]*x2)
	tp = sum(th[tau.inds]*x1)
	r.tgt = 0
	rmin = 0
	rmax = 200
	rg = seq(rmin,rmax,l=L)
	fg = 0*rg
	for(ir in 1:L){
		fg[ir] = f(rg[ir],ap,bp,tp)
	}
	ir = which(abs(fg)<1e-3)
	if(!length(ir)){
		ir = which.min(abs(fg))
	}
	r.tgt = rev(rg[ir])[1]	
	return(list(r.tgt=r.tgt,ap=ap,bp=bp,tp=tp))
}


render.meshing <- function(cutoff,tit="",...){
	phig = seq(phia,phib,l=Jphi)
	hg = seq(ab[1],ab[2],l=Jh)
	rs1 = matrix(0,nc=Jphi,nr=Jh)
	ap = matrix(NA,nc=Jphi,nr=Jh)
	bp = matrix(NA,nc=Jphi,nr=Jh)
	tp = matrix(NA,nc=Jphi,nr=Jh)
	for(j in 1:Jphi){
		for(k in 1:Jh){
			phi = phig[j]
			h = hg[k]
			bnd.1 = bnd.r(theta1,h,phi)
			rs1[k,j] = bnd.1$r.tgt
			ap[k,j] = bnd.1$ap
			bp[k,j] = bnd.1$bp
			tp[k,j] = bnd.1$tp
		}
	}

	out = eroi$msg.out
	zs  = out$z
	ms  = out$m
	md  = 0*ms
	ths = out$th
	bgam1 = bpts(zs,ms,ths,rs1)
	ngam1 = npts(zs,ms,ths,rs1)
	bn1 = prune.xv(bgam1,ngam1)
	bgam1 = bn1$bgam
	ngam1 = bn1$ngam

	# final Gamma volume
	xv=c(0,1,0);rrx=1.1
	viewxv(xv,cbind(bgam1[,1],bgam1[,2],bgam1[,3]),cbind(ngam1[,1],ngam1[,2],ngam1[,3]),
			"x3'","x1'",tit,rrx,...) 
	xv=c(1,0,0);rrx=1.1
	viewxv(xv,cbind(bgam1[,1],bgam1[,2],bgam1[,3]),cbind(ngam1[,1],ngam1[,2],ngam1[,3]),
			"x3'","x2'"," ",rrx,...)
	xv=c(0,0,1);rrx=1.1
	viewxv(xv,cbind(bgam1[,1],bgam1[,2],bgam1[,3]),cbind(ngam1[,1],ngam1[,2],ngam1[,3]),
			"x2'","x1'"," ",rrx,...)	
	return(list(rs=rs1,ap=ap,bp=bp,tp=tp))
}


sided.skinning <- function(bgam,phi.res=1){
	phig = seq(phia,phib,l=Jphi)
	hg = seq(ab[1],ab[2],l=Jh)
	hgi = hg
	if(length(bn$ilist)){ # remove outer slices if they have been removed elsewhere
		hgi = hg[-bn$ilist] 
	}
	hs = unique(bgam[,3])
	LH = length(hgi)
	xl = yl = numeric(LH)
	dst = dst.acc = matrix(NA,nr=LH,nc=(Jphi-1))
	for(hi in 1:LH){ # slice by slice
		slki = which(bgam[,3]==hs[hi])
		slk = bgam[slki,]
		xlk = length(unique(slk[,1]))
		ylk = length(unique(slk[,2]))
		if((xlk>1)&(ylk>1)){
			for(hj in 1:(Jphi-1)){
				dst[hi,hj] = c(dist(slk[hj:(hj+1),1:2]))
				# sqrt(sum((slk[1,1:2]-slk[2,1:2])^2))
				dst.acc[hi,hj] = sum(dst[hi,1:hj])+c(dist(slk[hj:(hj+1),1:2]))
			}
		}
	}
	# 1st image of unfolded boundary:
	dst.acc = cbind(c(NA,numeric(nrow(dst.acc)-2),NA),dst.acc)
	# ... make up images of relevant intensities
	up = tp+rs1/bp
	zm = ap*(up^(alpha-1)*exp(-up)) # matrix of fitted uptakes
	gm = -(up^(alpha-2)*exp(-up)*(alpha-1-up)) # matrix of growth rates
	zmi = zm
	gmi = gm
	if(length(bn$ilist)){ # remove outer slices if they have been removed elsewhere
		zmi = zm[-bn$ilist,]
		gmi = gm[-bn$ilist,]
	}
	dim(gmi); dim(dst.acc)
	# ... distortion mapping
	LL = phi.res*Jphi
	lg = seq(0,max(dst.acc,na.rm=T),len=LL)
	mp.out = matrix(NA,nr=LH,nc=LL)
	for(hi in 1:LH){ # slice by slice
		if((length(unique(gmi[hi,]))>1)&(length(dst.acc[hi,])>1)){
			mp.out[hi,] = approx(dst.acc[hi,],gmi[hi,],xout=lg)$y
		}
	}
	return(mp.out)
}

skinning <- function(bgam,phi.res=1){
	# first remove outer slices if they have been removed elsewhere
	hgi = hg
	if(length(bn$ilist)){ # remove outer slices if they have been removed elsewhere
		hgi = hg[-bn$ilist] 
	}
	hs = unique(bgam[,3])
	LH = length(hgi)
	LP = length(phig)
	xl = yl = numeric(LH)
	# slice-by-slice analysis of boundary lengthes
	dst = dst.acc = matrix(NA,nr=LH,nc=(Jphi))
	for(hi in 1:LH){ # slice by slice
		# subset
		slki = which(bgam[,3]==hs[hi])
		slk = bgam[slki,]
		xlk = length(unique(slk[,1]))
		ylk = length(unique(slk[,2]))
		# compute distances (chords) and cumulative distances; 
		# here we swipe in both angular directions from some center point (hjz)
		if((xlk>1)&(ylk>1)){
			hjz = floor((Jphi+1)/2)  # angular coordinate of central 'geodesic'
			dstt = as.matrix(dist(slk[,1:2],diag=T,upper=T))[,hjz]
			dst[hi,] = dstt
			dst.acc[hi,] = dstt
			for(hjd in c((hjz+1):Jphi)){ # +pi segment
				dst.acc[hi,hjd] = sum(dstt[hjz:hjd])
			}
			for(hjd in c((hjz-1):1)){ # -pi segment
				dst.acc[hi,hjd] = sum(dstt[hjz:hjd])
			}
		}
	}
	# ... make up images of relevant intensities
	up = tp+rs1/bp
	zm = ap*(up^(alpha-1)*exp(-up)) # matrix of fitted uptakes
	gm = -(up^(alpha-2)*exp(-up)*(alpha-1-up)) # matrix of growth rates
	zmi = zm
	gmi = gm
	if(length(bn$ilist)){ # remove outer slices if they have been removed elsewhere
		zmi = zm[-bn$ilist,]
		gmi = gm[-bn$ilist,]
	}
	# zmi and gmi have values wrt sectors spanning (0,2pi), rather than (-pi,pi); 
	# so we interpolate (periodically) on regular grid centered on central 'geodesic'
	phig2 = seq(-pi,pi,len=LP)
	gmi.c = NA*gmi
	zmi.c = NA*zmi
	ilh = c(1:LH)
	ilh = ilh[!is.na(apply(dst.acc,1,sum))]
	for(i in ilh){
		zmi.c[i,] = spline(x=phig,y=zmi[i,],xout=phig2,method="periodic")$y # fitted uptakes
		gmi.c[i,] = spline(x=phig,y=gmi[i,],xout=phig2,method="periodic")$y # growth gradients
	}
	# ... distortion mapping
	LP = phi.res*Jphi
	# left-right split of sampling
	i.r = c((hjz+1):Jphi)
	i.l = c(1:hjz)
	lg.r = seq(0,max(dst.acc[,i.r],na.rm=T),len=length(i.r))
	lg.l = rev(seq(0,max(dst.acc[,i.l],na.rm=T),len=(length(i.r)+1))[-1])
	lg0 = c(lg.l,lg.r)
	dg = sort(abs(diff(lg0)))[1]
	lg.r = seq(0,max(dst.acc[,i.r],na.rm=T),by=dg)[-1] # (remove 0)
	lg.l = rev(seq(0,max(dst.acc[,i.l],na.rm=T),by=dg))
	lg = c(lg.l,lg.r)
	# left-right split of analysis
	dst.acc.l = dst.acc[,c(1:hjz)]
	dst.acc.r = dst.acc[,c((hjz+1):Jphi)]
	zmi.c.l = zmi.c[,c(1:hjz)]
	zmi.c.r = zmi.c[,c((hjz+1):Jphi)]
	gmi.c.l = gmi.c[,c(1:hjz)]
	gmi.c.r = gmi.c[,c((hjz+1):Jphi)]
	up.out = matrix(NA,nr=LH,nc=length(lg))
	up.out.r = matrix(NA,nr=LH,nc=length(lg.r))
	up.out.l = matrix(NA,nr=LH,nc=length(lg.l))
	mp.out = matrix(NA,nr=LH,nc=length(lg))
	mp.out.r = matrix(NA,nr=LH,nc=length(lg.r))
	mp.out.l = matrix(NA,nr=LH,nc=length(lg.l))
	for(hi in ilh){
		up.out.r[hi,] = approx(x=dst.acc.r[hi,],y=zmi.c.r[hi,],xout=lg.r)$y # +pi segment
		up.out.l[hi,] = approx(x=dst.acc.l[hi,],y=zmi.c.l[hi,],xout=lg.l)$y # -pi segment
		up.out[hi,]	= c(up.out.l[hi,],up.out.r[hi,]) # concatenate
		mp.out.r[hi,] = approx(x=dst.acc.r[hi,],y=gmi.c.r[hi,],xout=lg.r)$y # +pi segment
		mp.out.l[hi,] = approx(x=dst.acc.l[hi,],y=gmi.c.l[hi,],xout=lg.l)$y # -pi segment
		mp.out[hi,]	= c(mp.out.l[hi,],mp.out.r[hi,]) # concatenate
	}
	if(0){ # test for symmetry
		image(t(mp.out),col=gray(c(0:255)/255),axes=T,
			main="growth gradients",xlab="sector (phi)",ylab="transverse elevation (h)")
	}
	# output
	return(list(up.out=up.out,mp.out=mp.out,lg=lg,dst.acc=dst.acc,gmi=gmi,gmi.c=gmi.c))
}
