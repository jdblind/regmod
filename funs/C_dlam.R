# Module designed to split computations between R and C code (May 2014).

# ........................................ R functions ........................................

ijk.lkup.tab <- function(xt,gxh,gyh,gzh){
# In: C_dlam
# Generates a lookup table for xt...
# To be passed onto C-coded function eval.dlam.C
	tab = matrix(0,nc=3,nr=nrow(xt))
	for(i in 1:nrow(xt)){
		tab[i,] = lookup.ijk(xt[i,],gxh,gyh,gzh)
	}
	return(tab)
}

xyz.lkup.tab <- function(ijk.tab,xth,gxh,gyh,gzh){
# In: C_dlam
# Generates a lookup table for ijk...
# To be passed onto C-coded function eval.dlam.C
	tab = matrix(0,nc=1,nr=nrow(xt))
	for(i in 1:nrow(ijk.tab)){
		ijk = ijk.tab[i,]
		tab[i] = lookup.xyz(ijk,xth,gxh,gyh,gzh)[1]
	}
	return(tab)
}

eval.dlam.h <- function(xyz.tab,
							xyz.v1p.tab,xyz.v1n.tab,
							xyz.v2p.tab,xyz.v2n.tab,
							xyz.v3p.tab,xyz.v3n.tab,
							LX){
# In: C_dlam
# Evaluates discrete laplacian of lambda at all points in xx;
# xh is the augmented VOI used for computation of finite differences.
# NB: the pmax(0,LX[...]) is used to trick possible numeric(0) values
# which may arise if the point is not found in LX (this can be the 
# case e.g. when voxels inside the volume are missing).
	nx = nrow(xyz.tab)
	dval = numeric(nx)
	z0 = 0
	vrai = TRUE
	for(i in 1:length(xyz.tab)){
		v   = LX[xyz.tab[i]]
		v1p = max(z0,LX[xyz.v1p.tab[i]],na.rm=vrai)
		v1n = max(z0,LX[xyz.v1n.tab[i]],na.rm=vrai)
		v2p = max(z0,LX[xyz.v2p.tab[i]],na.rm=vrai)
		v2n = max(z0,LX[xyz.v2n.tab[i]],na.rm=vrai)
		v3p = max(z0,LX[xyz.v3p.tab[i]],na.rm=vrai)
		v3n = max(z0,LX[xyz.v3n.tab[i]],na.rm=vrai)
		dval[i] = - 6*v + v1p+v1n+v2p+v2n+v3p+v3n
	}
	return(dval)
}

dlam.wrap <- function(xx,xxh,gxh,gyh,gzh,LXH){
# In: C_dlam
	ijk.tab = ijk.lkup.tab(xx,gxh,gyh,gzh)
	# new code:
	xyz.tab = matrix(0,nc=1,nr=nrow(xx))
	xyz.v1p.tab = xyz.v1n.tab = matrix(0,nc=1,nr=nrow(xx))
	xyz.v2p.tab = xyz.v2n.tab = matrix(0,nc=1,nr=nrow(xx))
	xyz.v3p.tab = xyz.v3n.tab = matrix(0,nc=1,nr=nrow(xx))
	for(i in 1:nrow(xx)){
		ijk = ijk.tab[i,]
		xyz.tab[i] = lookup.xyz(ijk,xxh,gxh,gyh,gzh)[1]
		xyz.v1p.tab[i] = lookup.xyz(ijk+matrix(c(1,0,0),nc=3,nr=nrow(xx),byr=T),xxh,gxh,gyh,gzh)[1]
		xyz.v1n.tab[i] = lookup.xyz(ijk-matrix(c(1,0,0),nc=3,nr=nrow(xx),byr=T),xxh,gxh,gyh,gzh)[1]
		xyz.v2p.tab[i] = lookup.xyz(ijk+matrix(c(0,1,0),nc=3,nr=nrow(xx),byr=T),xxh,gxh,gyh,gzh)[1]
		xyz.v2n.tab[i] = lookup.xyz(ijk-matrix(c(0,1,0),nc=3,nr=nrow(xx),byr=T),xxh,gxh,gyh,gzh)[1]
		xyz.v3p.tab[i] = lookup.xyz(ijk+matrix(c(0,0,1),nc=3,nr=nrow(xx),byr=T),xxh,gxh,gyh,gzh)[1]
		xyz.v3n.tab[i] = lookup.xyz(ijk-matrix(c(0,0,1),nc=3,nr=nrow(xx),byr=T),xxh,gxh,gyh,gzh)[1]
	}
	#
	res =  eval.dlam.h(xyz.tab,xyz.v1p.tab,xyz.v1n.tab,
						xyz.v2p.tab,xyz.v2n.tab,
						xyz.v3p.tab,xyz.v3n.tab,
						LXH)					
	return(res)
}

# ---------------------------------------------------------------------------------------------------------

dlam.C <- function(xyz.tab, xyz.v1p.tab,xyz.v1n.tab,
							xyz.v2p.tab,xyz.v2n.tab,
							xyz.v3p.tab,xyz.v3n.tab,
							LX){
# In: C_dlam
# Evaluates discrete laplacian of lambda at all points in xx;
# xh is the augmented VOI used for computation of finite differences.
# NB: the pmax(0,LX[...]) is used to trick possible numeric(0) values
# which may arise if the point is not found in LX (this can be the 
# case e.g. when voxels inside the volume are missing).
	nx = nrow(xyz.tab)
	dval = numeric(nx)
	out <- .C("c_dlam_old",xyztab=as.integer(xyz.tab),
						v1ptab=as.integer(xyz.v1p.tab),
						v1ntab=as.integer(xyz.v1n.tab),
						v2ptab=as.integer(xyz.v2p.tab),
						v2ntab=as.integer(xyz.v2n.tab),
						v3ptab=as.integer(xyz.v3p.tab),
						v3ntab=as.integer(xyz.v3n.tab),
						LX=as.double(LX),
						n=as.integer(nx),
						dval=as.double(dval))
	return(out$dval)
}

eval.phase.C <- function(rphs,x1s,x2s,tau,b){
# In: C_dlam
#
	vals = numeric(nrow(rphs))
	out <- .C("c_phase", rs = as.double(rphs[,1]), n = as.integer(nrow(rphs)),
				x1s = as.double(c(x1s)), n1s = as.integer(ncol(x1s)), 
				x2s = as.double(c(x2s)), n2s = as.integer(ncol(x2s)), 
				tau = as.double(tau), b = as.double(b),
				vals = as.double(vals))
	return(out$vals)
}

eval.g.C <- function(rphs,x1s,x2s,tau,b,alpha){
# In: C_dlam
#
	vals = numeric(nrow(rphs))
	out <- .C("c_g", rs = as.double(rphs[,1]), n = as.integer(nrow(rphs)),
				x1s = as.double(c(x1s)), n1s = as.integer(ncol(x1s)), 
				x2s = as.double(c(x2s)), n2s = as.integer(ncol(x2s)), 
				tau = as.double(tau), b = as.double(b),
				alpha = as.integer(alpha), vals = as.double(vals))
	return(out$vals)
}

eval.galam.C <- function(rphs,x1s,x2s,tau,a,b,alpha){
# In: C_dlam
#
	# vals = c(matrix(0,nr=n,nc=ncol(x2s)))
	vals = numeric(n*ncol(x2s))
	out <- .C("c_galam", rs = as.double(rphs[,1]), n = as.integer(nrow(rphs)),
				x1s = as.double(c(x1s)), n1s = as.integer(ncol(x1s)), 
				x2s = as.double(c(x2s)), n2s = as.integer(ncol(x2s)), 
				tau = as.double(tau), b = as.double(b),
				alpha = as.integer(alpha), vals = as.double(vals))
	return(matrix(out$vals,nr=n))
}
eval.gblam.C <- function(rphs,x1s,x2s,tau,a,b,alpha){
# In: C_dlam
#
	# vals = c(matrix(0,nr=n,nc=ncol(x2s)))
	vals = numeric((n*ncol(x2s)))
	out <- .C("c_gblam", rs = as.double(rphs[,1]), n = as.integer(nrow(rphs)),
				x1s = as.double(c(x1s)), n1s = as.integer(ncol(x1s)), 
				x2s = as.double(c(x2s)), n2s = as.integer(ncol(x2s)), 
				tau = as.double(tau), a = as.double(a), b = as.double(b),
				alpha = as.integer(alpha), vals = as.double(vals))
	return(matrix(out$vals,nr=n))
}
eval.gtlam.C <- function(rphs,x1s,x2s,tau,a,b,alpha){
# In: C_dlam
#
	# vals = c(matrix(0,nr=n,nc=ncol(x1s)))
	vals = numeric((n*ncol(x1s)))
	out <- .C("c_gtlam", rs = as.double(rphs[,1]), n = as.integer(nrow(rphs)),
				x1s = as.double(c(x1s)), n1s = as.integer(ncol(x1s)), 
				x2s = as.double(c(x2s)), n2s = as.integer(ncol(x2s)), 
				tau = as.double(tau), a = as.double(a), b = as.double(b),
				alpha = as.integer(alpha), vals = as.double(vals))
	return(matrix(out$vals,nr=n))
}

eval.lam.C <- function(rph,tau,a,b,alpha){
# In: C_dlam
# New implementation: now incorporates body of dlam.C().
# All preliminary loops are run within C code.
#
	if((modcon$axs%in%c(1,3))|(modcon$bxs==1)){
		eval.lam.R(rph,tau,a,b,alpha)
	} else {
		rphs=rph$rphs
		x1s=rph$x1s
		x2s=rph$x2s
		vals = numeric(nrow(rphs))
		out <- .C("c_lam", rs = as.double(rphs[,1]), n = as.integer(nrow(rphs)),
					x1s = as.double(c(x1s)), n1s = as.integer(ncol(x1s)), 
					x2s = as.double(c(x2s)), n2s = as.integer(ncol(x2s)), 
					tau = as.double(tau), a = as.double(a), b = as.double(b),
					alpha = as.integer(alpha), vals = as.double(vals))
		if(sum(is.na(out$vals))||(max(abs(out$vals))==Inf)){
			warning("In C_dlam::eval.lam.C: *** resorting to eval.lam.R due to NA or infinite u-values! ***")
			return(eval.lam.R(rph,tau,a,b,alpha))
		}
		return(out$vals)
	}
}

eval.dlam.C <- function(nghbr,lx){
# In: C_dlam
# Incorporates body of dlam.C().
# Calls on C code c_dlam_fast.
# void c_dlam_fast(int *nghbr, double *lx, int *n, double *dvals)
#
	n = length(lx)
	dvals = 0*numeric(n)
	out <- .C("c_dlam_fast",
				nghbr=as.integer(c(t(nghbr))),lx=as.double(lx),
				n=as.integer(n),dvals=as.double(dvals))
	return(out$dvals)
}
eval.dlam.C.OBSOLETE <- function(xx,xxh,gxh,gyh,gzh,LXH,retmet=0){
# In: C_dlam
# New implementation: now incorporates body of dlam.C().
# All preliminary loops are run within C code.
#
	nx = nrow(xx)
	gx = sort(unique(xx[,2]))
	gy = sort(unique(xx[,1]))
	gz = sort(unique(xx[,3]))
	dval = 0*numeric(nx)
	ijklup = 0*numeric(3*nx)
	indlup = 0*numeric(7*nx)
	vis = 0*numeric(7*nx)
	xyzlup = 0*numeric(7*nx)
	out <- .C("c_dlam",
				x1=as.double(xx[,1]),x2=as.double(xx[,2]),x3=as.double(xx[,3]),
				n=as.integer(nx),
				x1h=as.double(xxh[,1]),x2h=as.double(xxh[,2]),x3h=as.double(xxh[,3]),
				nh=as.integer(nrow(xxh)),
				gx=as.double(gx),gy=as.double(gy),gz=as.double(gz),
				nx=as.integer(length(gx)),ny=as.integer(length(gy)),nz=as.integer(length(gz)), #14
				gxh=as.double(gxh),gyh=as.double(gyh),gzh=as.double(gzh),
				nxh=as.integer(length(gxh)),nyh=as.integer(length(gyh)),nzh=as.integer(length(gzh)), #20
				LX=as.double(LXH),dval=as.double(dval), #22
				ijklup=as.integer(ijklup),xyzlup=as.double(xyzlup),indlup=as.integer(indlup),vis=as.integer(vis))
	if(retmet>0){
		return(list(dval=out$dval,
				ijklup=matrix(out$ijklup+1,nc=3,byr=T),
				xyzlup=matrix(out$xyzlup,nc=7,byr=T),
				indlup=matrix(out$indlup,nc=7,byr=T),
				vis=matrix(out$vis+1,nc=7,byr=T)))
	} else {
		return(out$dval)
	}
}

gradblock.a <- function(orph,th,nghbr,dx){
	n = nrow(orph$rphs)
	lvals = numeric(n)
	dvals = numeric(n)
	out <- .C("c_grad_block_a",rs=as.double(orph$rphs[,1]),n=as.integer(nrow(orph$rphs)),
					x1s=as.double(c(orph$x1s)),n1s=as.integer(ncol(orph$x1s)), 
					x2s=as.double(c(orph$x2s)),n2s=as.integer(ncol(orph$x2s)), 
					tha=as.double(th[a.inds]),thb=as.double(th[b.inds]),tht=as.double(th[tau.inds]),
					p=as.integer(length(a.inds)),
					mdx=as.double(mean(dx)),alpha=as.integer(alpha),nghbr=as.integer(c(t(nghbr))),
					lvals=as.double(lvals),dvals=as.double(dvals))
	return(list(lvals=out$lvals,dvals=out$dvals))
}

eval.grads.rd.R <- function(xx,nghbr,xinds,exinds,th,Jh,Jphi,alpha,gx=NULL,gy=NULL,gz=NULL) {
# In: C_dlam
# nghbr must be the neighborhood of xxh, rather than of xx
# rd - 'real data'. theta = (a,b,tau)
	inds = get.inds(Jphi,Jh)
	if(is.null(gx)){gx = sort(unique(xx[,2]))}
	if(is.null(gy)){gy = sort(unique(xx[,1]))}
	if(is.null(gz)){gz = sort(unique(xx[,3]))}
	n = nrow(xx)
	dx = c(unique(abs(diff(gx)))[1],unique(abs(diff(gy)))[1],unique(abs(diff(gz)))[1])
	
	GL = DGL = matrix(0,nr=nrow(xx),nc=(length(unlist(all.inds))-length(exinds)))
	mark = 1

	if(!sum(is.element(a.inds,exinds))){
		for(j in 1:length(inds$a)){
			eps=.01*mean(dx)/8 
			bp=bn=th[a.inds]
			bp[j]=bp[j]+eps/2 
			bn[j]=bn[j]-eps/2 
			lx.bp = eval.lam.C(orph,th[inds$tau],bp,th[inds$b],alpha)
			dlx.bp = eval.dlam.C(nghbr,lx.bp)
			lx.bn = eval.lam.C(orph,th[inds$tau],bn,th[inds$b],alpha)
			dlx.bn = eval.dlam.C(nghbr,lx.bn)
			#
	 		GL[,c(mark+j-1)] = (lx.bp[xinds]-lx.bn[xinds])/eps
			DGL[,c(mark+j-1)] = (dlx.bp[xinds]-dlx.bn[xinds])/eps
		}	
		mark = mark+length(inds$a)
	}
	
	if(!sum(is.element(b.inds,exinds))){
		for(j in 1:length(inds$b)){
			eps=.01*mean(dx)/8 
			bp=bn=th[inds$b]
			bp[j]=bp[j]+eps/2 
			bn[j]=bn[j]-eps/2 
			lx.bp = eval.lam.C(orph,th[inds$tau],th[inds$a],bp,alpha)
			dlx.bp = eval.dlam.C(nghbr,lx.bp)
			lx.bn = eval.lam.C(orph,th[inds$tau],th[inds$a],bn,alpha)
			dlx.bn = eval.dlam.C(nghbr,lx.bn)
			#
	 		GL[,(mark+j-1)] = (lx.bp[xinds]-lx.bn[xinds])/eps
			DGL[,(mark+j-1)] = (dlx.bp[xinds]-dlx.bn[xinds])/eps
		}	
		mark = mark+length(inds$b)
	}
	
	if(!sum(is.element(tau.inds,exinds))){
		for(j in 1:length(inds$tau)){
			eps=.01*mean(dx)/8 
			# eps=.001*mean(tau)
			bp=bn=th[tau.inds]
			bp[j]=bp[j]+eps/2 
			bn[j]=bn[j]-eps/2 
			lx.bp = eval.lam.C(orph,bp,th[inds$a],th[inds$b],alpha)
			dlx.bp = eval.dlam.C(nghbr,lx.bp)
			lx.bn = eval.lam.C(orph,bn,th[inds$a],th[inds$b],alpha)
			dlx.bn = eval.dlam.C(nghbr,lx.bn)
			#
	 		GL[,(mark+j-1)] = (lx.bp[xinds]-lx.bn[xinds])/eps
			DGL[,(mark+j-1)] = (dlx.bp[xinds]-dlx.bn[xinds])/eps
		}	
		mark = mark+length(inds$tau)
	}
		
	return(list(DG=GL,DGL=DGL))
}

eval.grads.rd.Cpp <- function(xx,nghbr,xinds,exinds,th,Jh,Jphi,alpha,gx=NULL,gy=NULL,gz=NULL) {
# alias...
	eval.grads.rd(xx,nghbr,xinds,exinds,th,Jh,Jphi,alpha,gx,gy,gz)
}
eval.grads.rd <- function(xx,nghbr,xinds,exinds,th,Jh,Jphi,alpha,gx=NULL,gy=NULL,gz=NULL) {
# nghbr=nghbr.xxh
# exinds=extauinds
# th=theta0
# In: C_dlam [C++ version]
# nghbr must be the neighborhood of xxh, rather than of xx
	
	if(sum(is.element(a.inds,exinds))&&
		sum(is.element(b.inds,exinds))&&
		(!sum(is.element(tau.inds,exinds)))){
		return(eval.grads.rd.R(xx,nghbr,xinds,exinds,th,Jh,Jphi,alpha,gx,gy,gz))
	}
	inds = get.inds(Jphi,Jh)
	if(is.null(gx)){gx = sort(unique(xx[,2]))}
	if(is.null(gy)){gy = sort(unique(xx[,1]))}
	if(is.null(gz)){gz = sort(unique(xx[,3]))}
	n = nrow(xx)
	dx = c(unique(abs(diff(gx)))[1],unique(abs(diff(gy)))[1],unique(abs(diff(gz)))[1])
	
	GL = DGL = matrix(0,nr=nrow(xx),nc=(length(unlist(all.inds))-length(exinds)))
	mark = 1

	if(!sum(is.element(a.inds,exinds))){
		p = length(a.inds)
		lvals = numeric(n*p)
		dvals = numeric(n*p)
		blocka(orph$rphs[,1], c(orph$x1s), c(orph$x2s), modcon$axs, modcon$bxs,
				th[a.inds],th[b.inds],th[tau.inds], 
				mean(dx),alpha,as.integer(c(t(nghbr))),lvals,dvals) 
		lvals = matrix(lvals,nr=n,byr=F)
		dvals = matrix(dvals,nr=n,byr=F)
 		GL[,c(mark:(mark+p-1))] = lvals
 		DGL[,c(mark:(mark+p-1))] = dvals
		mark = mark+p
	}
	if(!sum(is.element(b.inds,exinds))){
		p = length(b.inds)
		lvals = numeric(n*p)
		dvals = numeric(n*p)
		blockb(orph$rphs[,1], c(orph$x1s), c(orph$x2s), modcon$axs, modcon$bxs,
			th[a.inds],th[b.inds],th[tau.inds], 
			mean(dx),alpha,as.integer(c(t(nghbr))),lvals,dvals)
		lvals = matrix(lvals,nr=n,byr=F)
		dvals = matrix(dvals,nr=n,byr=F)
 		GL[,c(mark:(mark+p-1))] = lvals
 		DGL[,c(mark:(mark+p-1))] = dvals
		mark = mark+p
	}
	if(!sum(is.element(tau.inds,exinds))){
		p = length(tau.inds)
		lvals = numeric(n*p)
		dvals = numeric(n*p)
		blockt(orph$rphs[,1], c(orph$x1s), c(orph$x2s), modcon$axs, modcon$bxs,
			th[a.inds],th[b.inds],th[tau.inds], 
			mean(dx),alpha,as.integer(c(t(nghbr))),lvals,dvals)
		lvals = matrix(lvals,nr=n,byr=F)
		dvals = matrix(dvals,nr=n,byr=F)
 		GL[,c(mark:(mark+p-1))] = lvals
 		DGL[,c(mark:(mark+p-1))] = dvals
		mark = mark+p
	}
		
	return(list(DG=GL,DGL=DGL))
}

