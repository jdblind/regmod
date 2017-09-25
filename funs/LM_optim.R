snap.to.grid <- function(x1){
# In: LM_optim
# Needed for "real" VOIs with rounding off issues...
	g1 = sort(unique(x1))
	fg1 = seq(min(g1),max(g1),by=diff(g1)[1])
	fx1 = x1*0
	for(i in 1:length(x1)){
		fx1[i] = fg1[which.min(abs(fg1-x1[i]))]
	}
	return(fx1)
}

voi.to.grid <- function(xx){
# In: LM_optim
# Needed for "real" VOIs with rounding off issues...
	xx[,1] = snap.to.grid(xx[,1])
	xx[,2] = snap.to.grid(xx[,2])
	xx[,3] = snap.to.grid(xx[,3])
	return(xx)
}

rebin <- function(v,Jv,vmin=min(v),vmax=max(v)){
# In: LM_optim
# Re-bins v within Jv bins, within (vmin,vmax).
# Returns the corresponding rebinned vector.
# for some reason... does not work as well as voi.to.grid...
	bins = seq(vmin,vmax,length=Jv)
	ov = outer(v,bins,FUN="-")
	return(bins[apply(abs(ov),1,which.min)])
}	

rebin.roi <- function(yy,nx,ny,nz){
# In: LM_optim
# Re-bins VOI yy...
# Returns the corresponding rebinned volume.
	yy[,1] = rebin(yy[,1],ny)
	yy[,2] = rebin(yy[,2],nx)
	yy[,3] = rebin(yy[,3],nz)
	return(yy)
}	

fprint <- function(ch,method=1,app=TRUE,fil=NA){
# In: LM_optim
	if(method==1){
		print(ch)
	} else {
		if(is.na(fil)){
			fil = paste(getwd(),"/out/last_run.txt",sep="")
		}
		print(ch) # echo
		write(ch,file=fil,append=app)
	}
}

s2 <- function(yn,d,tau) {
# In: LM_optim
# Used in tauv().
# Scaled residual sum of squares, where the scale is d+tau (controlled for numerical stability).
	dtau = d+tau 
	# dm = max(d)*.0000001
	dm = max(d)*.0001
	dtau = ifelse(dtau<dm,dm,dtau)
	ss = sum( (yn/dtau)^2 ) 
	return(ss)
}

tauv <- function(yn,d,delt) { 
# In: LM_optim
# Used in steps().
# Computes step tau (used as scale within scaled RSS SRSS(yn,tau)).
# tau is picked so as to minimize 100*|SRSS(yn,tau)-delt|/delt
	# FALSE POSITION
	taua = max(d)/2 
	fa   = s2(yn,d,taua)-delt 
	if(fa<0) { taub=taua ; fb=fa ; while(fb<0) { taub=taub/100; fb=s2(yn,d,taub)-delt }}
	if(fa>0) { taub=taua ; fb=fa ; while(fb>0) { taub=taub*100; fb=s2(yn,d,taub)-delt }}
	# now fa and fb are of opposite signs
	#
	tol=.01 ; kk=0 ; crit=tol+1
	while((kk<30)&(crit>tol)) {   
		kk=kk+1
		tauc = taub 
		fc=fb
		if(abs(fb-fa)>.1e-9) { 
			tauc = taub - fb*(taub-taua)/(fb-fa) 
			fc = s2(yn,d,tauc)-delt 
		}	
		crit = 100*abs(fc)/delt
		if(fa*fc<0) 		{ taub=tauc ; fb=fc }
		if(fa*fc>0) 		{ taua=tauc ; fa=fc }
		if(abs(fa*fb)<=0) 	{ tau=tauc; crit=0 }
	}	
	if(abs(fa)<abs(fb)) { tau=taua }
	if(abs(fa)>abs(fb)) { tau=taub }
	return(tau)
}

steps <- function(ur,d,dels) {
# In: LM_optim
# Used in opt.LM().
# ur are model fit residuals, d = diag(D) in SVD decomposition.
	yn=ur
	md=.0001*max(d)
	y=yn/ifelse(d<md,md,d)
	mx=sum(y^2)
	nk=length(dels)
	tau=dels*0 
	for(k in 1:nk) { 
		delt = dels[k]*mx 
		tau[k]=tauv(yn,d,delt)
	 }
	return(tau)
}

lm.inverse <- function(M,lm.tau,method=1){
# In: LM_optim 
# Marquardt matrix, i.e. inverse of [M+lm.tau*DM], where
# DM = diag(sqrt(diag(M))) or DM = diag(diag(crossprod(M)))
# (DM is passed on as input to save on computation time).
# Test:
if(0){
	M  = matrix(round(runif(9,1,10)),nc=3)
	M %*% lm.inverse(M,lm.tau=.001,method=0)
	M %*% lm.inverse(M,lm.tau=.001,method=1)
	M %*% lm.inverse(M,lm.tau=.001,method=2)
}
#
	OK = FALSE
	if((method!=1)&(method!=2)){
		method=1
		warning("In LM_optim:lm.inverse(): coercing argument 'method' to 1...")
	}
	if(method==1){ 
		# solve (M + tau dM), ie apply tau to all diagonal elements of M
		dM = diag(M)
		DM = diag(dM) # for direct addition of matrices
		Md = M+lm.tau*DM
		svd.out = svd(Md) # Md = U%*%diag(L)%*%t(V)
		L = svd.out$d
		U = svd.out$v %*% tcrossprod(diag(1/L),svd.out$u)
		OK = TRUE
	} 
	if(method==2){ 
		# solve (M + tau I), ie apply maximum tau to all diagonal elements of M
		dM = diag(M)
		DM = diag(rep(1,ncol(M))) # for direct addition of matrices
		Md = M+lm.tau*max(dM)*DM
		svd.out = svd(Md) # Md = U%*%diag(L)%*%t(V)
		L = svd.out$d
		U = svd.out$v %*% tcrossprod(diag(1/L),svd.out$u)
		OK = TRUE
	} 
	if(!OK){
		warning("In LM_optim:lm.inverse(): first inverse method tried failed...")
		dM = diag(M)
		DM = diag(rep(1,ncol(M))) # for direct addition of matrices
		Md = M+lm.tau*max(dM)*DM
		spd = eigen(Md)
		L = spd$values
		if(sum(abs(L)<1e-16)){stop("lm.inverse(): 0-valued eigenvalues, Md is not invertible...")}
		V = spd$vectors 	# i.e. Mm = V %*% diag(L) %*% t(V)
		U = V%*%{diag(1/L)%*%solve(V)} # the inverse
	}
	return(U)
}

#-------------------------- OPTIMISER --------------------------

select.gam.R <- function(X1,X2,z,ftgt=4,gammas=c(seq(0.5,5,l=12)),method=0,gamval=2){
# In: LM_optim
# Selection of regularization parameter:
# - if method==0, via generalized cross-validation;
# - if method>0, gamma is set according to a target number of dof's, 
# 	namely p/ftgt, where p=ncol(X1).
# - if method<0 then gamma is fixed (= gamval).
#
	if(method<0){
		return(list(method=method, target.dof=4, 
				gamma=gamval, crit.val=c(1,1,1),
				gammas=c(1,1,1)*gamval, 
				# crit.redf=c(1,1,1),
				crit.alt=c(1,1,1),
				crit=c(1,1,1)))
	}
	crit = critg = 0*gammas
	p = ncol(X1)
	tgt = p/ftgt
	gg = 1
	X1T = crossprod(X1)
	X2T = crossprod(X2)
	for(gg in 1:length(gammas)){
		P = X1T + gammas[gg]*X2T
		Psvd = svd(P)
		# hat matrix [NxN]:
		# S = X1 %*% {{tcrossprod(Psvd$v,diag(1/(Psvd$d+1e-10)))} %*% crossprod(Psvd$u,t(X1))}
		S = tcrossprod(Psvd$v,diag(1/(Psvd$d+1e-10)))
		S = tcrossprod(tcrossprod(S,Psvd$u),X1)
		S = X1 %*% S
		if(method!=0){ # criterion using regression effective dof's
			crit[gg] = abs(sum(diag(S))-tgt)
			critg[gg] = abs(sum(1-diag(S))-tgt) # or using residual effective dof's
			mini = which.min(crit)
			gam = gammas[mini]
		} else { # [GolubHeathWahba79]
			n = length(z)
			crit[gg] = (sum(({diag(1,n)-S}%*%z)^2)/n) / ((sum(1-diag(S)))/n)^2 
			mini = which.min(crit)
			gam = gammas[mini]
		}
	}
	# plot(gammas,crit); abline(v=gammas[which.min(crit)],col='orange')
	return(list(method=method, target.dof=tgt, gamma=gam, crit.val=crit[mini],
				gammas=gammas, 
				crit.alt=critg,
				crit=crit))
}

select.gam <- function(X1,X2,z,ftgt=4,gammas=c(seq(0.5,5,l=12)),method=1,gamval=2,doRcpp=T){
# In: LM_optim
# Selection of regularization parameter:
# - if method==0, via generalized cross-validation;
# - if method>0, gamma is set according to a target number of dof's, 
# 	namely p/ftgt, where p=ncol(X1).
# - if method<0 then gamma is fixed (= gamval).
#
	if(method<0){
		return(list(method=method, target.dof=4, 
				gamma=gamval, crit.val=c(1,1,1),
				gammas=c(1,1,1)*gamval, 
				# crit.redf=c(1,1,1),
				crit.alt=c(1,1,1),
				crit=c(1,1,1)))
	}
	crit = critg = 0*gammas
	p = ncol(X1)
	tgt = p/ftgt
	gg = 1
	if(doRcpp){
		X1T = fxuprod(X1)
		X2T = fxuprod(X2)
	} else {
		X1T = crossprod(X1)
		X2T = crossprod(X2)
	}
	for(gg in 1:length(gammas)){
		P = X1T + gammas[gg]*X2T
		Psvd = svd(P)
		# hat matrix [NxN]:
		if(!doRcpp){
			S = X1 %*% {{tcrossprod(Psvd$v,diag(1/(Psvd$d+1e-10)))} %*% crossprod(Psvd$u,t(X1))}
			S = tcrossprod(Psvd$v,diag(1/(Psvd$d+1e-10)))
			S = tcrossprod(tcrossprod(S,Psvd$u),X1)
			S = X1 %*% S #matmult.C(X,S) #
		} else {
			S = fprod(X1,ftxprod(ftxprod(ftxprod(Psvd$v,diag(1/(Psvd$d+1e-10))),Psvd$u),X1))
		}
		if(method!=0){ # criterion using regression effective dof's
			crit[gg] = abs(sum(diag(S))-tgt)
			critg[gg] = abs(sum(1-diag(S))-tgt) # or using residual effective dof's
			mini = which.min(crit)
			gam = gammas[mini]
		} else { # [GolubHeathWahba79]
			n = length(z)
			crit[gg] = (sum(({diag(1,n)-S}%*%z)^2)/n) / ((sum(1-diag(S)))/n)^2 
			mini = which.min(crit)
			gam = gammas[mini]
		}
	}
	return(list(method=method, target.dof=tgt, gamma=gam, crit.val=crit[mini],
				gammas=gammas, 
				crit.alt=critg,
				crit=crit))
}

get.inds <- function(J1,J2){
# In: LM_optim
# J1 = Jphi, J2 = Jh
	c.inds = c(1:3)
	s.inds = c(4:6)
	xi.inds = c(7:9)
	mu1.inds = c(10:(J2+9))
	mu2.inds = c(mu1.inds+J2)
	if(modcon$axs==2){
		a.inds = c((mu2.inds[J2]+1):(mu2.inds[J2]+J1*J2))
	} else {
		if(modcon$axs==1){
			a.inds = c((mu2.inds[J2]+1):(mu2.inds[J2]+J2))
		} else {
			a.inds = c((mu2.inds[J2]+1):(mu2.inds[J2]+J1))
		}
	}
	if(modcon$bxs==2){
		b.inds = c((a.inds[length(a.inds)]+1):(a.inds[length(a.inds)]+J2*J1))
	} else {
		b.inds = c((a.inds[length(a.inds)]+1):(a.inds[length(a.inds)]+J2))
	}
	tau.inds = c((b.inds[length(b.inds)]+1):(b.inds[length(b.inds)]+Jh))
	il = list(c=c.inds,s=s.inds,xi=xi.inds,mu1=mu1.inds,
				mu2=mu2.inds,a=a.inds,b=b.inds,tau=tau.inds)
	return(il)
}

ul <- function(v,ss=1.6){
# In: LM_optim
# Used in set.bounds()
	w = smooth.spline(c(1:length(v)),v,df=5)$y
	bp = w+ss*max(abs(v-w))
	bm = w-ss*max(abs(v-w))
	return(list(l=bm,u=bp))
}
#
set.bounds <- function(th,inds,z=1,pc1=0.15,ss=1.6,sab=8){
# In: LM_optim
	dor = modcon$dorep
	bm = bp = th
	# ----- c
	bm[inds$c] = th[inds$c]-max(abs(pc1*th[inds$c]))
	bp[inds$c] = th[inds$c]+max(abs(pc1*th[inds$c]))
	# ----- s
	bm[inds$s] = 1e-8
	if(dor==1){
		bp[inds$s] = 1
	} else {
		bp[inds$s] = th[inds$s]+max(abs(2*th[inds$s]))
	}
	# ----- xi
	bm[inds$xi] = th[inds$xi]-max(abs(pc1*th[inds$xi]))
	bp[inds$xi] = th[inds$xi]+max(abs(pc1*th[inds$xi]))
	# ----- mu's
	if(0){ # spline-hull 
		bb = ul(th[inds$mu1],ss=ss)
		bm[inds$mu1] = bb$l; bp[inds$mu1] = bb$u
		bb = ul(th[inds$mu2],ss=ss)
		bm[inds$mu2] = bb$l; bp[inds$mu2] = bb$u
	} else {
		ii = c(inds$mu1,inds$mu2)
		pcii = .4
		bm[ii] = th[ii]-max(abs(pcii*th[ii]))
		bp[ii] = th[ii]+max(abs(pcii*th[ii]))
	}
	# ----- a: with positivity constraint
	bb = ul(th[inds$a],sab)
	bm[inds$a] = pmax(1e-12,bb$l); bp[inds$a] = bb$u
	# ----- tau: constrain for the Gamma(alpha=3) distribution
	bm[inds$tau] = 0.2; bp[inds$tau] = 8.7
	# ----- b
	bb = ul(th[inds$b],sab)
	bm[inds$b] = pmax(1e-3,bb$l)
	bp[inds$b] = bb$u
	return(list(lower=bm,upper=bp))
}	

# ****************************************************************************************** updated algorithms

find.xinds <- function(xx,xxh){
# In: LM_optim
# Companion function to add.hull.to.roi...
# xx: original VOI
# xh: output of add.hull.to.roi(xx)
# returns indices xinds such that xh[xinds,] == xx
	nh = nrow(xxh)
	nx = nrow(xx)
	if(nx==nh){
		return(c(1:nx))
	}
	zkh = sort(unique(xxh[,3]))
	xinds = numeric(nx)
	for(i in 1:nx){
		xp = xx[i,]
		isl = which(xxh[,3]==xp[3])
		slj = which(xxh[isl,2]==xp[2])
		xinds[i] = isl[slj][which(xxh[isl[slj],1]==xp[1])]
	}
	return(xinds)
}

add.hull.to.roi <- function(xx){
# In: LM_optim
# Sweeps input ROI xx voxel by voxel and adds all adjacent voxel coordinates in
# the 3 generic directions (i+1,i-1,j+1,j-1,k+1,k-1); then unique() this VOI.
# NOTE: we assume that the xx coordinates are "binned".
# Returns augmented ROI.
	###***
	# warning("Now add.hull.to.roi (RD) returns input function!")
	# return(xx)
	###***
	su <- function(v){sort(unique(v))}
	gy <- su(xx[,1])
	gx <- su(xx[,2])
	gz <- su(xx[,3])
	dy <- unique(diff(gy))[1]
	dx <- unique(diff(gx))[1]
	dz <- unique(diff(gz))[1]
	gyh <- c(min(gy)-dy,gy,max(gy)+dy)
	gxh <- c(min(gx)-dx,gx,max(gx)+dx)
	gzh <- c(min(gz)-dz,gz,max(gz)+dz)
	tmp = NULL
	blk = matrix(0,nr=6,nc=3)
	for(i in 1:nrow(xx)){
		ijk = lookup.ijk(xx[i,],gxh,gyh,gzh)
		blk[1,] = c(gyh[ijk[1]-1],gxh[ijk[2]],gzh[ijk[3]])
		blk[2,] = c(gyh[ijk[1]+1],gxh[ijk[2]],gzh[ijk[3]])
		blk[3,] = c(gyh[ijk[1]],gxh[ijk[2]-1],gzh[ijk[3]])
		blk[4,] = c(gyh[ijk[1]],gxh[ijk[2]+1],gzh[ijk[3]])
		blk[5,] = c(gyh[ijk[1]],gxh[ijk[2]],gzh[ijk[3]-1])
		blk[6,] = c(gyh[ijk[1]],gxh[ijk[2]],gzh[ijk[3]+1])
		tmp = rbind(tmp,blk)
	}
	return(unique(rbind(xx,tmp)))
}


eval.x.rph <- function(x,th,ha,hb,Jh,phia,phib,Jphi,alpha,phic=pi,
						lam=eroi$lams,S=eroi$S) {
# In: LM_optim
# Evaluates radial info and splines for all ROI points x.
# x is expected to be in the PAS domain.
# Assumes boundary points have been added beforehand with
#		x = add.hull.to.roi(xt)
# phic: phi-correction to achieve desired phi-range; 
#		atan2 has values in (-pi/2,pi/2] so default correction
#		phic=pi adds pi to achieve range (0,2*pi]
	all.inds = get.inds(Jphi,Jh)
	mu1 = th[all.inds$mu1]
	mu2 = th[all.inds$mu2]
	ce = th[all.inds$c]
	s = th[all.inds$s]
	xi = th[all.inds$xi]
	n=nrow(x)
	val=NULL
	x = unlist(x)
	rphs=0*x
	x1s=matrix(0,nr=n,nc=Jh)
	x2s=matrix(0,nr=n,nc=(Jh*Jphi))
	x3s=matrix(0,nr=n,nc=Jphi)		
	for(i in 1:n) {
		h = x[i,3]
		x1 = ecbs(h,ha,hb,Jh)
		x1s[i,]=x1
	 	xp1 = x[i,1]-sum(mu1*x1)
	 	xp2 = x[i,2]-sum(mu2*x1)
		r = sqrt( xp1*xp1+xp2*xp2 ) 
		phi = atan2(xp2,xp1)+phic
		rphs[i,] = c(r,phi,h)
		x2 = etpb12(c(phi,h),phia,phib,Jphi,ha,hb,Jh) 
		x2s[i,]=x2
		x3 = epcbs(phi,phia,phib,Jphi)
		x3s[i,]=x3
	}
	re.ee = NULL
	return(list(rphs=rphs,x1s=x1s,x2s=x2s,x3s=x3s,re.ee=re.ee))
}

lookup.ijk <- function(xyz,gx,gy,gz){
# In: LM_optim
# Returns c(i,j,k) matching xyz=(y,x,z) based on grids gy,gx,gz.
# xyz is a 3D point.
# NB: a simplified code works here because gx, gy, gz are assumed
# to be sorted and unique.
	return(c(which(gy==xyz[1]),which(gx==xyz[2]),which(gz==xyz[3])))
}

lookup.xyz <- function(ijk,x,gx,gy,gz){
# In: LM_optim
# Finds entry in ROI x matching (i,j,k) position based on grids gy,gx,gz.
# ijk = c(i,j,k).
# Assumes the point exists in x, i.e. that it can be found!
# Example:
# gx = sort(unique(x[,2]))
# gy = sort(unique(x[,1]))
# gz = sort(unique(x[,3]))
# ijk = c(2,5,4)
# c(gy[2],gx[5],gz[4])
# x[lookup.xyz(ijk,x,gx,gy,gz),]
	inds = c(1:nrow(x))
	py = gy[ijk[1]]
	px = gx[ijk[2]]
	pz = gz[ijk[3]]
	if(1){	 # faster
		ik = which(x[,3]==pz)
		xk = x[ik,]
		inds.k = inds[ik]
		ikj = which(xk[,1]==py)	
		xkj = matrix(xk[ikj,],nc=3)
		inds.kj = inds.k[ikj]
		res = inds.kj[which(xkj[,2]==px)]
	} else {
		res = which( (x[,3]==pz)&(x[,2]==px)&(x[,1]==py) )
	}
	if((!length(res))|sum(is.na(res))){
		res = NA
		# stop("error in lookup.xyz: no point found")
	}
	return(res)
}

eval.lam <- function(rph,tau,a,b,alpha){
# In: LM_optim
# Is now a wrapper for corresponding C code
	return(eval.lam.C(rph,tau,a,b,alpha))
}
eval.lam.old <- function(rphs,x1s,x2s,tau,a,b,alpha){
# In: LM_optim
# Is now a wrapper for corresponding C code
	return(eval.lam.C.old(rphs,x1s,x2s,tau,a,b,alpha))
}

eval.lam.R <- function(rph.in,tau,a,b,alpha){
# In: LM_optim
# Includes a safeguard for unusally and improperly large u values...
# ... as a secure backup for eval.dlam.C
	rphs = rph.in$rphs
	x1s=rph.in$x1s
	x2s=rph.in$x2s
	x3s=rph.in$x3s
	n = nrow(rphs)
	vals = numeric(n)
	for(id in 1:n){
		aterm = ifelse(modcon$axs==3,sum(x3s[id,]*a),
					ifelse(modcon$axs==2,sum(x2s[id,]*a),sum(x1s[id,]*a)))
		if(MUTE.G){
			vals[id] = aterm 
		} else {
			bterm = ifelse(modcon$bxs==2,sum(x2s[id,]*b),sum(x1s[id,]*b))
			u = sum(x1s[id,]*tau) + rphs[id,1]/bterm
			u = min(max(u,0),15)
			vals[id] = aterm * u^(alpha-1) * exp(-u)
		}
	}
	return(vals)
}

eval.blam <- function(rphs,x1s,x2s,tau,b,alpha){
# In: LM_optim
# Evaluates lambda profile function at all ROI points, *omitting a*
# (with tubular coordinates rphs[i,]=(r,phi,h)_i) 
# using B-spline components x1s (h-spline) and x2s ((h,phi)-spline)
	n = nrow(rphs)
	vals = numeric(n)
	for(id in 1:n){
		u = sum(x1s[id,]*tau) + rphs[id,1]/sum(x2s[id,]*b)
		vals[id] = u^(alpha-1) * exp(-u)
	}
	return(vals)
}

eval.dlam <- function(xt,xth,gxh,gyh,gzh,LX){
# In: LM_optim
# Evaluates discrete laplacian of lambda at all points in xx;
# xh is the augmented VOI used for computation of finite differences.
# NB: the pmax(0,LX[...]) is used to trick possible numeric(0) values
# which may arise if the point is not found in LX (this can be the 
# case e.g. when voxels inside the volume are missing).
	return(dlam.wrap(xt,xth,gxh,gyh,gzh,LX))
}

eval.dlap <- function(xt,xth,gxh,gyh,gzh,G){
# In: LM_optim
# Evaluates Discrete Laplacian of tabulated function G
# for the case where G is an [NxP] matrix approximating G 
# over a spline basis
	nx = nrow(xt)
	P  = ncol(G)
	dval = matrix(0,nr=nx,nc=P)
	z0 = 0
	for(i in 1:nx){
		ijk = lookup.ijk(xt[i,],gxh,gyh,gzh)
		v   = pmax(z0,matrix(G[lookup.xyz(ijk,xth,gxh,gyh,gzh),],nc=P)[1,],na.rm=T)
		v1p = pmax(z0,matrix(G[lookup.xyz(ijk+c(1,0,0),xth,gxh,gyh,gzh),],nc=P)[1,],na.rm=T)
		v1n = pmax(z0,matrix(G[lookup.xyz(ijk-c(1,0,0),xth,gxh,gyh,gzh),],nc=P)[1,],na.rm=T)
		v2p = pmax(z0,matrix(G[lookup.xyz(ijk+c(0,1,0),xth,gxh,gyh,gzh),],nc=P)[1,],na.rm=T)
		v2n = pmax(z0,matrix(G[lookup.xyz(ijk-c(0,1,0),xth,gxh,gyh,gzh),],nc=P)[1,],na.rm=T)
		v3p = pmax(z0,matrix(G[lookup.xyz(ijk+c(0,0,1),xth,gxh,gyh,gzh),],nc=P)[1,],na.rm=T)
		v3n = pmax(z0,matrix(G[lookup.xyz(ijk-c(0,0,1),xth,gxh,gyh,gzh),],nc=P)[1,],na.rm=T)
		dval[i,] = - 6*v + v1p+v1n+v2p+v2n+v3p+v3n
	}
	return(dval)
}

eval.galam <- function(rphs,x1s,x2s,tau,a,b,alpha,gxh,gyh,gzh){
# In: LM_optim
	eval.galam.C(rphs,x1s,x2s,tau,a,b,alpha)
}
eval.galam.R <- function(rphs,x1s,x2s,tau,a,b,alpha,gxh,gyh,gzh){
# In: LM_optim
# Evaluates gradient of lambda wrt a
	n = nrow(rphs)
	vals = matrix(0,nr=n,nc=ncol(x2s))
	for(id in 1:n){
		u = sum(x1s[id,]*tau) + rphs[id,1]/sum(x2s[id,]*b) 
		vals[id,] = x2s[id,]*u^(alpha-1)*exp(-u)
	}
	return(vals)
}

eval.gblam <- function(rphs,x1s,x2s,tau,a,b,alpha,gxh,gyh,gzh){
# In: LM_optim
	eval.gblam.C(rphs,x1s,x2s,tau,a,b,alpha)
}
eval.gblam.R <- function(rphs,x1s,x2s,tau,a,b,alpha,gxh,gyh,gzh){
# In: LM_optim
# Evaluates gradient of lambda wrt b
	n = nrow(rphs)
	vals = matrix(0,nr=n,nc=ncol(x2s))
	for(id in 1:n){
		u = sum(x1s[id,]*tau) + rphs[id,1]/sum(x2s[id,]*b) 
		amp = sum(x2s[id,]*a)
		ax = ((alpha-1)*u^(alpha-2)*exp(-u)-u^(alpha-1)*exp(-u))
		vals[id,] = amp*ax*(-rphs[id,1]/sum(x2s[id,]*b)^2)*x2s[id,]
	}
	return(vals)
}

eval.gtlam <- function(rphs,x1s,x2s,tau,a,b,alpha,gxh,gyh,gzh){
# In: LM_optim
	eval.gtlam.C(rphs,x1s,x2s,tau,a,b,alpha)
}
eval.gtlam.R <- function(rphs,x1s,x2s,tau,a,b,alpha,gxh,gyh,gzh){
# In: LM_optim
# Evaluates gradient of lambda wrt tau
	n = nrow(rphs)
	vals = matrix(0,nr=n,nc=ncol(x1s))
	for(id in 1:n){
		u = sum(x1s[id,]*tau) + rphs[id,1]/sum(x2s[id,]*b) 
		amp = sum(x2s[id,]*a)
		vals[id,] = amp*((alpha-1)*u^(alpha-2)*exp(-u)-u^(alpha-1)*exp(-u))*x1s[id,]
	}
	return(vals)
}

# ******************************************************************************************
opt.LM <- function(theta0,thetaL,thetaU,z,xx,w,alpha,maxit,dx,a1,b1,J1,ha,hb,J2,gam=0) {
#
# Arguments:
# 	theta0: 	initial value for theta (rescaled beta)
# 	thetaL: 	Lower bound for optimisation wrt theta
# 	thetaU: 	Higher bound for optimisation wrt theta
# 	z: 			uptake information
# 	xx: 		VOI listing, Cartesian coordinates of scaled spherical object (e.g. scaled PCA output)
# 	w: 			weights ***(not used for the moment)***
# 	alpha: 		Gamma distribution parameter (for shape function g(u))
# 	maxit: 		maximum number of iterations
# 	dx:			voxel dimensions (used in the finite differences)
# 	a1,b1,J1: 	phi-spline parameters
# 	ha,hb,J2: 	h-spline parameters
# 	gam:		Regularization parameter (overall model smoothing control)
# Values:
# 	returns estimated transformation theta using nonlinear least squares
#####********* add trace of steps taken, etc.
	# place code here when it works!
	return(list(theta.hat=thetam))
}

list.neighborhood <- function(xx,gx=NULL,gy=NULL,gz=NULL){
	if(is.null(gx)){gx = sort(unique(xx[,2]))}
	if(is.null(gy)){gy = sort(unique(xx[,1]))}
	if(is.null(gz)){gz = sort(unique(xx[,3]))}
	n = nrow(xx)
	nghbr = matrix(0,nr=n,nc=6)
	dx = unique(abs(diff(gx)))[1]
	dy = unique(abs(diff(gy)))[1]
	dz = unique(abs(diff(gz)))[1]
	for(i in 1:n){
		xyz = xx[i,]
		ix = which((xx[,1]==xyz[1])&(xx[,3]==xyz[3]))
		sx = matrix(xx[ix,],nc=3)
		iy = which((xx[,2]==xyz[2])&(xx[,3]==xyz[3]))
		sy = matrix(xx[iy,],nc=3)
		iz = which((xx[,1]==xyz[1])&(xx[,2]==xyz[2]))
		sz = matrix(xx[iz,],nc=3)
		# v1 = ix[which(sx[,1]==xyz[1]+dx)], etc.
		# alternatively, and also deals with ocnditions at borders:
		v1 = ix[which.min(abs(sx[,2]-(xyz[2]+dx)))] 
		v2 = ix[which.min(abs(sx[,2]-(xyz[2]-dx)))] 
		v3 = iy[which.min(abs(sy[,1]-(xyz[1]+dy)))] 
		v4 = iy[which.min(abs(sy[,1]-(xyz[1]-dy)))] 
		v5 = iz[which.min(abs(sz[,3]-(xyz[3]+dz)))] 
		v6 = iz[which.min(abs(sz[,3]-(xyz[3]-dz)))] 
		nghbr[i,] = c(v1,v2,v3,v4,v5,v6)
	}
	return(nghbr)
}

eval.dlam.R <- function(nghbr,lx){
# In: LM_optim
# Evaluates spatial Laplacian of lx at all points in VOI xx.
# xx can be any shape, not just a raster.
	n = length(lx)
	dvals = numeric(n)
	for(i in 1:n){
		dvals[i] = - 6*lx[i] + sum(lx[c(nghbr[i,])])
	}	
	return(dvals)
}

mcrit <- function(th,gam=0){
# In: crit_nls.R
	ab = hab(xt)
	LXHc = eval.lam.C(orph,th[tau.inds],th[a.inds],th[b.inds],alpha)
	return(sum((z-LXHc[xinds])^2))
}
