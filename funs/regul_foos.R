# --------------------------------------------------------

spline.omega <- function(hab,phiab=c(0,2*pi),Lh=100,Lphi=100){
	ha = hab[1]
	hb = hab[2]
	phia = phiab[1]
	phib = phiab[2]
	phig = seq(phia,phib,l=Lphi)
	hg = seq(ha,hb,l=Lh)
	dh=mean(diff(hg)) 
	dphi=mean(diff(phig))
	nh=length(hg)      
	nphi=length(phig)
	xF = matrix(0,nr=(Lphi*Lh),nc=(Jh*Jphi))
	for(j in 1:Lh){
		for(i in 1:Lphi){
			xF[(i+(j-1)*Lphi),] = etpb12(c(phig[i],hg[j]),phia,phib,Jphi,ha,hb,Jh) 
		}
	}
	return(omegaF(dphi,nphi,dh,nh,xF))
}

spline.omega.W <- function(hab,phiab=c(0,2*pi),Lh=100,Lphi=100){
	ha = hab[1]
	hb = hab[2]
	phia = phiab[1]
	phib = phiab[2]
	phig = seq(phia,phib,l=Lphi)
	hg = seq(ha,hb,l=Lh)
	dh=mean(diff(hg)) 
	dphi=mean(diff(phig))
	nh=length(hg)      
	nphi=length(phig)
	xF = matrix(0,nr=(Lphi*Lh),nc=(Jh*Jphi))
	for(j in 1:Lh){
		for(i in 1:Lphi){
			xF[(i+(j-1)*Lphi),] = etpb12(c(phig[i],hg[j]),phia,phib,Jphi,ha,hb,Jh) 
		}
	}
	return(omegaFW(dphi,nphi,dh,nh,xF))
}

omegaF <- function(dphi,nphi,dh,nh,X) {
# Given Basis elements evaluated on a fine grid
# compute LX = Laplacian of X and Omega = (LX)' (LX)
	p=ncol(X) ; 
	LX=X*0
	for(j in 1:p) { 
		LX[,j] = laplacef(dphi,nphi,dh,nh,X[,j]) 
	}
	Omega = t(LX) %*% LX
	return(Omega)
}

omegaFW <- function(dphi,nphi,dh,nh,X) {
# Given Basis elements evaluated on a fine grid
# compute LX = Laplacian of X and Omega = (LX)' (LX)
	p=ncol(X) ; 
	LX=X*0
	for(j in 1:p) { 
		LX[,j] = laplacef(dphi,nphi,dh,nh,X[,j]) 
	}
	Omega = t(LX) %*% LX
	return(list(Omega=Omega,LX=LX))
}

laplacef <- function(dx,nx,dy,ny,f) {
# Given a function f evaluated on a fine grid
# compute Lf = f_xx + f_yy
	f= matrix(f,ncol=ny)
	lf=f 
	dx2=dx*dx 
	dy2=dy*dy
	for(j in 1:ny) {
		for(i in 1:nx) {
			# Use 2'nd-Difference Approximation
			fij = f[i,j]
			if(i==1) { fij1 = f[i+1,j]  } 
			if(i>1) { fij1=f[i-1,j]}
			if(i==nx){ fij2 = f[nx-1,j] } 
			if(i<nx){ fij2=f[i+1,j]}
			if(j==1) { fij3 = f[i,j+1]  } 
			if(j>1) { fij3=f[i,j-1]}
			if(j==ny){ fij4 = f[i,ny-1] } 
			if(j<ny){ fij4=f[i,j+1]}
			lf[i,j] = (fij1-2*fij+fij2)/dx2 + (fij3-2*fij+fij4)/dy2
		}
	}
	return(c(lf))
}

getGM <- function(X,Omega,gam=0,XD=NULL) {
# In: regul_foos.R
# Smoothing Matrix
	if(is.null(XD)){
		o = svd(X) ; u=o$u ; v=o$v ; d=o$d
		xxh = fprod(v,ftxprod(diag(1/sqrt(d)),v))  # [X'X]^{-1/2}
	} else {
		# M = crossprod(X)+gam*crossprod(XD)
		M = fxuprod(X)+gam*fxuprod(XD)
		o = eigen(M)
		v = o$vectors
		d = o$values
		if((sum(d<0))&&(abs(d[d<0])>1e-3)){
			stop("In getGM(): some seriously negative eigenvalues...")
		}
		d = abs(d) + (1e-6)
		xxh = fprod(v,ftxprod(diag(1/sqrt(d)),v))  # M^{-1/2}
		# xxh = v %*% diag(1/sqrt(d)) %*% t(v)  # M^{-1/2}
	}
	M = xxh %*% Omega %*% xxh 
	M = (M+t(M))/2 # enforce symmetry
	o = eigen(M)
	gam = o$vec 
	lam = o$val+1e-8
	gm = gam %*% diag(lam) %*% t(gam)
	# check: gm should match M
	# e<- new.env()
	# e$gam=gam 
	# e$gm=gm 
	# e$lam=lam
	e = list(gam=gam,gm=gm,lam=lam)
	return(e)
}

dfv <- function(lam,tau){
# Find tau so that
# df = sum_i ( 1/(1+tau *lam_i ) 
	return(sum(1/(1+lam*tau)))
}

taudf <- function(dff,lam){
# Finds optimal regularisation parameter tau such that
# 	dff = sum_i ( 1/(1+tau*lam_i ) 
# Arguments:
# 	dff: desired achieved dof
# 	lam: eigenvalues of (X'X)^{-1}%*%Omega
#
	lo=log(min(lam[lam>0])) 
	hi=log(max(lam))
	tau = exp((lo+hi)/2)
	tau1=tau 
	df1=dfv(lam,tau1)
	maxit=50
	iter=0
	while((df1<dff)&(iter<maxit)){
		tau1=tau1/2 
		df1=dfv(lam,tau1)
		iter=iter+1 
	}
	tau0=tau 
	df0=dfv(lam,tau0)
	maxit=30 
	iter=0
	while((df0>dff)&(iter<maxit)){
		tau0=tau0*2
		df0=dfv(lam,tau0)
		iter=iter+1 
	}
	done=0
	if(dff<df0) { tauv=tau0      ; done=1 }  # Can't Do any better
	if(dff>df1) { tauv=tau1*100  ; done=1 }  # Can't Do any better 
	taun = tauv
	if(done==0) { # Look for value between tau0 and tau1
		maxit=100
		err=.01 
		iter=0 
		cgce=err+100
		while( (iter<maxit) & (cgce>err) ){
			taun = tau1 - (df1-dff)*(tau1-tau0)/(df1-df0)
			dft  = dfv(lam,taun)
			if(dft>dff) { df1 = dft ; tau1=taun }
			if(dft<dff) { df0 = dft ; tau0=taun }
			iter=iter+1
			cgce = (dff-dft)
		}
	}
	return(taun)
}

regul.beta.wrap <- function(X,z,Omega,dff=ncol(X)/2,bnd,gam=0,XD=NULL,G0=NULL,b0=NULL,do.TPS=TRUE){
# In: regul_foos.R	
	lam = 0
	gm = getGM(X,Omega,gam,XD)
	if(do.TPS){ # otherwise leave lam set at 0
		lam = taudf(dff,gm$lam)
	}
	oo = regul.beta(X,z,Omega,gm,lam,bnd,gam=gam,XD=XD,G0=G0,b0=b0)
	return(list(beta.ols=oo$beta.ols,beta.orls=oo$beta.orls,lam=lam,cond.nb=oo$cond.nb))
}

regul.beta <- function(X,z,Omega,gm,lam,bnd,gam=0,XD=NULL,G0=NULL,b0=NULL,do.TPS=TRUE,do.orlm=TRUE){
	lam=lam/1000
	if(is.null(XD)){
		# M = crossprod(X)+lam*Omega
		M = fxuprod(X)+lam*Omega
		svdM = svd(M)
		Mi = fprod(svdM$v,ftxprod(diag(1/svdM$d),svdM$u))
		beta.lam = c(Mi %*% fxprod(X,z))
		# (I+lam*gm$gm)1/2:
		igh = (gm$gam) %*% diag(1+lam*gm$lam) %*% t(gm$gam)
		# (X'X)1/2:
		o = svd(X) ; u=o$u ; v=o$v ; d=o$d # X = uDv'
		xxh = v %*% diag(d) %*% t(v)  # [X'X]^{1/2}
	} else {
		condnb = condnbm = 1e7
		maxiter=10
		iter=0
		while(iter<maxiter){
			condnbm = condnb
			iter = iter+1
			M = fxuprod(X)+lam*Omega
			if(gam){M = M+gam*fxuprod(XD)}
			A = M #Mi%*%Omega
			A = (A+t(A))/2
			o = eigen(A)
			od = o$values
			od = od[od>0]
			condnb = max(od)/min(od)
		}
		Mi = o$vec %*% ftxprod(diag(1/o$values),o$vec)
		yD = G0 - XD%*%b0
		beta.lam = c(Mi %*% (fxprod(X,z) - gam*fxprod(XD,yD)))
		# (I+lam*gm$gm)1/2:
		igh = (gm$gam) %*% diag(1+lam*gm$lam) %*% t(gm$gam)
		# (X'X)1/2:
		o = eigen(fxuprod(X)+gam*fxuprod(XD)) # M = UDU'
		od = o$values
		if((sum(od<0))&&(abs(od[od<0])>1e-3)){
			stop("In regul.beta(): some seriously negative eigenvalues...")
		}
		od = abs(od) + (1e-6)
		xxh = o$vectors %*% ftxprod(diag(sqrt(od)),o$vectors)  # M^{1/2}
	}
	D = igh %*% xxh
	D = (D+t(D))/2
	# now to OLS and ORLS
	yD = D%*%beta.lam
	o = lm(yD~D+0)
	if(sum(is.na(o$coefficients))){
		o$coefficients[is.na(o$coefficients)]=1e-3
	}
	bols = as.numeric(coef(o))
	if(do.orlm){
		ui = diag(1,ncol(X))
		ci = bnd
		bet = o$coefficients
		o$coef[!((ui%*%bet-ci>0))]
		check = sum(ui%*%bet-ci>0)/ncol(X)
		icheck = which((ui%*%bet-ci>0))
		oo = orlm(o,ui=diag(1,ncol(X)),ci=bnd,index=c(1:ncol(X)))
		borls = as.numeric(coef(oo))
	} else {
		borls = bols
	}
	oas=list(beta.ols=bols,beta.orls=borls,cond.nb=cond.nb(D),check=check,icheck=icheck)
	return(oas)
}

marquardt.beta <- function(z,theta0,LXH.0,GLXH.0,X1.0,X2.0,gam,lmtau0,lm.method=1,inds=tau.inds){
	rr1 = LXH.0[xinds]
	rr2 = GLXH.0
	# (projected) pseudo-values:
	rv = crossprod(X1.0,(z-rr1))-gam*crossprod(X2.0,rr2)
	# matrix to ML-pseudo-inverse:
	Ma  = crossprod(X1.0)+gam*crossprod(X2.0) 
	tau.s = unique(pmin(lmtau0*.9*10^c(-5:3),1))
	LXHK = matrix(0,nr=nrow(xxh),nc=length(tau.s))
	GLXHK = matrix(0,nr=nrow(xx),nc=length(tau.s))
	crit.tauk = normkt = numeric(length(tau.s))
	thetas = NULL
	kt=1
	for(kt in 1:length(tau.s)){
		# Marquardt matrix, i.e. we inverse [M+tau*DM]
		if(length(Ma)>1){
			U = lm.inverse(Ma,tau.s[kt],lm.method)
		} else {
			U = 1/(Ma+tau.s[kt])
		}
		thupdate = 0*theta0
		thupdate[inds] = c(U%*%rv)
		thetak = theta0+thupdate
		# need to project structure
		xt  = project.voi(xx,thetak[c.inds],thetak[s.inds],thetak[xi.inds])
		xth  = project.voi(xxh,thetak[c.inds],thetak[s.inds],thetak[xi.inds])
		ab = hab(xt)
		rph.k = eval.x.rph(xth,thetak,ab[1],ab[2],Jh,phia,phib,Jphi,alpha)
		LXHK[,kt] = eval.lam(rph.k,thetak[tau.inds],thetak[a.inds],thetak[b.inds],alpha)
		GLXHK[,kt] = eval.dlam(xx,xxh,gxh,gyh,gzh,LXHK[,kt])
		crit.tauk[kt] = sum((z-LXHK[xinds,kt])^2 + gam*GLXHK[,kt]^2)
		normkt[kt] = sqrt(sum(c(U%*%crossprod(X1.0,diag(1,nrow(X1.0))))^2))					
		thetas = cbind(thetas,thetak)
	}
	kt = which.min(crit.tauk)
	tau.k = tau.s[kt]
	theta1 = thetas[,kt]
	return(list(theta=theta1,tau=tau.k,crit=crit.tauk,taus=tau.s))	
}

marquardt.beta.rd <- function(z,theta0,lx,glx,X1.0,X2.0,gam,lmtau0,lm.method=1,inds=tau.inds){
# In: regul_foos.R
	if(0){
		lx=LXH.0
		glx=GLXH.0
		lm.method=1
		inds=tau.inds
	}
	rr1 = lx[xinds]
	rr2 = glx[xinds]
	# (projected) pseudo-values:
	rv = crossprod(X1.0,(z-rr1))-gam*crossprod(X2.0,rr2)
	# matrix to ML-pseudo-inverse:
	Ma  = crossprod(X1.0)+gam*crossprod(X2.0) 
	tau.s = unique(pmin(lmtau0*.9*10^c(-5:3),1))
	LXHK = matrix(0,nr=nrow(xxh),nc=length(tau.s))
	GLXHK = matrix(0,nr=nrow(xxh),nc=length(tau.s))
	crit.tauk = normkt = numeric(length(tau.s))
	thetas = NULL
	kt=1
	for(kt in 1:length(tau.s)){
		# Marquardt matrix, i.e. we inverse [M+tau*DM]
		if(length(Ma)>1){
			U = lm.inverse(Ma,tau.s[kt],lm.method)
		} else {
			U = 1/(Ma+tau.s[kt])
		}
		thupdate = 0*theta0
		thupdate[inds] = c(U%*%rv)
		thetak = theta0+thupdate
		# update structure
		LXHK[,kt] = eval.lam.C(orph,thetak[tau.inds],thetak[a.inds],thetak[b.inds],alpha)
		GLXHK[,kt] = eval.dlam.R(nghbr.xxh,LXHK[,kt])
		crit.tauk[kt] = sum((z-LXHK[xinds,kt])^2 + gam*GLXHK[xinds,kt]^2)
		normkt[kt] = sqrt(sum(c(U%*%crossprod(X1.0,diag(1,nrow(X1.0))))^2))					
		thetas = cbind(thetas,thetak)
	}
	kt = which.min(crit.tauk)
	tau.k = tau.s[kt]
	theta1 = thetas[,kt]
	return(list(theta=theta1,tau=tau.k,crit=crit.tauk,taus=tau.s))	
}

# ------------------------------------------------------------------------ tau

spline.omega1 <- function(hab,Jh,Lh=100){
	ha = hab[1]
	hb = hab[2]
	hg = seq(ha,hb,l=Lh)
	dh=mean(diff(hg)) 
	xF = matrix(0,nr=Lh,nc=Jh)
	for(i in 1:Lh){
		xF[i,] = ecbs(hg[i],ha,hb,Jh) 
	}
	return(omegaF1(dh,Lh,xF))
}

omegaF1 <- function(dh,nh,X){
# Given Basis elements evaluated on a fine grid
# compute LX = Laplacian of X and Omega = (LX)' (LX)
	p=ncol(X) ; 
	LX=X*0
	for(j in 1:p) { 
		LX[,j] = laplacef1(dh,nh,X[,j]) 
	}
	Omega = t(LX) %*% LX
	return(Omega)
}

laplacef1 <- function(dx,nx,f) {
# Given a 1D function f evaluated on a fine grid
# compute Lf = f_xx 
	L=numeric(nx)
	dx2=dx*dx 
	for(i in 1:nx) {
		# Use 2'nd-Difference Approximation
		fi = f[i]
		if(i==1) { fi1 = f[i+1]  } 
		if(i>1) { fi1=f[i-1]}
		if(i==nx){ fi2 = f[nx-1] } 
		if(i<nx){ fi2=f[i+1]}
		L[i] = (fi1-2*fi+fi2)/dx2 
	}
	return(c(L))
}


regul.tau.wrap <- function(X,z,Omega.tau,dff=ncol(X)/2,bnd,gam=0,XD=NULL,G0=NULL,b0=NULL,do.TPS=TRUE){
if(0){	
X=X1.0
dff=ncol(X)/2
bnd=thL[tau.inds]
gam=gam
XD=X2.0
G0=GLXH.0
b0=theta0[tau.inds]
}
	gm = getGM(X,Omega.tau,gam,XD)
	lam = 0
	if(do.TPS){ # otherwise leave lam set at 0
		lam = taudf(dff,gm$lam)
	}
	oo = regul.beta(X,z,Omega.tau,gm,lam,bnd,gam=gam,XD=XD,G0=G0,b0=b0)
	return(list(beta.ols=oo$beta.ols,beta.orls=oo$beta.orls,lam=lam,cond.nb=oo$cond.nb))
}

# ------------------------------------------------------------------------ condition...

cond.nb <- function(X){
	Ssvd = svd((X))
	return(max(Ssvd$d)/min(Ssvd$d))
}
cond.nb.j <- function(X){
	cnj = numeric(ncol(X))
	for(j in 1:ncol(X)){
		Xj = X[,-j]
		cnj[j] = cond.nb(Xj)
	}
	j = which.min(cnj)
	return(list(j=j,cn=cnj[j],Xj=X[,-j[1]]))
}
reduced.X <- function(X){
	j0 = ncol(X)
	vj = cns = inds = NULL
	Xj = X
	cn = cond.nb(X)
	while((ncol(Xj)>1)&(cn>30)){
		out = cond.nb.j(Xj)
		inds = c(inds,out$j[1])		 # store index 
		vj = cbind(vj,Xj[,out$j[1]]) # store removed vector
		cns = c(cns, out$cn)
		cn = out$cn
		Xj = out$Xj
	}
	return(list(cn=cn,cns=cns,inds=inds,Xj=Xj,j0=j0,j1=ncol(Xj)))		
}
