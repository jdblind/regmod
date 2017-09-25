# Functions for Model and Gradients of Tubular Update Model

library(splines)

#-------------------------- LOCAL FUNS --------------------------
######## Functions for Evaluation of Model and Discrete Laplacian ####

project.voi <- function(xx,ce,s,xi,muk=NULL,doscale=TRUE){
# In: tubemod-funs
# Maps xx |--> xt using (ce,s,xi)
# xx is expected to be of dimension [Nx3]
# muk is the vector of (mu1,mu2) spine coordinates for each voxel
	if(!dosynthetic){
		if(is.null(muk)){muk=matrix(0,nc=2,nr=nrow(xx))}
		roi = cbind(xx[,1]-eroi$c[1],xx[,2]-eroi$c[2],xx[,3]-eroi$c[3])%*%eroi$Gxi
		return(cbind(roi[,1]-muk[,1],roi[,2]-muk[,2],roi[,3]))
	} else {
		project.voi.synthetic(xx,ce,s,xi,doscale)	
	}
}

back.project.voi <- function(xt,ce,s,xi,muk=NULL,doscale=TRUE){
# In: tubemod-funs
# Maps xt |--> xx using (ce,s,xi)
# xt is expected to be of dimension [Nx3]
	if(is.null(muk)){muk=matrix(0,nc=2,nr=nrow(xx))}
	xtm = cbind(xt[,1]+muk[,1],xt[,2]+muk[,2],xt[,3])
	if(doscale){
		xx = as.matrix(xtm,ncol=3) %*% {diag(s) %*% pmat(xi[1],xi[2],xi[3])}
	} else {
		xx = as.matrix(xtm,ncol=3) %*% solve(pmat(xi[1],xi[2],xi[3]))
	}
	return(xx+matrix(ce,nc=3,nr=nrow(xt),byrow=T))
}

project.voi.synthetic <- function(xx,ce,s,xi,doscale=TRUE){
# In: tubemod-funs
# Maps xx |--> xt using (ce,s,xi)
# xx is expected to be of dimension [Nx3]
	if(!dosynthetic){
		return(cbind(xx[,1]-eroi$c[1],xx[,2]-eroi$c[2],xx[,3]-eroi$c[3])%*%eroi$Gxi)		
	} else {
		xxc = xx-matrix(rep(ce,nrow(xx)),byrow=T,ncol=3)
		if(doscale){
			return(xxc %*% crossprod(pmat(xi[1],xi[2],xi[3]),diag(1/s)))		
		} else {
			return(xxc %*% pmat(xi[1],xi[2],xi[3]))
		}
	}
}

back.project.voi.synthetic <- function(xt,ce,s,xi,doscale=TRUE){
# In: tubemod-funs
# Maps xt |--> xx using (ce,s,xi)
# xt is expected to be of dimension [Nx3]
	if(doscale){
		xx = as.matrix(xt,ncol=3) %*% {diag(s) %*% pmat(xi[1],xi[2],xi[3])}
	} else {
		xx = as.matrix(xt,ncol=3) %*% solve(pmat(xi[1],xi[2],xi[3]))
	}
	return(xx+matrix(ce,nc=3,nr=nrow(xt),byrow=T))
}

rph <- function(x,mu1,mu2,ha,hb,Jh,ce,xi,s) {
# In: tubemod-funs
# Given x=(x1,x2,x3) and mapping to tubular parameters,
# get tubular co-ordinates (r,phi,h).
# x is now expected to be in the PAS domain - projection must be done beforehand.
 	v = unlist(c(x))
	h = v[3]
	x1 = ecbs(h,ha,hb,Jh)
 	xp1=v[1]-sum(mu1*x1) 
 	xp2=v[2]-sum(mu2*x1)
	r = sqrt( xp1*xp1+xp2*xp2 ) 
	phi = atan2(xp2,xp1) 
	return(c(r,phi,h))
}

lam <- function(x,a,b,tau,mu1,mu2,ha,hb,Jh,phia,phib,Jphi,ce,xi,s,alpha) {
# In: tubemod-funs
# Evaluate Model 
	x=matrix(x,ncol=3) 
	n=nrow(x) 
	val=NULL
	for(i in 1:n) {
		#print(i)
		v = rph(x[i,],mu1,mu2,ha,hb,Jh,ce,xi,s)   # r,phi,h
		r=v[1]
		phi=v[2]
		h=v[3]
		x1 = ecbs(h,ha,hb,Jh) 
		x2 = etpb12(c(phi,h),phia,phib,Jphi,ha,hb,Jh) 
		u = sum(x1*tau) + r/sum(x2*b)
		amp = sum(x2*a)
		val = c(val,amp * u^(alpha-1) * exp(-u)) 
	}
	return(val)
}

dlam <- function(dx,x,a,b,tau,mu1,mu2,ha,hb,Jh,phia,phib,Jphi,ce,xi,s,alpha,
					v=lam(x,a,b,tau,mu1,mu2,ha,hb,Jh,phia,phib,Jphi,ce,xi,s,alpha)) {
# In: tubemod-funs
# Evaluate Discrete Laplacian of Model 
	x=matrix(x,ncol=3) 
	n=nrow(x)	
	xx=x ; xx[,1] = x[,1]+dx[1]; v1p= lam(xx,a,b,tau,mu1,mu2,ha,hb,Jh,phia,phib,Jphi,ce,xi,s,alpha)
	xx=x ; xx[,1] = x[,1]-dx[1]; v1n= lam(xx,a,b,tau,mu1,mu2,ha,hb,Jh,phia,phib,Jphi,ce,xi,s,alpha)
	xx=x ; xx[,2] = x[,2]+dx[2]; v2p= lam(xx,a,b,tau,mu1,mu2,ha,hb,Jh,phia,phib,Jphi,ce,xi,s,alpha)
	xx=x ; xx[,2] = x[,2]-dx[2]; v2n= lam(xx,a,b,tau,mu1,mu2,ha,hb,Jh,phia,phib,Jphi,ce,xi,s,alpha)
	xx=x ; xx[,3] = x[,3]+dx[3]; v3p= lam(xx,a,b,tau,mu1,mu2,ha,hb,Jh,phia,phib,Jphi,ce,xi,s,alpha)
	xx=x ; xx[,3] = x[,3]-dx[3]; v3n= lam(xx,a,b,tau,mu1,mu2,ha,hb,Jh,phia,phib,Jphi,ce,xi,s,alpha)
	dval = 0 - 6*v + v1p+v1n+v2p+v2n+v3p+v3n
	return(dval)
}

#########  Functions  for Spline Basis Evaluation  #########################

epcbs <- function(xp,a1,b1,J1,verbose=FALSE) {
# In: tubemod-funs
# Evaluate at a vector xp (= phi(x)) a Periodic Cubic B-spline Basis 
#	with J1 elements and equi-spaced knots over the interval [a1,b1]
# xp must elements be in the interval [a1,b1]
	# Must have at least 4 Basis elements
	# xp must be in [a1,b1]
	if(verbose){
		if(min(xp)<a1|max(xp)>b1){
			print("epcbs: xp not in domain...") 
		}
	}
	if(J1<4){
		print("epcbs: J1 must be at least 4") 
	}
	nk=max(4,J1)+1
	knots=(b1-a1)*c(-3,-2,-1,0:(nk-1),nk,nk+1,nk+2)/(nk-1)+a1
	m=splineDesign(knots,xp,outer=T)
	if(length(xp)>1){
		bx = cbind(m[,1]+m[,nk],m[,2]+m[,(nk+1)],m[,3]+m[,(nk+2)],m[,4:(nk-1)])
	} else {
		bx = matrix(c(m[,1]+m[,nk],m[,2]+m[,(nk+1)],m[,3]+m[,(nk+2)],m[,4:(nk-1)]),nrow=1)
	}
	return(bx)
}

ecbs <- function(x,a,b,J,verbose=FALSE) {
# In: tubemod-funs
# Evaluate at a vector x a Cubic B-spline Basis 
#	with J elements and equi-spaced knots over the interval [a,b]
# x is restricted to the nearest point in the interval [a,b]
	# Must have at least 7 Basis elements
	# x restricted to [a,b]
	xx=ifelse(x<a,a,x) 
	xx=ifelse(xx>b,b,xx)
	if(verbose){
		if(min(x)<a|max(x)>b){ 
			print(paste("ecbs: x ",round(x,2)," not in domain (a,) =",a,b,sep="")) 
		}
		if(J<7){ 
			print("ecbs: J must be at least 7") 
		}
	}
	nk=max(J,7)-2
	knots =(b-a)*c(-.1,-.05,-.025,0:(nk-1),nk-1+.025,nk-1+.05,nk-1+.1)/(nk-1)+a
	bx=splineDesign(knots,xx,outer=T)
	return(bx)
}

etpb12 <- function(xph,a1,b1,J1,a2,b2,J2) {
# In: tubemod-funs
# Evaluation at xph (x = c(phi,h)) of a tensor Product Cubic B-spline Basis 
# (periodic in 1 and a-periodic in 2)
# with J1>=4   Periodic B-spline Elements and equi-spaced knots [a1,b1]
#      J2>=7 A-periodic B-spline Elements and equi-spaced knots [a2,b2]
	xph=matrix(c(xph),ncol=2)
	bx1=epcbs(xph[,1],a1,b1,J1)
	bx2=ecbs(xph[,2],a2,b2,J2)
	n=nrow(xph)
	bx=NULL
	for(i in 1:n) {
		bx = rbind(bx,c(bx1[i,] %o% bx2[i,]))
	}
	return(bx)
}

################################# Derivatives ##################################

depcbs <- function(x,a,b,J) {
# In: tubemod-funs
# Evaluate Derivative at x of a Periodic Cubic B-spline Basis 
#	with J elements and equi-spaced knots over the interval [a,b]
# x must elements be in the interval [a,b]
	# Must have at least 4 Basis elements
	# x must be in [a,b]
	if(min(x)<a|max(x)>b){ 
		print("epcbs: x not in domain")
	}
	if(J<4){
		print("epcbs: J must be at least 4") 
	}
	nk = max(4,J)+1
	knots = (b-a)*c(-3,-2,-1,0:(nk-1),nk,nk+1,nk+2)/(nk-1)+a
	nx = length(x)
	m = splineDesign(knots,x,derivs=rep(1,nx),outer=T)
	dbx = cbind(m[,1]+m[,nk],m[,2]+m[,(nk+1)],m[,3]+m[,(nk+2)],m[,4:(nk-1)])	
	return(dbx)
}

decbs <- function(x,a,b,J) {
# In: tubemod-funs
# Evaluate Derivative at x of a Cubic B-spline Basis 
#	with J elements and equi-spaced knots over the interval [a,b]
# x is restricted to the nearest point in the interval [a,b]
	# Must have at least 7 Basis elements
	# x restricted to [a,b]
	if(min(x)<a|max(x)>b){ print("ecbs: x not in domain      ") }
	if(J<7)    		    { print("ecbs: J must be at least 7") }	
	xx=ifelse(x<a,a,x) ; xx=ifelse(xx>b,b,xx) ; nx=length(xx)
	nk=max(J,7)-2; knots =(b-a)*c(-.1,-.05,-.025,0:(nk-1),nk-1+.025,nk-1+.05,nk-1+.1)/(nk-1)+a	
	dbx=splineDesign(knots,xx,derivs=rep(1,nx),outer=T)
	return(dbx)
}

detpb12 <- function(x,a1,b1,J1,a2,b2,J2) {
# In: tubemod-funs
# Evaluate gradient at x of a tensor Product Cubic B-spline Basis (periodic in 1 and a-periodic in 2)
#  with J1>=4   Periodic B-spline Elements and equi-spaced knots [a1,b1]
#       J2>=7 A-periodic B-spline Elements and equi-spaced knots [a2,b2]
	x=matrix(c(x),ncol=2)
	bx1 = epcbs(x[,1],a1,b1,J1) 
	bx2 = ecbs(x[,2],a2,b2,J2)
	dbx1=depcbs(x[,1],a1,b1,J1) 
	dbx2=decbs(x[,2],a2,b2,J2)	
	n=nrow(x) 
	dbx=NULL
	for(i in 1:n) {
		dbx = rbind(dbx,c(c(dbx1[i,]%o% bx2[i,]),c(bx1[i,]%o% dbx2[i,]))) 
	}
	return(dbx)
}

################################ Support Functions

init <- function(xx,z,w,J1,J2) {
#  theta value given xx,z,w and discretization J1 and J2
#  Initialization with Weighted Mean and Covariance Matrix
 	n = nrow(xx)
 	ce = apply(xx*ifelse(z>0,z*w,0),2,sum)/sum(ifelse(z>0,z*w,0))
 	xxc = (xx-matrix(rep(ce,n),ncol=3,byrow=T))
 	xc = xxc*sqrt(ifelse(z>0,z*w,0))
 	SS = as.matrix(t(xc),nrow=3) %*% as.matrix(xc,ncol=3) / sum(ifelse(z>0,z*w,0))
	e = eigen(SS)  
	s = sqrt(e$values)[c(3,2,1)] 
	G = e$vectors[,c(3,2,1)]	
	# Get xi-representation of G 
	xi = euler.angles(G,method=0)
	xi.o = xi # old xi (input)
	# Check : pmat(xi1,xi2,xi3)-G
     
	# Set initial values for Model Parameters
 	xt = project.voi(xx,ce,s,xi)
 	u = hab(xt)
 	a2=u[1] 
 	b2=u[2] 
 	a1=0-pi 
 	b1=pi
 	n1=16
 	n2=16
 	phi=rep(seq(a1,b1,length=n1),n2) 
 	h=rep(seq(a2,b2,length=n2),n1) 
 	h=c(matrix(h,ncol=n2,byrow=T)) 
 	xs=cbind(phi,h) 
 	bx = etpb12(xs,a1,b1,J1,a2,b2,J2)
 	av=1+0*h
 	a = lm(av ~0+ bx)$coef         # Amplitude
 	hm=(a2+b2)/2 
 	hr = (b2-a2)
 	bv=cos(pi*(h-hm)/hr)
 	b = lm(bv ~0+ bx)$coef         # Ellipsoidal Contour
 	tau = rep(max(0,alpha-1),J2)   # Phase
    e=NULL
    e$ce=ce ; e$s=s ; e$xi=xi ; e$xi.o=xi.o ; e$a=a ; e$b=b ; e$tau=tau 
    e$a1=a1 ; e$b1=b1 ;	e$a2=a2 ; e$b2=b2 
    return(e)
}

hab <- function(v,ce,xi,s) {
# In: tubemod-funs
# gets the h-range; when theta is changed this changes
	# v =  S^{-1} P(xi1,xi2,xi3) (x-ce)
 	ha=min(v[,3]) 
 	hb=max(v[,3])
 	return(c(ha,hb)) 	
}

pmat <- function(xi1,xi2,xi3)  {
# In: tubemod-funs
# Create an Orthogonal Matrix from the parameters xi
	
	c1=cos(xi1);s1=sin(xi1);c2=cos(xi2);s2=sin(xi2);c3=cos(xi3);s3=sin(xi3)
	e1= matrix(c( c1*s2,s1*s2,c2),ncol=1);q=qr.Q(qr(e1),complete=T)[,2:3]
 	e2= cos(xi3)*q[,1]+sin(xi3)*q[,2] 
	e3=-sin(xi3)*q[,1]+cos(xi3)*q[,2]
	e2=e2/sign(sum(e2))  ; e3=e3/sign(sum(e3))
		
	P=cbind(e1,e2,e3)
	
	return(P)	
}


eval.rho <- function(xp,mu1,mu2,ha,hb,J2,a1,b1,J1,ce,xi,s,tau,b){
# In: tubemod-funs
# Given (x_i,theta), returns 
#	rho_i = tau_i+r_i/b_i
# i.e. the argument in model shape function g(rho_i).
# Note: tau_i and b_i are voxel-specific as they are evaluated 
# using B-splines for (h_i,phi_i).
#
	v = rph(xp,mu1,mu2,ha,hb,J2,ce,xi,s)   # r,phi,h
	r = v[1]
	phi = v[2]
	h = v[3]
	x1 = ecbs(h,ha,hb,J2)
	x2 = etpb12(c(phi,h),a1,b1,J1,ha,hb,J2) 
	return( sum(x1*tau) + r/sum(x2*b) )
}

eval.amp <- function(xp,mu1,mu2,ha,hb,J2,a1,b1,J1,ce,xi,s,a){
# In: tubemod-funs
# Given (x_i,theta), returns voxel-specific amplitude a_i 
# i.e. the scale in model shape function a_i*g(rho_i).
# Note: a_i is voxel-specific as it is evaluated 
# using B-splines for (h_i,phi_i).
#
	v = rph(xp,mu1,mu2,ha,hb,J2,ce,xi,s)   # r,phi,h
	r = v[1]
	phi = v[2]
	h = v[3]
	x2 = etpb12(c(phi,h),a1,b1,J1,ha,hb,J2) 
	return( sum(x2*a) )
}
