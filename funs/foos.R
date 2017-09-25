rowMatchnum <- function(A,B) { 
# adaptation of rowMatch()...
# Rows in A that match the rows in B 
# The row indexes correspond to A 
    f <- function(...) paste(..., sep=":") 
    if(!is.matrix(B)) B <- matrix(B, 1, length(B)) 
    a <- do.call("f", as.data.frame(A)) 
    b <- do.call("f", as.data.frame(B)) 
    match(b, a) 
} 

rowMatch <- function(A,B) { 
# Author: Ravi Varadhan
# Source: http://r.789695.n4.nabble.com/Matching-a-vector-with-a-matrix-row-td3467355.html
# Rows in A that match the rows in B 
# The row indexes correspond to A 
    f <- function(...) paste(..., sep=":") 
    if(!is.matrix(B)) B <- matrix(B, 1, length(B)) 
    a <- do.call("f", as.data.frame(A)) 
    b <- do.call("f", as.data.frame(B)) 
    match(b, a) 
} 

# ------------------------------------------------------------------------ imaging funs...

ew.blur <- function(img,sblur=0.6){
# In foos.R
# sblur is set "ad hoc"
# require(spatstat)
	return(as.matrix(blur(as.im(img),sigma=sblur)))
}

crosshair <- function(type=1,slcs,dims,col='yellow',lwd=1,lty=1,MUTE=FALSE){
# In foos.R
	if(!MUTE){
		if(type==1){
			abline(v=c(c(1:dims[3])/dims[3])[slcs[3]],col=col,lwd=lwd,lty=lty)
			abline(h=c(c(1:dims[2])/dims[2])[slcs[2]],col=col,lwd=lwd,lty=lty)
		}
		if(type==2){
			abline(v=c(c(1:dims[1])/dims[1])[slcs[1]],col=col,lwd=lwd,lty=lty)
			abline(h=c(c(1:dims[3])/dims[3])[slcs[3]],col=col,lwd=lwd,lty=lty)
		}
		if(type==3){
			abline(v=c(c(1:dims[1])/dims[1])[slcs[1]],col=col,lwd=lwd,lty=lty)
			abline(h=c(c(1:dims[2])/dims[2])[slcs[2]],col=col,lwd=lwd,lty=lty)
		}
	}
}

# -----------------------------------------------------------------------------------
gauss.to.gam <- function(v,alpha,gquant=.98,doplot=FALSE){
# In foos.R	
# Adjusts pseudo-gaussian phase to pseudo-Gamma phase...
# Returns
# 	gv = v*sc.adj+tau.adj (phase for Gamma-profile)
# 	tau.adj
# 	sc.adj
#	fac, such that range(Gaussian.fit/fac) = range(gam.fit)
#
	# 1st pass: initialise Gamma profile features
	xs = seq(0,10,l=10000)
	Gfit = abs(xs)^(alpha-1)*exp(-abs(xs))
	Gmode = abs(xs[which.max(Gfit)])
	qg = qgamma(gquant,shape=alpha)
	# adjust Gaussian phase location 
	tau.adj = Gmode-min(v)-0.5
	# adjust scale of Gaussian phase
	sc.adj = qg/max(v+tau.adj)
	# new phase:
	xg = v*sc.adj+tau.adj
	# Gaussian profile
	gfit = exp(-v^2/2)
	# Gamma profile
	Gfit0 = xg^(alpha-1)*exp(-xg)
	fac = 1.1*max(gfit)/max(Gfit0)
	
	# 2nd pass: refine offset estimates
	Gf <- function(s,d){
		xs = v*s+d
		return(abs(xs)^(alpha-1)*exp(-abs(xs)))
	}
	Nfit = exp(-v^2/2)
	G2G <- function(p){
	# p=(d,s)
		return(sum((Nfit-p[3]*Gf(p[2],p[1]))^2))
	}
	d0 = tau.adj
	s0 = sc.adj
	a0 = fac
	p = c(d0,s0,a0)
	oo = optim(p,G2G)
	p1 = oo$par
	# input and refined initial Gamma fits:
	G0 = Gf(s0,d0)
	G1 = Gf(p1[2],p1[1])
	# Now compare previous and new init Gamma profiles wrt Gauss one:
	if(doplot){
		yl = range(c(Nfit,a0*G0,p1[3]*G1))
		plot(v,Nfit,ylim=yl)
		par(new=TRUE)
		plot(v,a0*G0,col=4,ylim=yl)
		par(new=TRUE)
		plot(v,p1[3]*G1,col=2,ylim=yl)
	}
	xg = c(v*p1[1]+p1[2])
	Gfit = xg^(alpha-1)*exp(-xg)
	return(list(
			# initial objects:
			gv0=xg,gam.fit0=Gfit0,tau.adj0=tau.adj,sc.adj0=sc.adj,fac0=fac,
			# final initial Gamma profile objects:
			gv=xg,gam.fit=Gfit,tau.adj=p1[2],sc.adj=p1[1],fac=p1[3]))
}

# -----------------------------------------------------------------------------------
init.gpars <- function(eroi,g2g){
# In foos.R	
	a0 = eroi$a*g2g$fac
	b0 = 0*eroi$sg
	for(i in 1:dim(eroi$sg)[2]){
		b0[,i] = matrix(eroi$sg[,i]/(eroi$b*g2g$sc.adj),nc=1)
	}
	# then revert b0 [J2xJ1] so as to match x2s[i,] = c([J1xJ2])
	b0 = c(t(b0))
	t0 = g2g$tau.adj-(eroi$b*eroi$tau*g2g$sc.adj)
	return(list(a0=a0,b0=b0,t0=t0))
}

# -----------------------------------------------------------------------------------

mwhich <- function(sample,obj){
# In foos.R	
# Assumes length(sample)>=length(obj)...
	inds = NULL
	for(k in 1:length(obj)){
		inds = c(inds, which(sample==obj[k])[1])
	}
	return(inds[!is.na(inds)])
}

mmatch <- function(sample,obj){
# In foos.R	
	inds = mwhich(sample,obj)
	return(length(inds)==length(obj))
}

# -----------------------------------------------------------------------------------
agrad <- function(rph.in,tau,a,b,alpha){
# In foos.R	
# Exact a-gradient - compares with finite-difference evaluation
	rphs = rph.in$rphs
	x1s=rph.in$x1s
	x2s=rph.in$x2s
	x3s=rph.in$x3s
	n = nrow(rphs)
	vals = numeric(n)
	gvals = matrix(0,nr=nrow(x1s),nc=length(a))
	for(id in 1:n){
		aterm = ifelse(modcon$axs==3,sum(x3s[id,]*a),
				ifelse(modcon$axs==2,sum(x2s[id,]*a),sum(x1s[id,]*a)))
		if(MUTE.G){
			vals[id] = aterm
		} else {
			bterm = ifelse(modcon$bxs==2,sum(x2s[id,]*b),sum(x1s[id,]*b))
			u = sum(x1s[id,]*tau) + rphs[id,1]/bterm
			vals[id] = aterm * u^(alpha-1) * exp(-u)
		}
		if(MUTE.G){
			uval = 1
		} else {
			uval = u^(alpha-1) * exp(-u)
		}
		if(modcon$axs==2){
			gvals[id,] = uval * x2s[id,]
		} else {
			if(modcon$axs==1){
				gvals[id,] = uval * x1s[id,]
			} else {
				gvals[id,] = uval * x3s[id,]
			}
		}
	}
	return(list(vals=vals,gvals=gvals))
}

# -----------------------------------------------------------------------------------
bgrad <- function(rph.in,tau,a,b,alpha){
# In foos.R	
# Exact b-gradient - compares with finite-difference evaluation
#
	rphs = rph.in$rphs
	x1s=rph.in$x1s
	x2s=rph.in$x2s
	x3s=rph.in$x3s
	n = nrow(rphs)
	vals = matrix(0,nr=n,nc=length(b))
	for(id in 1:n){
		aterm = ifelse(modcon$axs==3,sum(x3s[id,]*a),
				ifelse(modcon$axs==2,sum(x2s[id,]*a),sum(x1s[id,]*a)))
		bterm = ifelse(modcon$bxs==2,sum(x2s[id,]*b),sum(x1s[id,]*b))
		u = sum(x1s[id,]*tau) + rphs[id,1]/bterm 
		ax = ((alpha-1)*u^(alpha-2)*exp(-u)-u^(alpha-1)*exp(-u))
		if(modcon$bxs==2){
			vals[id,] = aterm*ax*(-rphs[id,1]/bterm^2)*x2s[id,]
		} else {
			vals[id,] = aterm*ax*(-rphs[id,1]/bterm^2)*x1s[id,]
		}
	}
	return(vals)
}

# -----------------------------------------------------------------------------------
taugrad <- function(rph.in,tau,a,b,alpha){
# In foos.R	
# Exact tau-gradient
#
	rphs = rph.in$rphs
	x1s=rph.in$x1s
	x2s=rph.in$x2s
	x3s=rph.in$x3s
	n = nrow(rphs)
	vals = matrix(0,nr=n,nc=length(tau))
	for(id in 1:n){
		aterm = ifelse(modcon$axs==3,sum(x3s[id,]*a),
				ifelse(modcon$axs==2,sum(x2s[id,]*a),sum(x1s[id,]*a)))
		bterm = ifelse(modcon$bxs==2,sum(x2s[id,]*b),sum(x1s[id,]*b))
		u = sum(x1s[id,]*tau) + rphs[id,1]/bterm 
		ax = ((alpha-1)*u^(alpha-2)*exp(-u)-u^(alpha-1)*exp(-u))
		if(modcon$bxs==2){
			vals[id,] = aterm*ax*x1s[id,]
		} else {
			vals[id,] = aterm*ax*x1s[id,]
		}
	}
	return(vals)
}
