csum <- function(x){
	summary(c(x))
}

elam <- function(){
	rphs = orph$rph
	n = nrow(rphs)
	vals = numeric(n)
	for(id in 1:n){
		bterm = ifelse(modcon$bxs==2,sum(x2s[id,]*thT[b.inds]),sum(x1s[id,]*thT[b.inds]))
		u = sum(x1s[id,]*thT[tau.inds]) + rphs[id,1]/bterm
		vals[id] = u^(alpha-1) * exp(-u)
	}
	return(vals)
}

# ------------------------------------------------------------------------ do.it()
do.it <- function(){
	tic = Sys.time()

	lm.method = modcon$lm.method
	maxit = modcon$maxit
	gammet = modcon$gammet
	gam = modcon$gamval
	nocritchange = 0
	RSSs.nog  = RSSs = RSSs.loops = added = ks = ns = ginds = NULL
	theta.all = theta.vec = trhat = diffths = NULL
	gamma.grid = gamma.crit = rss.critg  = NULL
	lap.vec = normk = NULL
	grid.gammas = c(seq(5e-3,0.995,.05),1.25)
	theta.init = th0
	theta1 = theta.init
	
	rph.i = orph
	rphs = orph$rphs[xinds,]
	x1s = orph$x1s[xinds,]
	x2s = orph$x2s[xinds,]
	x3s = orph$x3s[xinds,]

	Omega.tau = spline.omega1(hab(xth),Jh)
	if(modcon$bxs==1){
		Omega.b = Omega.tau
	} else {
		Omega.b = Omega
	}
	if(modcon$axs==1){
		Omega.a = Omega.tau
	} else {
		Omega.a = Omega
	}
	
	n = nrow(xt)
	M0s = M1s = M01s = RSS.all.steps = NULL
	gamma.grid = gamma.crit = gams = NULL
	thetak.a = thetak.b = thetar.a = thetar.b = thetar.t = thetak.t = NULL
	thetas.all = NULL
	thetas.all.b = NULL
	thetas.all.mu = NULL
	thetas.all.tau = NULL
	
	pe.curr = 100
	kk=0
	nochange = FALSE
	
	while((kk<modcon$maxit)&&(abs(pe.curr)>modcon$pe.min)&&(!nochange)){
		kk = kk+1

		theta00 = theta1 # keep this in case we decide to reject the update
		theta0 = theta1

		LXH.0 = eval.lam.C(rph.i,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)
		GLXH.0 = eval.dlam.C(nghbr.xxh,LXH.0)
		if((kk==1)){
			if(gammet>=0){
				eg = eval.grads.rd(xx,nghbr.xxh,xinds,exinds,theta0,Jh,Jphi,alpha)
			} else {
				eg = list(DG=NULL,DGL=NULL)
			}
		}
		M0s = c(M0s, (sum((z-LXH.0[xinds])^2 + gam*GLXH.0^2)))

		gamg = sort(c( gam, gam*(seq(0.9,0.1,l=4)), gam*(seq(1.1,5,l=5))))
		sgam.o1 = select.gam(eg$DG,eg$DGL,z,gammas=gamg,method=gammet,gamval=gam)
		gam = sgam.o1$gamma
		gamma.grid = rbind(gamma.grid, sgam.o1$gammas)
		gamma.crit = rbind(gamma.crit, sgam.o1$crit)
		if(gam){
			ccrit = (sum((z-LXH.0[xinds])^2 + gam*GLXH.0^2))
		} else {
			ccrit = (sum((z-LXH.0[xinds])^2))
		}
		RSS.all.steps = c(RSS.all.steps, ccrit)

		### a-step
		if(!mmatch(exinds,a.inds)){
			theta0 = theta1
			if(!gam){ # OLS!!
				lx = elam()
				if(modcon$axs==1){lxs = data.frame(lx*x1s)}
				if(modcon$axs==2){lxs = data.frame(lx*x2s)}
				if(modcon$axs==3){lxs = data.frame(lx*x3s)}
				ahat <- as.numeric(lm(z~.+0,lxs)$coef)				
				theta1 = theta0
				theta1[a.inds] = ahat
			} else {	
				exainds = sort(as.numeric(unlist(all.inds)[-a.inds]))
				LXH.0 = eval.lam.C(orph,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)
				GLXH.0 = eval.dlam.C(nghbr.xxh,LXH.0)[xinds]
				if(!gam){
					ok=0
					try({
					Xa = agrad(orph,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)$gvals[xinds,];
					X2.0 = Xa; # dummy setting
					ok=1;
					})
				} 
				if((gam!=0)||(sum(is.na((c(Xa)))))){
					eg = eval.grads.rd(xx,nghbr.xxh,xinds,exainds,theta0,Jh,Jphi,alpha)
					Xa = eg$DG
					X2.0 = eg$DGL
				}
				# matrix to ML-pseudo-inverse:
				# Ma  = crossprod(Xa)+gam*crossprod(X2.0)
				Ma  = fxuprod(Xa)
				if(gam){Ma = Ma+gam*fxuprod(X2.0)}
				oas = regul.beta.wrap(Xa,z,Omega.a,bnd=thL[a.inds],gam=gam,
							XD=X2.0,G0=GLXH.0,b0=theta0[a.inds],do.TPS=doTPS)
				thetak.a = rbind(thetak.a, oas$beta.ols)
				thetar.a = rbind(thetar.a, oas$beta.orls)
				theta1[a.inds] = oas$beta.orls
			}
		}
		
		### b-step
		if(!mmatch(exinds,b.inds)){
			for(i in 1:modcon$nbloops.b){
				theta0 = theta1			
				exbinds = sort(as.numeric(unlist(all.inds)[-b.inds]))
				LXH.a = eval.lam.C(orph,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)
				GLXH.a = eval.dlam.C(nghbr.xxh,LXH.a)[xinds]
				if(gam==0){
					ok=0
					try({
					Xb = bgrad(orph,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)[xinds,];
					X2.0 = NULL; # dummy
					ok=1;
					})
				}
				if((gam!=0)||(sum(is.na(c(Xb))))){
					eg = eval.grads.rd(xx,nghbr.xxh,xinds,exbinds,theta0,Jh,Jphi,alpha)
					Xb = eg$DG
					X2.0 = eg$DGL
				}
				# pseudo-values:
				zstar = z-LXH.a[xinds]+c(Xb%*%theta0[b.inds])
				# matrix to ML-pseudo-inverse:
				# Ma  = crossprod(Xb)+gam*crossprod(X2.0) 
				Ma = fxuprod(Xb)
				if(gam){Ma=Ma+gam*fxuprod(X2.0)}
				ob = regul.beta.wrap(Xb,zstar,Omega.b,bnd=thL[b.inds],
					gam=gam,XD=X2.0,G0=GLXH.a,b0=theta0[b.inds],do.TPS=doTPS)
				thetak.b = rbind(thetak.b, ob$beta.ols)
				thetar.b = rbind(thetar.b, ob$beta.orls)
				theta1[b.inds] = ob$beta.orls
				# theta1 = pmin(thU,theta1)
				thetas.all.b = rbind(thetas.all.b,theta1)
			}
			LXH.1 = eval.lam.C(orph,theta1[tau.inds],theta1[a.inds],theta1[b.inds],alpha)
			GLXH.1 = eval.dlam.C(nghbr.xxh,LXH.1)[xinds]
			if(gam){
				ccrit = (sum((z-LXH.1[xinds])^2 + gam*GLXH.1^2))
			} else {
				ccrit = (sum((z-LXH.1[xinds])^2))
			}
			RSS.all.steps = c(RSS.all.steps, ccrit)
		}
		
		### tau-step (finite-diffs, usual Marquardt step)
		if(!mmatch(exinds,tau.inds)){
			extauinds = sort(as.numeric(unlist(all.inds)[-tau.inds]))
			for(i in 1:modcon$nbloops.tau){
				theta0 = theta1
				LXH.0 = eval.lam.C(orph,theta0[tau.inds],theta0[a.inds],theta0[b.inds],alpha)
				GLXH.0 = eval.dlam.C(nghbr.xxh,LXH.0)[xinds]
				eg = eval.grads.rd(xx,nghbr.xxh,xinds,extauinds,theta0,Jh,Jphi,alpha)
				X1.0 = eg$DG
				X2.0 = eg$DGL
				# pseudo-values:
				zstar = z-LXH.0[xinds]+c(X1.0%*%theta0[tau.inds])
				# matrix to ML-pseudo-inverse:
				# Ma  = crossprod(X1.0)+gam*crossprod(X2.0) 
				Ma  = fxuprod(X1.0)
				if(gam){Ma = Ma+gam*fxuprod(X2.0)}
				ob = regul.tau.wrap(X1.0,zstar,Omega.tau,bnd=thL[tau.inds],
					gam=gam,XD=X2.0,G0=GLXH.0,b0=theta0[tau.inds],do.TPS=doTPS)
				thetak.t = ob$beta.orls
				thetak.t = rbind(thetak.t, ob$beta.ols)
				thetar.t = rbind(thetar.t, ob$beta.orls)				
				theta1[tau.inds] = ob$beta.orls
				theta1[tau.inds] = pmin(thU[tau.inds],pmax(thL[tau.inds],theta1[tau.inds]))
				thetas.all.tau = rbind(thetas.all.tau,theta1)
			}
			LXH.1 = eval.lam.C(orph,theta1[tau.inds],theta1[a.inds],theta1[b.inds],alpha)
			GLXH.1 = eval.dlam.C(nghbr.xxh,LXH.1)[xinds]
			if(gam){
				ccrit = (sum((z-LXH.1[xinds])^2 + gam*GLXH.1^2))
			} else {
				ccrit = (sum((z-LXH.1[xinds])^2))
			}
			RSS.all.steps = c(RSS.all.steps, ccrit)
		}
		
		LXH.1 = eval.lam.C(orph,theta1[tau.inds],theta1[a.inds],theta1[b.inds],alpha)
		GLXH.1 = eval.dlam.C(nghbr.xxh,LXH.1)[xinds]
		if(gam){
			ccrit = (sum((z-LXH.1[xinds])^2 + gam*GLXH.1^2))
		} else {
			ccrit = (sum((z-LXH.1[xinds])^2))
		}
		RSS.all.steps = c(RSS.all.steps, ccrit)
		theta1.a = theta1

		if((kk>1)&&(ccrit>(rev(M1s)[1]))){
			nochange=TRUE
			theta1 = theta00
		} else {
			gams = c(gams, gam)
			M1s = c(M1s, ccrit)
			pe.curr = ((M1s[kk]-M0s[kk])/M0s[kk]*100)
			M01s = c(M01s, pe.curr)
			### record current estimate 
			thetas.all = rbind(thetas.all,theta1)
		}
	} # end of kk-for-loop
	((kk<modcon$maxit)&&(abs(pe.curr)>modcon$pe.min)&&(!nochange))
	if(nochange){
		flagout = 2
	} else {
		if(kk==modcon$maxit){
			flagout = 1
		}
	}
	if(abs(pe.curr)<modcon$pe.min){
		flagout = 0
	}

	M1 = ccrit
	zhat = eval.lam.C(orph,theta1[tau.inds],theta1[a.inds],theta1[b.inds],alpha)[xinds]
	tictoc = (Sys.time()-tic)
	out = list(thT=theta.true,thU=thU,thL=thL,th0=th0,th1=theta1,orph=orph,thetas.all=thetas.all,
			thetak.a=thetak.a,thetak.b=thetak.b,
			M0=M0s[1],M1=M1,M01=(M1-M0s[1])/M0s[1]*100,M0s=M0s,M1s=M1s,M01s=M01s,zhat=zhat,
			tictoc=tictoc,flag=flagout)
	return(out)
}
