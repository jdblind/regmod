#! /usr/bin/Rscript 
while( dev.next()>1 ){ dev.off() }
rm(list=ls(pos=.GlobalEnv), pos=.GlobalEnv)

library(splines)
library(gplots)
library(ic.infer)
library(Matrix)
library(plot3D)
library(Rcpp)
library(inline)
library(RcppArmadillo)
dyn.load("./Cpack/libdlam.so")
dyn.load("./Cpack/libmaths.so")
options(CBoundsCheck = TRUE)

source("utils/impro.R")
source("utils/amide.R")
source("utils/tictoc.r")
source("funs/tubemod-funs.R")
source("funs/LM_optim.R")
source("funs/C_dlam.R")
source("funs/C_maths.R")
source("funs/exsarcroi.R")
source("funs/foos.R")
source("funs/plot_phases.R")
source("funs/regul_foos.R")
source("funs/Rcpp_mat.R")
sourceCpp("Cpack/grad_Rcpp.cpp")
source("funs/image_toolbox.R")
source("./optim_routine_exsarcroi.R")

# ------------------------------------------------------------------------ inits...

print(Sys.time())
dosynthetic=TRUE
docoldcore=TRUE
donoise=TRUE
dotrueinit=FALSE
doTPS=FALSE
doMC=TRUE
dopeel=FALSE
dohull=FALSE
MUTE.G = FALSE
n1  = 12
# g-shape parameter (Gamma distribution)
alpha = 3
gam0 = 0.3
modcon = list(astep=TRUE,axs=2,bxs=2,dorep=0,
				gammet = -1, # fixed gamma, no selection
				gamval = gam0, lm.method = 2,
				pe.min = .01, # minimum percent error to be achieved
				maxit = 25, nbloops.b = 1, nbloops.tau = 1)
M = ifelse(doMC,100,1) # number of Monte Carlo runs

# -------------------------------------------------------------------------
# ------------------ create 3D structure and uptake data ------------------
# -------------------------------------------------------------------------

# ------------------
# Initial values
set.seed(1)
# Npoints
nx1 = 16
nx2 = 16
nx3 = 12
(n = nx1*nx2*nx3)
# Nsplines
n3 = n1
n2 = nx3
(n1*n2)
# Values and Gradients of Tensor-Product Basis
Jphi = n1
Jh = n2
phia = 0
phib = 2*pi

# ------------------
# Spatial structure
xxo = NULL
x1 = rep(seq(-12,12,length=nx1),nx2) 
x2 = rep(seq(-12,12,length=nx2),nx1)
x2 = c(matrix(x2,ncol=nx2,byrow=T))
for(k in 1:nx3) { 
	xxo = rbind(xxo,cbind(x1,x2,rep(k-nx3/2,nx1*nx2))) 
}
xx = xxo
xxh = xxo
xinds = c(1:nrow(xx))

# Structural Parameters
xi = 1+0*c( .2,-.3,.9)
ce = 0*c(0.1,-0.2,0.15)
s = 1+0*c(1.1,.9,1.05) 
# xs = cbind(xx[,1]-ce[1],xx[,2]-ce[2],xx[,3]-ce[3]) 
gam = matrix(c(1,0,0,0,1,0,0,0,1),nc=3)

# projected structure
xto = xxo 
xt  = xx 
xth = xxh
ab = hab(xt)
abh = hab(xth)
ha = ab[1]
hb = ab[2]
# Core- Center var1 and var2
mu1=mu2=numeric(Jh)
# -----------------------------------------------------------------------------------
print(paste("N_xx =",nrow(xx)," ; N_xxh =",nrow(xxh)))
gx = sort(unique(xx[,2]))
gy = sort(unique(xx[,1]))
gz = sort(unique(xx[,3]))
gxh = sort(unique(xxh[,2]))
gyh = sort(unique(xxh[,1]))
gzh = sort(unique(xxh[,3]))
nghbr.xx = list.neighborhood(xx,gx,gy,gz)
nghbr.xxh = nghbr.xx	
Omega = spline.omega(ab)

# optimal parameters
if(modcon$axs==1){
	avv = (1+1/(1+(c(1:n2)-(n2/2))^2/10))/2;
	a = 2.1*avv
}
if(modcon$axs==2){
	a = rep(c(seq(15,20,len=Jh/2),seq(20,15,len=Jh/2)),each=Jphi)
}
if(modcon$axs==3){
	avv = (1+1/(1+(c(1:n1)-(n1/2))^2/10))/2;
	a = 2.1*avv
}
##
if(modcon$bxs==1){
	b = c(seq(5.6,4.8,len=Jh/2),seq(4.8,5.6,len=Jh/2))
} else {
	# smooth volume:
	b = rep(c(seq(5.6,4.8,len=Jh/2),seq(4.8,5.6,len=Jh/2)),each=Jphi)
}
##
if(docoldcore){ # cold core
	tau = c(seq(1.2,0.85,len=Jh/2),seq(.8,1.2,len=Jh/2))
} else { # hot core
	tau = c(seq(2,1.8,len=Jh/2),rev(seq(2,1.8,len=Jh/2)))
}

# ------------------
# evaluate splines

orph = eval.x.rph(xth,c(ce,s,xi,mu1,mu2),ha,hb,Jh,phia,phib,Jphi,alpha)
x1s = orph$x1s[xinds,]
x2s = orph$x2s[xinds,]
x3s = orph$x3s[xinds,]

# optimal parameter
theta.opt = c(ce,s,xi,mu1,mu2,a,b,tau)
theta.true = c(ce,s,xi,mu1,mu2,a,b,tau)
thT = theta.true

# indices...
ha = ab[1]-.1; hb = ab[2]+.1
rep.inds = get.inds(Jphi,Jh)
all.inds = get.inds(Jphi,Jh)
a.inds = all.inds$a
tau.inds = all.inds$tau
b.inds = all.inds$b
c.inds = all.inds$c
s.inds = all.inds$s
xi.inds = all.inds$xi
mu1.inds = all.inds$mu1
mu2.inds = all.inds$mu2

# ************************************************************************************************

ii = sort(c(c.inds,s.inds,xi.inds,mu1.inds,mu2.inds))
exinds = c(1:length(unlist(all.inds)))[ii]
pis = c(1:length(unlist(all.inds)))[-ii]
fixstruct = prod(is.element(c(c.inds,s.inds,xi.inds,mu1.inds,mu2.inds),exinds))

theta.true = thT
# Estimation bounds on theta
thUL = set.bounds(thT,all.inds)#,method=0,dor=dorep)
thL = thUL$lower
thU = thUL$upper
# reset bounds for non-optimised params
thL[exinds] = thT[exinds]-1e-12
thU[exinds] = thT[exinds]+1e-12
thL[-exinds] = pmin(thL[-exinds],c(thT[-exinds]-1.5))
thU[-exinds] = pmax(thU[-exinds],c(thT[-exinds]+1.5))
# mu's
mu.inds = sort(c(mu1.inds,mu2.inds))
th0 = thT
# a
if(!mmatch(exinds,a.inds)){
	thL[a.inds] = 1e-3
	thU[a.inds] = max(theta.true[a.inds])*2
	th0[a.inds] = thT[a.inds]/4
}
# b
if(!mmatch(exinds,b.inds)){
	thU[b.inds] = theta.true[b.inds]*2
	thL[b.inds] = pmax(1e-3,theta.true[b.inds]/10)
	th0[b.inds] = min(thT[b.inds])+(max(thT[b.inds])-thT[b.inds])-2
	# th0[b.inds] = thT[b.inds]/5
}
# tau
if(!mmatch(exinds,tau.inds)){
	thL[tau.inds] = rep(.2,length(tau.inds))
	thU[tau.inds] = rep(5,length(tau.inds))
	th0[tau.inds] = thT[tau.inds]/4
}
if(dotrueinit){
	th0 = pmax(thL,pmin(thT,thU))
} else {
	th0 = pmax(thL,pmin(th0,thU))
}
sum(thL>thT)
sum(thT>thU)

# ************************************************************************************************

orph.t = orph
LX = eval.lam.C(orph.t,
				theta.true[tau.inds],theta.true[a.inds],theta.true[b.inds],alpha)
zexact = LX
w = zexact*0 + 1
z = zexact

LX.R = eval.lam.R(orph,th0[tau.inds],th0[a.inds],th0[b.inds],alpha)
LX.C = eval.lam.C(orph,th0[tau.inds],th0[a.inds],th0[b.inds],alpha)
summary(LX.R-LX.C)
LX = LX.C
r1 = orph$rphs[xinds,1]

rphs = orph$rph
n = nrow(rphs)
lx = numeric(n)
for(id in 1:n){
	bterm = ifelse(modcon$bxs==2,sum(x2s[id,]*thT[b.inds]),sum(x1s[id,]*thT[b.inds]))
	u = sum(x1s[id,]*thT[tau.inds]) + rphs[id,1]/bterm
	lx[id] = u^(alpha-1) * exp(-u)
}

zT = zhat0 = numeric(n)
for(id in 1:n){
	aterm = ifelse(modcon$axs==3,sum(x3s[id,]*thT[a.inds]),
			ifelse(modcon$axs==2,sum(x2s[id,]*thT[a.inds]),sum(x1s[id,]*thT[a.inds])))
	zT[id] = aterm*lx[id]
	aterm = ifelse(modcon$axs==3,sum(x3s[id,]*th0[a.inds]),
			ifelse(modcon$axs==2,sum(x2s[id,]*th0[a.inds]),sum(x1s[id,]*th0[a.inds])))
	zhat0[id] = aterm*lx[id]
}
summary(zT-zhat0)
th0[a.inds]-thT[a.inds]
th0[b.inds]-thT[b.inds]
th0[tau.inds]-thT[tau.inds]

NPA = length(unlist(all.inds))-length(exinds) # nb of params being optimised
(fixstruct = as.logical(prod(is.element(c(c.inds,s.inds,xi.inds,mu1.inds,mu2.inds),exinds))))

# ************************************************************************************************
mnoise = 0
snoise = c(1,.75,.5,.25)
if(!doMC){
	snoise = 0
}
K = length(snoise)
noisek = 1

voi.noise <- function(z,m=0,s=.5){
	rnorm(length(z),m=m,s=s)#*sqrt(z)
}
voi.noise.bg <- function(zexact,z){
	abs(rnorm(sum(z<=0)))* zexact[z<=0]
}

thetas.mc = array(NA,dim=c(M,length(thT),K))
res.m = res.v = matrix(NA,nr=M,nc=K)
rss = rss0 = flag.mc = tictoc.mc = matrix(NA,nr=M,nc=K)
sp.a.m = sp.b.m = sp.t.m = matrix(NA,nr=M,nc=K)
sp.a.v = sp.b.v = sp.t.v = matrix(NA,nr=M,nc=K)
sp.a.mse = sp.b.mse = sp.t.mse = matrix(NA,nr=M,nc=K)

aspT = bspT = tspT = numeric(n)
for(id in 1:n){
	aspT[id] = ifelse(modcon$axs==2,sum(x2s[id,]*thT[a.inds]),sum(x1s[id,]*thT[a.inds]))
	bspT[id] = ifelse(modcon$bxs==2,sum(x2s[id,]*thT[b.inds]),sum(x1s[id,]*thT[b.inds]))
	tspT[id] = sum(x1s[id,]*thT[tau.inds])
}

noisek = mc = 1
for(noisek in 1:K){	
	for(mc in 1:M){
		if(donoise){ # add noise
			z = zexact + voi.noise(zexact,m=mnoise,s=snoise[noisek])
			z[z<=0] = zexact[z<=0] + voi.noise.bg(zexact,z)
		}
		summary(z)
		z.long = z
		z = z[xinds]
		
		if(!doMC){
			quartz()
			chut=hist(z,plot=FALSE)
			hist(zexact[xinds],col='navy',breaks=chut$breaks,main="data",xlab="uptake")
			hist(z,col=rgb(.5,.5,.5,alpha=.5),breaks=chut$breaks,add=T)
			legend("topright",col=c('navy','gray'),legend=c("true","noisy"),pch=15)
			cg = gray(c(0:255)/255)
			roi = cbind(z,rep(1,nrow(xx)),xx)
			voi = rasterize.voi(roi,def=NA)
			rr = range(roi[,1])
			quartz()
			par(mfrow=c(3,4))
			for(k in 1:dim(voi)[3]){
				mimage(t(voi[,,k]), rr, main=paste(k))
			}
			quartz()
			par(mfrow=c(3,1))
			plot(a,pch=20,t='b')
			plot(b,pch=20,t='b')
			plot(tau,pch=20,t='b')
		}
		
		# ----------------------------------------------------------------------------
		# run main optim code
		out = do.it()
		save.image(paste("out_simdata/synthetic_run_gam_",gam0,".Rdata",sep=""))
	
		thI = out$th0
		thF = out$th1
		# csum(out$zhat-z)
		# csum(thT[a.inds]-thI[a.inds])
		# csum(thT[a.inds]-thF[a.inds])
		# compute spline MSE's
		aspF = bspF = tspF = numeric(n)
		for(id in 1:n){
			aspF[id] = ifelse(modcon$axs==2,sum(x2s[id,]*thF[a.inds]),sum(x1s[id,]*thF[a.inds]))
			bspF[id] = ifelse(modcon$bxs==2,sum(x2s[id,]*thF[b.inds]),sum(x1s[id,]*thF[b.inds]))
			tspF[id] = sum(x1s[id,]*thF[tau.inds])
		}
		sp.a.mse[mc,noisek] = mean((aspF-aspT)^2)
		sp.b.mse[mc,noisek] = mean((bspF-bspT)^2)
		sp.t.mse[mc,noisek] = mean((tspF-tspT)^2)
		sp.a.m[mc,noisek] = mean((aspF-aspT))
		sp.b.m[mc,noisek] = mean((bspF-bspT))
		sp.t.m[mc,noisek] = mean((tspF-tspT))
		sp.a.v[mc,noisek] = var((aspF-aspT))
		sp.b.v[mc,noisek] = var((bspF-bspT))
		sp.t.v[mc,noisek] = var((tspF-tspT))
				
		thetas.mc[mc,,noisek] = thF
		rss0[mc,noisek] = out$M0
		rss[mc,noisek] = out$M1
		res = z-out$zhat
		res.m[mc,noisek] = mean(res)
		res.v[mc,noisek] = var(res)
		tictoc.mc[mc,noisek] = out$tictoc[[1]]
		flag.mc[mc,noisek] = out$flag
		print(data.frame(noise.k=noisek, 
					mc=mc,
					RSS0=round(out$M0,0), RSS1=round(out$M1,0), # RSSs
					pRSS=round((out$M0-out$M1)/out$M0,3)*100, # %-drop in RSS
					res.m=round(res.m[mc,noisek],6), # stats on residuals
					res.v=round(res.v[mc,noisek],3),
					dur=round(out$tictoc[[1]],2), flag=round(out$flag,0)))
	}
}

if(doMC){ # noiseless case does not require Monte Carlo
	z = zexact[xinds]
	out0 = do.it()
	thF0 = out0$th1
	# compute spline stats
	aspF0 = bspF0 = tspF0 = numeric(n)
	for(id in 1:n){
		aspF0[id] = ifelse(modcon$axs==2,sum(x2s[id,]*thF0[a.inds]),sum(x1s[id,]*thF0[a.inds]))
		bspF0[id] = ifelse(modcon$bxs==2,sum(x2s[id,]*thF0[b.inds]),sum(x1s[id,]*thF0[b.inds]))
		tspF0[id] = sum(x1s[id,]*thF0[tau.inds])
	}
	sp.a.mse = cbind(sp.a.mse, rep(mean((aspF0-aspT)^2),M))
	sp.b.mse = cbind(sp.b.mse, rep(mean((bspF0-bspT)^2)))
	sp.t.mse = cbind(sp.t.mse, rep(mean((tspF0-tspT)^2)))
	sp.a.m = cbind(sp.a.m, rep(mean((aspF0-aspT)),M))
	sp.b.m = cbind(sp.b.m, rep(mean((bspF0-bspT)),M))
	sp.t.m = cbind(sp.t.m, rep(mean((tspF0-tspT)),M))
	sp.a.v = cbind(sp.a.v, rep(var((aspF0-aspT)),M))
	sp.b.v = cbind(sp.b.v, rep(var((bspF0-bspT)),M))
	sp.t.v = cbind(sp.t.v, rep(var((tspF0-tspT)),M))
	# aggregate with previous stats
	snoise = c(snoise,0.)
	res0 = z-out0$zhat
	res.m = cbind(res.m, rep(mean(res0),M))
	res.v = cbind(res.v, rep(var(res0),M))
	rss0 = cbind(rss0, rep(out0$M0,M))
	rss = cbind(rss, rep(out0$M1,M))
}

# ************************************************************************************************
# output analysis

if(doMC){
	quartz()
	par(mfrow=c(2,2))
	boxplot(res.m,col=8,names=snoise,main="resid. means",ylim=range(c(0,res.m)))
	abline(h=0,lty=2,col=8)
	boxplot(res.v,col=8,names=snoise,main="resid. vars",ylim=range(c(0,res.v)))
	abline(h=0,lty=2,col=8)
	boxplot(rss0,col=8,names=snoise,main="RSS0")
	boxplot(rss,col=8,names=snoise,main="RSS1",ylim=range(c(0,rss)))
	abline(h=0,lty=2,col=8)
	
	quartz()
	par(mfrow=c(3,3))
	boxplot(sp.a.m,col=8,names=snoise,main="a-spline bias",ylim=range(c(0,c(sp.a.m)),na.rm=T))
	abline(h=0,lty=2,col=8)
	boxplot(sp.b.m,col=8,names=snoise,main="b-spline bias",ylim=range(c(0,c(sp.b.m)),na.rm=T))
	abline(h=0,lty=2,col=8)
	boxplot(sp.t.m,col=8,names=snoise,main="tau-spline bias",ylim=range(c(0,c(sp.t.m)),na.rm=T))
	abline(h=0,lty=2,col=8)
	boxplot(sp.a.v,col=8,names=snoise,main="a-spline error var",ylim=range(c(0,c(sp.a.v))))
	abline(h=0,lty=2,col=8)
	boxplot(sp.b.v,col=8,names=snoise,main="b-spline error var",ylim=range(c(0,c(sp.b.v))))
	abline(h=0,lty=2,col=8)
	boxplot(sp.t.v,col=8,names=snoise,main="tau-spline error var",ylim=range(c(0,c(sp.t.v))))
	abline(h=0,lty=2,col=8)
	boxplot(sp.a.mse,col=8,names=snoise,main="a-spline MSE",ylim=range(c(0,c(sp.a.mse))))
	abline(h=0,lty=2,col=8)
	boxplot(sp.b.mse,col=8,names=snoise,main="b-spline MSE",ylim=range(c(0,c(sp.b.mse))))
	abline(h=0,lty=2,col=8)
	boxplot(sp.t.mse,col=8,names=snoise,main="tau-spline MSE",ylim=range(c(0,c(sp.t.mse))))
	abline(h=0,lty=2,col=8)
} 
if(0){
	# MC stats
	b.a = apply(thetas.mc[,a.inds]-thT[a.inds],2,mean)
	b.b = apply(thetas.mc[,b.inds]-thT[b.inds],2,mean)
	b.t = apply(thetas.mc[,tau.inds]-thT[tau.inds],2,mean)
	v.a = apply(thetas.mc[,a.inds],2,var)
	v.b = apply(thetas.mc[,b.inds],2,var)
	v.t = apply(thetas.mc[,tau.inds],2,var)
	m.a = apply((thetas.mc[,a.inds]-thT[a.inds])^2,2,mean)
	m.b = apply((thetas.mc[,b.inds]-thT[b.inds])^2,2,mean)
	m.t = apply((thetas.mc[,tau.inds]-thT[tau.inds])^2,2,mean)
	stats = data.frame(
		bias.a=b.a, bias.b=b.b, bias.tau=b.t,
		var.a=v.a, var.b=v.b, var.tau=v.t,
		mse.a=m.a, mse.b=m.b, mse.tau=m.t)
	# round(stats,3)
}
if(!doMC){
	quartz()
	par(mfrow=c(3,3), font.lab=2, font.axis=2, font=2)
	##
	ths = a.inds
	matplot(ths,cbind(thT[ths],thI[ths],thF[ths]),pch=c(1,15,20),cex=c(1,.8,.8),col=c(1,4,2),main="a",t='b')
	abline(h=0,col=8)
	ths = b.inds
	matplot(ths,cbind(thT[ths],thI[ths],thF[ths]),pch=c(1,15,20),cex=c(1,.8,.8),col=c(1,4,2),main="b",t='b')
	abline(h=0,col=8)
	ths = tau.inds
	matplot(ths,cbind(thT[ths],thI[ths],thF[ths]),pch=c(1,15,20),cex=c(1,.8,.8),col=c(1,4,2),main="tau",t='b')
	abline(h=0,col=8)
	legend("bottomright",pch=20,col=2,legend="final")
	##
	aspT = bspT = tspT = numeric(n)
	aspI = bspI = tspI = numeric(n)
	aspF = bspF = tspF = numeric(n)
	for(id in 1:n){
		aspT[id] = ifelse(modcon$axs==2,sum(x2s[id,]*thT[a.inds]),sum(x1s[id,]*thT[a.inds]))
		aspI[id] = ifelse(modcon$axs==2,sum(x2s[id,]*thI[a.inds]),sum(x1s[id,]*thI[a.inds]))
		aspF[id] = ifelse(modcon$axs==2,sum(x2s[id,]*thF[a.inds]),sum(x1s[id,]*thF[a.inds]))
		bspT[id] = ifelse(modcon$bxs==2,sum(x2s[id,]*thT[b.inds]),sum(x1s[id,]*thT[b.inds]))
		bspI[id] = ifelse(modcon$bxs==2,sum(x2s[id,]*thI[b.inds]),sum(x1s[id,]*thI[b.inds]))
		bspF[id] = ifelse(modcon$bxs==2,sum(x2s[id,]*thF[b.inds]),sum(x1s[id,]*thF[b.inds]))
		tspT[id] = sum(x1s[id,]*thT[tau.inds])
		tspI[id] = sum(x1s[id,]*thI[tau.inds])
		tspF[id] = sum(x1s[id,]*thF[tau.inds])
	}
	plot(aspT,ylim=range(c(aspI,aspF,aspT)),t='b',pch=1,main="a-sp")
	points(aspI,col=4,pch=20,cex=.8,t='b')
	points(aspF,col=2,pch=20,cex=.8,t='b')
	abline(h=0,col=8)
	plot(bspT,ylim=range(c(bspI,bspF,bspT)),t='b',pch=1,main="b-sp")
	points(bspI,col=4,pch=20,cex=.8,t='b')
	points(bspF,col=2,pch=20,cex=.8,t='b')
	abline(h=0,col=8)
	plot(tspT,ylim=range(c(tspI,tspF,tspT)),t='b',pch=1,main="tau-sp")
	points(tspI,col=4,pch=20,cex=.8,t='b')
	points(tspF,col=2,pch=20,cex=.8,t='b')
	abline(h=0,col=8)
	##
	res0 = z-zhat0
	res = z-out$zhat
	plot(z,main="true and fitted",ylim=range(c(res,z,out$zhat)))
	points(out$zhat,pch=20,cex=.8,col=2)
	points(res,pch=20,col=8,cex=.8)
	abline(h=0,col=3)
	boxplot(cbind(res0,res),main="residuals",col=c(4,8))
	abline(h=0,col=8)
	plot(out$M1s,main="RSS",t='b',pch=20)
}
if(modcon$gamval){
	roi = cbind(z,rep(1,nrow(xx)),xx)
	voi = rasterize.voi(roi,def=NA)
	roiF = cbind(out$zhat,rep(1,nrow(xx)),xx)
	voiF = rasterize.voi(roiF,def=NA)
	rr = range(c(roi[,1],roiF[,1]))
	quartz()
	slcs = c(1,6,7,12)
	par(mfcol=c(2,4))
	for(k in slcs){
		mimage(t(voi[,,k]), rr, main=paste(k))
		mimage(t(voiF[,,k]), rr, main=paste(k))
	}
}
