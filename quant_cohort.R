#! /usr/bin/Rscript --vanilla
while( dev.next()>1 ){ dev.off() }
rm(list=ls(pos=.GlobalEnv), pos=.GlobalEnv)

library(splines)
library(ic.infer)
library(Matrix)
library(plot3D)
# library(spatstat)
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

# ------------------------------------------------------------------------ inits...
paste(R.Version()$major,R.Version()$minor)
print(Sys.time())
dosynthetic=FALSE

MUTE.G = FALSE
pcut=.25
alpha=xalpha=5
# nb=Jh=26 
# nres=Jphi=25 
dohull=FALSE

if(!file.exists("./out")){
	dir.create("./out")
} 

if(!file.exists("./data/outR")){
	dir.create("./data/outR")
} 

dirlist = dir("./data/alldatasetsJ")
if(sum(is.element(dirlist,"output_extractroi"))){
	dirlist = dirlist[-1]
}

gam0 = 1
modcon = list(astep=TRUE,axs=2,bxs=2,dorep=0,
			gammet = 1, gamval = gam0, lm.method = 2,
			# gam.Om = .1,
			pe.min = 1, # minimum percent error to be achieved
			maxit = 10, nbloops.b = 1, nbloops.tau = 1)
if(modcon$gammet<0){
	patch = paste("_gam",modcon$gamval,"_",sep="")
} else {
	patch = paste("_gammet",modcon$gammet,sep="")
}

sizes = read.csv("data/sizes.csv")
is = order(sizes[order(sizes[,2]),1]) # proceed by increasing size 

for(i in is){
	modcon$gamval = gam0
	print("------------------------------------------")
	ptid = dirlist[i]

	# read in the data
	dirroot = paste("./data/outE/output_extractroi_",ptid,sep="")
	xr = matrix(scan(file=paste("./data/alldatasetsJ/",ptid,
		"/tumor-pre.tsv",sep=""),skip=8),nc=5,byrow=TRUE)
	load(file=paste(dirroot,"/eroi_init.Rdata",sep=""))
	thres = scan(paste(dirroot,"/thres.txt",sep=""))
	nb = Jh = eroi$Jh
	nres = Jphi = eroi$Jphi
	xp = matrix(scan(file=paste(dirroot,"/ROIN.txt",sep="")),nc=5,byrow=TRUE)
	# read in larger VOI and adjust for mu's (as done for xp)
	xt = matrix(scan(file=paste(dirroot,"/ROIT.txt",sep="")),nc=5,byrow=TRUE)
	head(xt) # (w,z,x1,x2,x3)
	xt = xt[,c(1,3:5,2)] # (w,x1,x2,x3,z)
	ab = hab(xt)
	zb = eroi$z.pasbins
	for(k in 1:nb){
		ii = which((xt[,4]<=zb[k+1])&(xt[,4]>zb[k]))
		if(length(ii)){
			xt[ii,2] = xt[ii,2]-as.numeric(eroi$mu1[k])
			xt[ii,3] = xt[ii,3]-as.numeric(eroi$mu2[k])
		}
	}
	if(0){ ### run once
		print("running spline.Omega()...")		
		Omega = spline.omega(ab)
		print("done.")		
		save(Omega,file=paste(dirroot,"/splineOmega.Rdata",sep=""))
	}
	load(file=paste(dirroot,"/splineOmega.Rdata",sep=""))

	# -----------------------------------------------------------------------------------
	# info from PAS domain...
	xr.init = xr
	xa = xr[eroi$icut,]
	nrow(xr)
	nrow(xt)
	nrow(xa)
	nrow(xp)
	
	nx = length(unique(xa[,4]))
	ny = length(unique(xa[,3]))
	nz = length(unique(xa[,5]))
	
	i0 = which(xa[,1]<=0)
	xa[i0,1] = abs(.001*rnorm(length(i0)))
	
	eroi$s = 1+0*eroi$s
	z = xr[,1]
	hv = xr[,5]
	w = xr[,2]
	xxo = xr[,3:5]
	xx = xxo
	xt = xt[,2:4]

	# ------------------------------------------------------------------------ VOI features
	xx.o = xx
	z.o = z
	xxh = xx

	print(paste("N_xx =",nrow(xx)," ; N_xxh =",nrow(xxh)))
	gx = sort(unique(xx[,2]))
	gy = sort(unique(xx[,1]))
	gz = sort(unique(xx[,3]))
	gxh = sort(unique(xxh[,2]))
	gyh = sort(unique(xxh[,1]))
	gzh = sort(unique(xxh[,3]))
	
	nghbr.xx = list.neighborhood(xx,gx,gy,gz)
	if(nrow(xx)==nrow(xxh)){
		xth = xt
		xinds = c(1:nrow(xx))
		nghbr.xxh = nghbr.xx
	} else {
		xth = project.voi(xxh,eroi$c,eroi$s,eroi$xi) 
		xinds = find.xinds(xx,xxh)
		nghbr.xxh = list.neighborhood(xxh,gxh,gyh,gzh)
	}		
	
	# ------------------------------------------------------------------------ spline stuff
	phia = 0
	phib = 2*pi
	alpha = 3
	
	# -----------------------------------------------------------------------------------
	# Convert initial profile to Gamma profile and init phases
	source("scripts/script_convert_to_gamma.R")
	source("scripts/script_init_phases.R")
	
	# -----------------------------------------------------------------------------------
	# Initialize new parametrization using priors from extractroi()
	gyt = sort(unique(xt[,1]))
	gxt = sort(unique(xt[,2]))
	gzt = sort(unique(xt[,3]))
	ntx = length(gxt)
	nty = length(gyt)
	ntz = length(gzt)
	ab = range(rphs[,3])
	
	# Initial params
	ztrue = z
	
	ii = sort(c(c.inds,s.inds,xi.inds,mu1.inds,mu2.inds)) # optim wrt a+b+tau
	exinds = c(1:length(unlist(all.inds)))[ii]
	pis = c(1:length(unlist(all.inds)))[-ii]
	fixstruct = prod(is.element(c(c.inds,s.inds,xi.inds,mu1.inds,mu2.inds),exinds))
	if((length(pis)==length(a.inds))&&(sum(abs(pis-a.inds))==0)){modcon$maxit=1}
	
	theta.true = th0
	# Estimation bounds on theta
	thUL = set.bounds(th0,all.inds)#,method=0,dor=dorep)
	thL = thUL$lower
	thU = thUL$upper
	# reset bounds for non-optimised params
	thL[exinds] = th0[exinds]-1e-12
	thU[exinds] = th0[exinds]+1e-12
	thL[-exinds] = pmin(thL[-exinds],c(th0[-exinds]-1.5))
	thU[-exinds] = pmax(thU[-exinds],c(th0[-exinds]+1.5))
	# mu's
	mu.inds = sort(c(mu1.inds,mu2.inds))
	# a
	if(!mmatch(exinds,a.inds)){
		thL[a.inds] = 1e-2 #pmax(1e-3,theta.true[a.inds]-10) 
		thU[a.inds] = theta.true[a.inds]+10
	}
	# b
	if(!mmatch(exinds,b.inds)){
		thU[b.inds] = theta.true[b.inds]+4*abs(theta.true[b.inds])
		thL[b.inds] = pmax(1e-3,theta.true[b.inds]/10)
	}
	# tau
	if(!mmatch(exinds,tau.inds)){
		thL[tau.inds] = rep(.2,length(tau.inds))
		thU[tau.inds] = rep(5,length(tau.inds))
	}	
	th0 = pmax(thL,pmin(th0,thU))
	theta.true = th0 
	
	# LX.R = eval.lam.R(orph,th0[tau.inds],th0[a.inds],th0[b.inds],alpha)
	LX.C = eval.lam.C(orph,th0[tau.inds],th0[a.inds],th0[b.inds],alpha)
	# summary(LX.R-LX.C)
	LX = LX.C
	r1 = orph$rphs[xinds,1]
	
	# ----------------------------------------------------------------------------
	NPA = length(unlist(all.inds))-length(exinds) # nb of params being optimised
	(fixstruct = as.logical(prod(is.element(c(c.inds,s.inds,xi.inds,mu1.inds,mu2.inds),exinds))))

	# ----------------------------------------------------------------------------
	# run main optim code
	print(paste("------>",ptid,"......."))
	source("optim_routine_exsarcroi.R")
	save.image(paste("./data/outR/doit_",ptid,patch,".Rdata",sep=""))
	print(Sys.time())
}
