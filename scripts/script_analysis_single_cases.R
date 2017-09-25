#! /usr/bin/Rscript --vanilla
while( dev.next()>1 ){ dev.off() }
rm(list=ls(pos=.GlobalEnv), pos=.GlobalEnv)

library(splines)
library(gplots)
library(ic.infer)
library(Matrix)
library(plot3D)
library(inline)

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
source("funs/image_toolbox.R")
source("funs/devfuns.R")

# ------------------------------------------------------------------------ load data

itgt = 2
ptid = c("Case_1","Case_2")[itgt]
outdir = paste("./data/",ptid,"/out_an/",sep="")
if(!file.exists(outdir)){dir.create(outdir)}
dfile = paste("data/",ptid,"/_doit_large_VOI.Rdata",sep="")
load(dfile)
ptdir = paste(outdir,"gammet",modcon$gammet,"_",sep="")
	
# ------------------------------------------------------------------------ model summaries

source("scripts/script_summaries.R")
source("scripts/script_summaries_gradients.R")	

# ------------------------------------------------------------------------ plots...
if(0){
	xcheck = matrix(scan(file="./data/Case_1/tumor_x3.tsv",skip=8),nc=5,byrow=TRUE)
	# nrow(xcheck)
	# nrow(xx)
	# head(xx)
	# head(xcheck[,3:5])
	xs = sort.roi(xcheck)
	zs = unique(xs[,5])
	pdf(file=paste(ptdir,"Tslices.pdf",sep=""),width=15,height=10)
	cgs = c(gray(c(0:255)/255))
	par(mfrow=c(3,5))
	tim = rasterize.voi(xs)
	for(k in 1:length(zs)){
		slc = xs[which(xs[,5]==zs[k]),]
		image(t(tim[,,k]),axes=F,col=cgs,main=paste(k," (x3=",round(zs[k],2),")",sep=""))
	}
	dev.off()
	
	doquartz=FALSE
	pdf(file=paste(ptdir,"sampling.pdf",sep=""),width=15,height=10)
	source("scripts/script_pas_sampling_analysis.R")
	dev.off()
}

pdf(file=paste(ptdir,"3D_meshings.pdf",sep=""),width=15,height=10)
source("scripts/script_viewxv_gam.R")
dev.off()

# ------------------------------------------------------------------------ profile and convergence checks
cg = gray(c(0:255)/255)
slc.c = round(length(unique(xx[,1]))/2,0) # coronal slice
slc.s = round(length(unique(xx[,2]))/2,0) # sagittal slice
slc.t = round(length(unique(xx[,3]))/2,0) # transverse slice
slc.t = which(sort(unique(xx[,3]))==0) # transverse slice
slcs = round(c(slc.c,slc.s,slc.t),0)
dims = dim(roi.zT)

# quartz()
pdf(file=paste(ptdir,"execution.pdf",sep=""),width=15,height=10)
par(mfrow=c(2,3),pch=20,lwd=2,font=2)
plot(M0s,t='b',xlab="iter",ylab="RSS",main=paste("RSS incl. init,",ptid))
plot(M1s,t='b',xlab="iter",ylab="RSS",main="RSS w/out init")
plot(-(M01s),t='b',xlab="iter",ylab="RSS",main="%-change in RSS")
abline(h=0,col="grey",lwd=1)
#
plot(ps1.r,gh1,main="output std fit vs phases")
plot(gams,t='b',xlab="iter",ylab="gamma",main="Regularisation parameter")
plot(z-zh0.r,col=2,cex=.8,ylab="unweighted residuals",
	ylim=range(c(z-zh0.r,z-zh1.r)),main="Residual errors")
points(z-zh1.r)
legend("topleft",col=c(2,1),pch=20,legend=c("init","final"))
abline(h=0,col=3)
dev.off()

# ------------------------------------------------------------------------ load data

# images of uptakes:
# quartz()
pdf(file=paste(ptdir,"uptakes.pdf",sep=""),width=15,height=15)
par(mfrow=c(3,3))
rr = rr[!is.na(rr)]
mimage(roi.zT[slc.c,,],rr,main="coronal - true uptake");
crosshair(1,slcs,dims)
mimage(roi.zT[,slc.s,],rr,main="sagittal - true uptake")
crosshair(2,slcs,dims)
mimage(roi.zT[,,slc.t],rr,main="transverse - true uptake")
crosshair(3,slcs,dims)
mimage(roi.z0[slc.c,,],rr,main="coronal - initial uptake")
crosshair(1,slcs,dims)
mimage(roi.z0[,slc.s,],rr,main="sagittal - initial uptake")
crosshair(2,slcs,dims)
mimage(roi.z0[,,slc.t],rr,main="transverse - initial uptake")
crosshair(3,slcs,dims)
mimage(roi.z1[slc.c,,],rr,main="coronal - fitted uptake")
crosshair(1,slcs,dims)
mimage(roi.z1[,slc.s,],rr,main="sagittal - fitted uptake")
crosshair(2,slcs,dims)
mimage(roi.z1[,,slc.t],rr,main="transverse - fitted uptake")
crosshair(3,slcs,dims)
dev.off()

# quartz()
pdf(file=paste(ptdir,"phase_keys.pdf",sep=""),width=15,height=15)
dims = dim(roi.zT)
rr = c(range(roi.pTa),range(roi.pmask))
par(mfrow=c(3,3))
#
mimage(roi.pTa[slc.c,,],rr,main="coronal - init phase key");
crosshair(1,slcs,dims)
mimage(roi.pTa[,slc.s,],rr,main="sagittal - init phase key")
crosshair(2,slcs,dims)
mimage(roi.pTa[,,slc.t],rr,main="transverse - init phase key")
crosshair(3,slcs,dims)
#
mimage(roi.pmask[slc.c,,],rr,main="coronal - fitted phase key");
crosshair(1,slcs,dims)
mimage(roi.pmask[,slc.s,],rr,main="sagittal - fitted phase key")
crosshair(2,slcs,dims)
mimage(roi.pmask[,,slc.t],rr,main="transverse - fitted phase key")
crosshair(3,slcs,dims)
#
image(roi.zT[slc.c,,],col=cg,axes=F,main="coronal - true uptake");
crosshair(1,slcs,dims)
image(roi.zT[,slc.s,],col=cg,axes=F,main="sagittal - true uptake")
crosshair(2,slcs,dims)
image(roi.zT[,,slc.t],col=cg,axes=F,main="transverse - true uptake")
crosshair(3,slcs,dims)
dev.off()

# ------------------------------------------------------------------------ quantitative analysis

# quartz()
pdf(file=paste(ptdir,"summary.pdf",sep=""),width=20,height=10)
par(mfcol=c(2,5),pch=20,lwd=2,font=2,font.lab=2,font.axis=2)
#
# sample transverse slice
image(roi.zT[,,slc.t],col=cg,axes=F,
	main=paste("true (T, ",ptid,", N=",nrow(xt),")",sep=""))
#
# overall profile
plot(gp1$phases,gp1$ghat,axes=F,
	xlab="phases",ylab="",main="uptake profile")
axis(1)
#
# gamma selection
plot(gams,t='b',main="gamma (regul. par.)",xlab="iteration",ylab="gamma")
#
# z vs zhat
zhat = zh1.r
rr = range(c(z,zhat))
plot(z,zhat,xlim=rr,ylim=rr,
	main=paste("fitted vs true (",round(median(100*(z-zhat)/z),1),"% median err)",sep=""))
abline(a=0,b=1,lwd=2,col='orange')
#
# radial uptake summary
image(t(-uptk),col=cg,axes=F,
	xlab="",ylab="spine",
	main="radial uptake")
axis(1,at=c(.005,.995),labels=c("C","B"))
#
# radial phase summary
image(t(-abs(2-phs)),col=cg,axes=F,
	xlab="",ylab="spine",
	main="radial phase keys")
axis(1,at=c(.005,.995),labels=c("C","B"))
# 
# hets
plot(hets1,hk,xlim=c(0,100),t='b',pch=20,
				main=paste("Het (",round(HET,1),")",sep=""),xlab="Heterogeneity (%)",
				ylab="transverse slice coordinate")
# 
# phases
plot(mp,c(1:Jh),pch=20,t='b',main="phases",
	xlim=c(0,8),xlab="phase",
	ylab="transverse slice coordinate");
points(th1[tau.inds],c(1:Jh),pch=20,cex=.5,t='b',col=2)
legend("topright",col=c(1,2),legend=c("median","core"),pch=20,lty=1)
#
# transverse phase image
image(roi.pmask[,,slc.t],col=cg,axes=F,
	main="estimated phase keys (transverse slice)")
# image((-abs(2-roi.pmask))[,,slc.t],col=cg,axes=F,
	# main="estimated phase keys (transverse slice)")
#
image(cph[,,slc.t],col=cg,axes=F,main="a+phase classification (raw)")
dev.off()

# --------------
pdf(file=paste(ptdir,"progress_map.pdf",sep=""),width=15,height=10)
par(mfrow=c(2,4))
plot(gp,gh,main="profile")
plot(gp,-grp,pch=20,main="neg. grads"); abline(h=0)
plot(gp,aas*(-grp),main="coupled a_i*ng_i"); abline(h=0)
#
image(cph.abs[,,slc.t],col=cg,axes=F,main="a*phase classification (raw)")
#
image(cphg[,,slc.t],col=cg,axes=F,main="a*grad classification (raw)")

cols = c(pal.cold(nln),pal.hot(nlh))
image(cphg[,,slc.t],col=cols,axes=F,main="a+grad classification (raw)")
# pie(rep(1, nl), labels=mids, col = pal(nl))
if(length(hbn)){
	pie(rep(1, nln), labels=mids[which(mids<0)], col = pal.cold(nln), main="palette for negative areas")
} else {
	plot(rep(1,nlh),rep(1,nlh),t="n",axes=F,xlab="",ylab="")
	legend("center",legend="No negative values")
}
pie(rep(1, nlh), labels=mids[which(mids>=0)], col = pal.hot(nlh), main="palette for positive areas")
dev.off()

ptids = c(ptids, ptid)
hets = c(hets, HET)
medcphs = c(medcphs, median(theta1[tau.inds]))
medphs = c(medphs, median(gp1$phases))
minphs = c(minphs, min(gp1$phases))
medcphgs = c(medcphgs, mean(c(cphg)))
q90cphgs = c(q90cphgs, quantile(c(cphg),.9))
vols = c(vols, nrow(xx))
vols.voi = c(vols.voi, nrow(xr))
