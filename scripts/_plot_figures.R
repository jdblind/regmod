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
source("funs/image_toolbox.R")
source("./optim_routine_exsarcroi.R")

fpath = "output_figs/"

# ************************************************************************************************ data plots
gam00 = 0.1
sn = 0.25
load(paste("out_simdata/synthetic_onerun_gam_",gam00,"_sn_",sn,".Rdata",sep=""))

pdf(file=paste(fpath,"simdata.pdf",sep=""),  width=12, height=5)
roi = cbind(zexact,rep(1,nrow(xx)),xx)
igrid = sort(unique(xx[,2]))
jgrid = sort(unique(xx[,1]))
kgrid = sort(unique(xx[,3]))
voi = rasterize.voi(roi,def=NA)
rr = range(c(roi[,1]))
sc = 1
layout(matrix(c(1,1,1,1,2,2,2,3,3,3,
				1,1,1,1,2,2,2,3,3,3,
				1,1,1,1,2,2,2,3,3,3,
				1,1,1,1,2,2,2,3,3,3), nr=4, byrow=T))
par(mar=c(2,2,2,2),
	font=2, font.lab=2, font.axis=2, cex=sc, cex.lab=sc, cex.axis=sc)
dd = dim(voi)
kslc = round(dd[3]/2,0)
jslc = round(dd[2]/2,0)
islc = round(dd[2]/2,0)
dd = 1/dd
lw = 3
midslice=voi[,,kslc]
mimage((voi[,,kslc]), rr, main=paste("Transverse"))
	abline(h=islc*dd[1]-dd[1]/2,col="yellow",lwd=lw)
	abline(v=jslc*dd[2]-dd[2]/2,col="yellow",lwd=lw)
mimage((voi[,jslc,]), rr, main=paste("Sagittal"))
	abline(h=kslc*dd[3]-dd[3]/2,col="yellow",lwd=lw)
	abline(v=islc*dd[1]-dd[1]/2,col="yellow",lwd=lw)
mimage((voi[islc,,]), rr, main=paste("Coronal"))
	abline(h=kslc*dd[3]-dd[3]/2,col="yellow",lwd=lw)
	abline(v=jslc*dd[2]-dd[2]/2,col="yellow",lwd=lw)
dev.off()

# picture of noise effect on mid-T slice 
pdf(file=paste(fpath,"noisy.pdf",sep=""),  width=10, height=5)
sc = .8
par(mfrow=c(1,2),
	mar=c(2,3,2,2),font=2,font.axis=2,font.lab=2,cex=1,cex.axis=sc,cex.lab=sc)
roi = cbind(zexact,rep(1,nrow(xx)),xx)
voi = rasterize.voi(roi,def=NA)
slcs = 6
snoise = c(1,.75,.5,.25)
set.seed(1)
sn = 0.5
zk = zexact + voi.noise(zexact,m=mnoise,s=sn)
roik = cbind(zk,rep(1,nrow(xx)),xx)
voik = rasterize.voi(roik,def=NA)
hist(zexact,col=8,xlab="",main="Uptake",xlim=range(c(zexact,zk)))
hist(zk,col=rgb(.5,0,1,alpha=.5),main="",xlab="",add=T)
hist(zexact,col=gray(.5,alpha=.5),xlab="",main="",add=T)
legend("topleft",legend=c("true","noisy"),col=c(gray(.25,alpha=.5),rgb(.5,0,1,alpha=.5)),pch=15)
rr = range(c(c(midslice),c(voik[,,slcs])))
mimage(t(voik[,,slcs]), rr, main=expression(paste(s[epsilon],"=0.5")))
dev.off()

# ************************************************************************************************ one-run example

theta.init=th0
theta.final=thF
source("scripts/wiremesh.R")
require(TeachingDemos)

pdf(file=paste(fpath,"OneRun_s025_gam01_RSS_p.pdf",sep=""), width=8, height=8)
sc=1.2
layout(matrix(c(1,1,2,2,3,3,
				1,1,2,2,3,3,
				1,1,2,2,3,3,
				0,4,4,4,4,0,
				0,4,4,4,4,0,
				0,4,4,4,4,0,
				0,4,4,4,4,0),nr=7,byrow=T))
par(mar=c(3,2,1,2),oma=c(0,0,0,0)+.1,
	font=2,font.axis=2,font.lab=2,cex=1,cex.axis=sc,cex.lab=sc)
##
xv=c(0,1,0);rrx=1.1
viewxv(xv,cbind(bxt[,1],bxt[,2],bxt[,3]),cbind(nxt[,1],nxt[,2],nxt[,3]),
		"x3'","x1'","True volume",rrx) 
viewxv(xv,cbind(bgam0[,1],bgam0[,2],bgam0[,3]),cbind(ngam0[,1],ngam0[,2],ngam0[,3]),
		"x3'","x1'","Initial volume",rrx) 
viewxv(xv,cbind(bgam1[,1],bgam1[,2],bgam1[,3]),cbind(ngam1[,1],ngam1[,2],ngam1[,3]),
		"x3'","x1'","Output volume",rrx) 
##
sc = 1
res0 = z-zhat0
res = z-out$zhat
boxplot(cbind(res0,res),ylab="Distributions of residuals",
	col=c(4,8),names=c("initial","final"))
abline(h=0,col=1)
ssc = sc-.1
subplot(
	plot(out$M1s,main="RSS",t='b',pch=20,ylab="",xlab="Iteration"),
	x=grconvertX(c(0.35,1),from='npc'),
	y=grconvertY(c(0.45,1),from='npc'),
	type='fig',	
	pars=list(mar=c(c(4,4,1,1)-0.1),cex.lab=ssc,cex.axis=ssc)
)
dev.off()

# ************************************************************************************************ data plots
gam00 = 0.1
sn = 0.25
load(paste("out_simdata/synthetic_onerun_gam_",gam00,"_sn_",sn,".Rdata",sep=""))

pdf(file=paste(fpath,"simdata.pdf",sep=""),  width=17, height=10)
roi = cbind(zexact,rep(1,nrow(xx)),xx)
igrid = sort(unique(xx[,2]))
jgrid = sort(unique(xx[,1]))
kgrid = sort(unique(xx[,3]))
voi = rasterize.voi(roi,def=NA)
rr = range(c(roi[,1]))
sc = 1
layout(matrix(c(0,0,2,2,2,2,3,3,4,4,
				1,1,2,2,2,2,3,3,4,4,
				1,1,2,2,2,2,3,3,4,4,
				1,1,2,2,2,2,3,3,4,4,
				5,5,6,6,7,7,8,8,9,9,
				5,5,6,6,7,7,8,8,9,9
				), nr=6, byrow=T))
par(mar=c(2,2,2,2),
	font=2, font.lab=2, font.axis=2, cex=sc, cex.lab=sc, cex.axis=sc)
hist(zexact,col=8,main="",xlab="Uptake")
dd = dim(voi)
kslc = round(dd[3]/2,0)
jslc = round(dd[2]/2,0)
islc = round(dd[2]/2,0)
dd = 1/dd
lw = 3
mimage((voi[,,kslc]), rr, main=paste("Transverse"))
	abline(h=islc*dd[1]-dd[1]/2,col="yellow",lwd=lw)
	abline(v=jslc*dd[2]-dd[2]/2,col="yellow",lwd=lw)
mimage((voi[,jslc,]), rr, main=paste("Sagittal"))
	abline(h=kslc*dd[3]-dd[3]/2,col="yellow",lwd=lw)
	abline(v=islc*dd[1]-dd[1]/2,col="yellow",lwd=lw)
mimage((voi[islc,,]), rr, main=paste("Coronal"))
	abline(h=kslc*dd[3]-dd[3]/2,col="yellow",lwd=lw)
	abline(v=jslc*dd[2]-dd[2]/2,col="yellow",lwd=lw)
# picture of noise effect on mid-T slice 
sc = .8
roi = cbind(zexact,rep(1,nrow(xx)),xx)
voi = rasterize.voi(roi,def=NA)
slcs = 6
rr = c(-1.8,9.6)
snoise = c(1,.75,.5,.25)
mimage(t(voi[,,slcs]), rr, main=paste("0"))
for(k in 4:1){
	set.seed(1)
	zk = zexact + voi.noise(zexact,m=mnoise,s=snoise[k])
	roik = cbind(zk,rep(1,nrow(xx)),xx)
	voik = rasterize.voi(roik,def=NA)
	mimage(t(voik[,,slcs]), rr, main=paste(paste(snoise[k])))
	print(c(min(c(voik)),max(c(voik))))
}
dev.off()

# ************************************************************************************************ one-run example

theta.true=thT
theta.init=th0
theta.final=thF
source("scripts/wiremesh.R")
require(TeachingDemos)
pdf(file=paste(fpath,"OneRun_s025_gam01_RSS.pdf",sep=""), width=16, height=6)
sc=1.2
# with true volume meshing:
layout(matrix(c(0,0,0,0,0,0,4,4,4,4,
				1,1,2,2,3,3,4,4,4,4,
				1,1,2,2,3,3,4,4,4,4),nr=3,byrow=T))
par(mar=c(2,3,2,2),font=2,font.axis=2,font.lab=2,cex=1,cex.axis=sc,cex.lab=sc)
##
xv=c(0,1,0);rrx=1.1
viewxv(xv,cbind(bxt[,1],bxt[,2],bxt[,3]),cbind(nxt[,1],nxt[,2],nxt[,3]),
		"x3'","x1'","True volume",rrx) 
viewxv(xv,cbind(bgam0[,1],bgam0[,2],bgam0[,3]),cbind(ngam0[,1],ngam0[,2],ngam0[,3]),
		"x3'","x1'","Initial volume",rrx) 
viewxv(xv,cbind(bgam1[,1],bgam1[,2],bgam1[,3]),cbind(ngam1[,1],ngam1[,2],ngam1[,3]),
		"x3'","x1'","Output volume",rrx) 
##
sc = 1
res0 = z-zhat0
res = z-out$zhat
boxplot(cbind(res0,res),ylab="Distributions of residuals",
	col=c(4,8),names=c("initial","final"))
abline(h=0,col=1)
ssc = sc-.1
subplot(
	plot(out$M1s,main="RSS",t='b',pch=20,ylab="",xlab="Iteration"),
	x=grconvertX(c(0.45,1),from='npc'),
	y=grconvertY(c(0.45,1),from='npc'),
	type='fig',	
	pars=list(mar=c(c(4,4,1,1)-0.1),cex.lab=ssc,cex.axis=ssc)
)
dev.off()

# ************************************************************************************************ Monte Carlo
rm(list=ls(pos=.GlobalEnv), pos=.GlobalEnv)
fpath = "output_figs/"

pdf(file=paste(fpath,"rss_vs_s_up.pdf",sep=""), width=10, height=8)
sc = 1.2
par(font=2,font.axis=2,font.lab=2,cex=sc,cex.axis=sc,cex.lab=sc,
		mar=c(4,4,3,3),oma=c(0,0,0,0))
###
gammas = c(0,0.1,.2)
rbox = array(NA,dim=c(100,5,length(gammas)))
for(gamk in c(1:length(gammas))){
	gam00 = gammas[gamk]
	load(paste("out_simdata/synthetic_noiseless_gam_",gam00,".Rdata",sep=""))
	out0 = out
	orss0 = out$M1
	load(paste("out_simdata/synthetic_run_gam_",gam00,".Rdata",sep=""))
	rbox[,,gamk] = cbind(rss, rep(orss0,100))
}
rbox=rbox/(n-length(c(a.inds,b.inds,tau.inds)))
#
cs = c(8,2,4)
gk = 2

meds = matrix(NA,nr=5,nc=length(gammas))
for(gk in c(1:length(gammas))){
	meds[,gk] = apply(rbox[,,gk],2,"median")
}
atk = c(-.25,0,.25)
bx = .4
boxplot(rbox[,,gk], col=cs[gk], ylim=range(c(rbox)), main="", 
	names=c(snoise,0), boxwex=bx, xlab=expression(s[epsilon]), ylab="MSE")
for(gk in c(1:length(gammas))){
	lines((c(1:5)+atk[gk]),meds[,gk],col=cs[gk],lwd=1.2,lty=2)
	boxplot(rbox[,,gk], col=cs[gk], add=T, names=NA, at=c(1:5)+atk[gk], xaxt="n", boxwex=bx)
}
abline(v=c(.5,1.5,2.5,3.5,4.5,5.5),lty=3,col=gray(c(1,1,1)/10,alpha=.3))
abline(h=0)

rect(grconvertX(0.45, from='npc'), grconvertY(0.45, from='npc'),
     grconvertX(.99, from='npc'), grconvertY(1, from='npc'), 
     col="white", border=NA)

rx = log(rbox[,1:4,])
meds = matrix(NA,nr=4,nc=length(gammas))
for(gk in c(1:length(gammas))){
	meds[,gk] = apply(rx[,,gk],2,"median")
}
xlocs = log(snoise)
cs = c(1,2,4)
atk = c(-.25,0,.25)*0
bx = .1
gk = 2
ssc=1
MSEline = lm(meds[,1]~xlocs)
subplot(
	{boxplot(rx[,,gk], ylim=range(c(rx)), col=cs[gk], main="", border=cs[gk],	
		at=xlocs, names=c(snoise), xlim=rev(range(c(xlocs))-c(.1,-.1)),
		bg="white",
		boxwex=bx, xlab=expression(paste("log(",s[epsilon],")")), ylab="log-MSE");
	abline(a=0,b=2,lwd=2,col=1,lty=1);
		lines((xlocs+atk[1]),meds[,1],col=cs[1],lwd=1.5,lty=2);
		boxplot(rx[,,1], col=cs[1], add=T, names=NA, border=cs[1],
			at=xlocs+atk[1], xaxt="n", boxwex=bx);
		lines((xlocs+atk[2]),meds[,2],col=cs[2],lwd=1.5,lty=2);
		boxplot(rx[,,2], col=cs[2], add=T, names=NA, border=cs[2],
			at=xlocs+atk[2], xaxt="n", boxwex=bx);
		lines((xlocs+atk[3]),meds[,3],col=cs[3],lwd=1.5,lty=2);
		boxplot(rx[,,3], col=cs[3], add=T, names=NA, border=cs[3],
			at=xlocs+atk[3], xaxt="n", boxwex=bx);
	},
	x=grconvertX(c(0.45,1),from='npc'),
	y=grconvertY(c(0.45,1),from='npc'),
	type='fig',	
	pars=list(mar=c(c(4,4,1,1)-0.1),cex.lab=ssc,cex.axis=ssc)
)
dev.off()

if(0){ ### MSE rate of CV to s_noise
quartz()
sns = c(snoise)
ins = c(1:K)
gk = 1
px = log(sns)
py = log(meds[1:4,])
plot(px,log(meds[ins,gk]),t='n',ylim=range(c(py)),xlim=rev(range(px)))
for(gk in c(1:3)){
	points(px,py[ins,gk],pch=20,t='b',col=gk)
}
legend("bottomleft",legend=paste("gamma =",gammas),pch=20,col=c(1:3))
}

# ************************************************************************************************ parameters
rm(list=ls(pos=.GlobalEnv), pos=.GlobalEnv)
fpath = "output_figs/"

# 1-study example of consistency
gam = 0.1
load(paste("out_simdata/synthetic_run_gam_",gam,".Rdata",sep=""))
# compute spline stats
sc = .9
quartz()
par(mfcol=c(3,2),font=2,font.axis=2,font.lab=2,cex=sc,cex.axis=sc,cex.lab=sc)
##
tis = which.min(rss[,4])
thF = c(thetas.mc[tis,,4])
ths = a.inds
matplot(ths,cbind(thT[ths],thI[ths],thF[ths]),
	pch=c(1,15,20),cex=c(1,.8,.8),col=c(1,4,2),
	main="a",xlab="",ylab="",t='b')
abline(h=0,col=8)
ths = b.inds
matplot(ths,cbind(thT[ths],thI[ths],thF[ths]),pch=c(1,15,20),cex=c(1,.8,.8),col=c(1,4,2),main="b",xlab="",ylab="",t='b')
abline(h=0,col=8)
ths = tau.inds
matplot(ths,cbind(thT[ths],thI[ths],thF[ths]),pch=c(1,15,20),cex=c(1,.8,.8),col=c(1,4,2),main="tau",xlab="",ylab="",t='b')
abline(h=0,col=8)
legend("bottomright",pch=20,col=2,legend="final")
###
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

# MC boxplots
sc = 1
pdf(file=paste(fpath,"MC_spline_boxplots_gam01_vs_gam05.pdf",sep=""), width=10, height=12)
# quartz(width=10,height=9)
par(mfrow=c(3,3),font=2,font.axis=2,font.lab=2,cex=sc,cex.axis=sc,cex.lab=sc)
gammas = c(0,0.1,0.5)
a.box = b.box = t.box = array(NA,dim=c(100,5,length(gammas)))
arange = c(0,0.6)
brange = c(0,1.4)
trange = c(0,0.12)
for(gamk in c(1:length(gammas))){
	gam00 = gammas[gamk]
	# load stats for noiseless case
	load(paste("out_simdata/synthetic_noiseless_gam_",gam00,".Rdata",sep=""))
	out0 = out
	amse0 = sp.a.mse
	bmse0 = sp.b.mse
	tmse0 = sp.t.mse
	am0 = sp.a.m
	bm0 = sp.b.m
	tm0 = sp.t.m
	avar0 = sp.a.v
	bvar0 = sp.b.v
	tvar0 = sp.t.v
	# aggregate with MC stats
	snoise = c(snoise,0.)
	load(paste("out_simdata/synthetic_run_gam_",gam00,".Rdata",sep=""))
	sp.a.mse = cbind(sp.a.mse, rep(amse0,M))
	sp.b.mse = cbind(sp.b.mse, rep(bmse0,M))
	sp.t.mse = cbind(sp.t.mse, rep(tmse0,M))
	sp.a.m = cbind(sp.a.m, rep(am0,M))
	sp.b.m = cbind(sp.b.m, rep(bm0,M))
	sp.t.m = cbind(sp.t.m, rep(tm0,M))
	sp.a.v = cbind(sp.a.v, rep(avar0,M))
	sp.b.v = cbind(sp.b.v, rep(bvar0,M))
	sp.t.v = cbind(sp.t.v, rep(tvar0,M))
	res0 = z-out0$zhat
	res.m = cbind(res.m, rep(mean(res0),M))
	res.v = cbind(res.v, rep(var(res0),M))
	rss0 = cbind(rss0, rep(out0$M0,M))
	rss = cbind(rss, rep(out0$M1,M))
	# create object for plot
	a.box[,,gamk] = sp.a.mse
	b.box[,,gamk] = sp.b.mse
	t.box[,,gamk] = sp.t.mse
	
	###
	if(gam00==0){
		atit=expression(paste("a, ",gamma,"=0"))
		btit=expression(paste("b, ",gamma,"=0"))
		ttit=expression(paste(tau,", ",gamma,"=0"))
	}
	if(gam00==0.1){
		atit=expression(paste("a, ",gamma,"=0.1"))
		btit=expression(paste("b, ",gamma,"=0.1"))
		ttit=expression(paste(tau,", ",gamma,"=0.1"))
	}
	if(gam00==0.5){
		atit=expression(paste("a, ",gamma,"=0.5"))
		btit=expression(paste("b, ",gamma,"=0.5"))
		ttit=expression(paste(tau,", ",gamma,"=0.5"))
	}
	boxplot(sp.a.mse,col=8,main=atit,names=c(snoise,0),ylim=arange,xlab=expression(s[epsilon]))
	abline(h=0,lwd=2,col=8)
	boxplot(sp.b.mse,col=8,main=btit,names=c(snoise,0),ylim=brange,xlab=expression(s[epsilon]))
	abline(h=0,lwd=2,col=8)
	boxplot(sp.t.mse,col=8,main=ttit,names=c(snoise,0),ylim=trange,xlab=expression(s[epsilon]))
	abline(h=0,lwd=2,col=8)
}
dev.off()
