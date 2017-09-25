#! /usr/bin/Rscript --vanilla
while( dev.next()>1 ){ dev.off() }
rm(list=ls(pos=.GlobalEnv), pos=.GlobalEnv)

# Paper plots
dsrc = "./"
setwd("./figs")
library(TeachingDemos)
library(splines)
library(fields)
source(paste(dsrc,"funs/image_toolbox.R",sep=""))
source(paste(dsrc,"funs/devfuns.R",sep=""))
load(file=paste(dsrc,"data/Case_1/_doit_large_VOI.Rdata",sep=""))
load(file=paste(dsrc,"data/Case_1/output_extractroi/eroi_init.Rdata",sep=""))
source(paste(dsrc,"utils/impro.R",sep=""))
source(paste(dsrc,"funs/mapping_tools.R",sep=""))
source(paste(dsrc,"scripts/script_summaries.R",sep=""))
source(paste(dsrc,"scripts/script_summaries_gradients.R",sep=""))

# ------------------------------------------------------------------------------------------

make.m.img  <- function(im.int.suv,im.x,im.s,m="",a=.5,...){
# In: script_summaries_gradients.R
	mcols1 = c(rgb(0,0,0,alpha=.05),rgb(1,0,0,alpha=a))
	mcols2 = c(rgb(0,0,0,alpha=.05),rgb(0,0,1,alpha=a))
	image(im.int.suv,col=cg,axes=F,
		main=m,
		cex.axis=cx.ax,cex.lab=cx.lab,cex.main=cx.main,...)
	image(im.x,col=mcols1,add=T)
	image(im.int.suv,col=gray(c(0:255)/255,alpha=.5),add=T)
	image(im.s,col=mcols1,add=T)
	image(im.x,col=mcols2,add=T)
}

# ------------------------------------------------------------------------------------------

# FIG: profile curves with gradient arrows
pdf("Case1_profile_and_gradients_new.pdf",height=6,width=6,compress=TRUE)
# quartz(height=6,width=6)
layout(matrix(c(1,1,1,1,1,1,2,2,2,2,2,2,
			1,1,1,1,1,1,2,2,2,2,2,2,
			1,1,1,1,1,1,2,2,2,2,2,2,
			3,3,3,3,3,3,4,4,4,4,4,4,
			3,3,3,3,3,3,4,4,4,4,4,4,
			3,3,3,3,3,3,4,4,4,4,4,4),nr=6,byr=T))
cx.main = 1.
cx.ax = 1.1
cx.lab = 1.1
par(font=2,font.lab=2,font.axis=2,
	cex=cx.main,cex.lab=cx.ax,cex.lab=cx.lab,
	mgp=c(2,0.5,0),
	mar=c(c(3,2,2.5,1.5)),oma=c(c(0,0,0,0)))
as.out = rasterize.voi(cbind(aas,c(1+0*z),xx))
v.e = gp.e*g.im.e  # weighting by voxel phase
v.e = as.out*g.im.e # weighting by voxel amplitude
p = c(gp.e)
gh = c(gh.e)
wg = c(g.im.e)
vg = c(v.e)
# order phases, g-curve, and gradients
is = order(p)
p.s = p[is]
g.s = gh[is]
w.s = wg[is]
v.s = vg[is]
#
# (a) uptake profile curve
# pick a point on curve and compute gradient arrow
plot(p.s,g.s,t='l',lwd=2.5,
	ylim=range(c(g.s,w.s),na.rm=T),
	axes=F,
	main="Uptake profile",
	xlab="phases",ylab="")
axis(1); axis(2)
is1 = 20
is = is1
a = p.s[is]+1
b = g.s[is]-w.s[is]
clr = 4
points(cbind(p.s[is],g.s[is]),pch=20,cex=1,col=clr)
arrow.plot(matrix(c(p.s[is],g.s[is]),nc=2),	
			matrix(c(a,b),nc=2), 
			true.angle=F,
			arrow.ex=abs(v.s[is]), length=.1, col=clr, lwd=3)
is2 = 1030
clr = 2
is = is2
a = p.s[is]+1
b = g.s[is]-w.s[is]
points(cbind(p.s[is],g.s[is]),pch=20,cex=1,col=clr)
arrow.plot(p.s[is],g.s[is],	
			u=1,v=-w.s[is], 
			arrow.ex=abs(v.s[is]), length=.1, col=clr, lwd=3)
#
# (b) gradient profile curve
points(p.s,w.s,t='l',col=8,lwd=4)#,
	# axes=F,
	# main="Gradient profile",
	# xlab="phases",ylab="Gradient")
abline(h=0,col=1,lwd=2,lty=1)
axis(1); axis(2)
points(cbind(p.s[is1],w.s[is1]),pch=20,cex=2,col=4)
points(cbind(p.s[is2],w.s[is2]),pch=20,cex=2,col=2)
#
# (c) voxel gradient values
plot(p.s,v.s,t='l',lwd=3,
	axes=F,
	main="Amplitude-weighted\n gradient profile",
	xlab="phases",ylab="Weighted gradient")
abline(h=0,col=1,lwd=2,lty=1)
points(cbind(p.s[is1],v.s[is1]),pch=20,cex=2,col=4)
points(cbind(p.s[is2],v.s[is2]),pch=20,cex=2,col=2)
axis(1); axis(2)
# (d) transverse slices - examples
# suv+g transverse slice (Case 1 sl 6)
slc.t = 7
im.int.suv = make.interp.im(roi.zT,slc.t)
im.s = make.interp.im(s.mask,slc.t)
im.g = make.interp.im(v.mask.g,slc.t)
make.m.img(im.int.suv,im.g,im.s,m="",xlab=paste("(T",slc.t,") raw gradients",sep=""))
# suv+v transverse slice (Case 1 sl 6)
im.ap = make.interp.im(v.mask.a,slc.t)
make.m.img(im.int.suv,im.ap,im.s,m="",xlab=paste("(T",slc.t,") v values",sep=""))	
dev.off()

# FIG: initial profile and phase conversion
pdf("input_profile.pdf",height=6,width=7)
sc=1.1
par(font=2,font.lab=2,font.axis=2,cex=sc,cex.lab=sc,cex.lab=sc,
		mar=c(c(4,0,4,0)))
yn=eroi$efit[,3] 
wy=eroi$efit[,2]
wy=wy/max(wy) 
wcut=sort(wy)[.4*length(wy)]
if(wcut==max(wy)){ 
	wcut=sort(wy)[.2*length(wy)]
}
i.s = which(wy>wcut)
v.s = v[i.s]; y.s = yn[i.s]
i.s = which(!is.na(y.s))	
ylimi = range(c(ghat[i.s]/g2g$fac,g2g$gam.fit))
plot(v,ghat/g2g$fac,ylim=ylimi,main="Initial input profile",
	xlim=c(-1.05,2.1),
	xlab="[<- Core]     Voxel phase    [Boundary ->]",
	ylab="",
	axes=F,pch=20,col=8,cex=.6)
is = order(v)
points(v[is],ghat2[is]/g2g$fac,col="navy",t='l',lwd=6)
axis(1)
text(x=-1.1,y=0.4,srt=90,labels="Uptake")
###
subplot(
	plot(v[is],g2g$gv[is],pch=20,cex=.6,xlab="Initial input phases",ylab="Final input phases"),
		x=grconvertX(c(0.75,1),from='npc'),
		y=grconvertY(c(0.6,1),from='npc'),
		type='plt',
		pars=list(mar=c(1.5,1.5,0,0)+0.5,cex=.8))
dev.off()

# FIG: initial analysis (het char. and vol rendering)
source("funs/_script_viewxv_gam.R")
pdf("het_and_mesh.pdf",height=4,width=10)
sc=1.1
sc.m=1
par(mfrow=c(1,3),
	font=2,font.lab=2,font.axis=2,
	mgp=c(1.8,0.5,0),
	mar=c(c(3,3,2.5,1.5)), oma=c(c(0,0,0,0)),
	cex=sc.m,cex.lab=sc,cex.axis=sc)
# uptake slice:
rr = c(range(roi.zT,na.rm=TRUE))
cg = gray(c(0:255)/255)
slc.c = round(length(unique(xx[,1]))/2,0) # coronal slice
slc.s = round(length(unique(xx[,2]))/2,0) # sagittal slice
slc.t = round(length(unique(xx[,3]))/2,0) # transverse slice
slc.t = which(sort(unique(xx[,3]))==0) # transverse slice
slcs = round(c(slc.c,slc.s,slc.t),0)
eepars = eroi$epars
hets0 = 100-100*eepars[,3]
# Overall het measures:
# initial tubular measure:
(HET0 = sum(eepars[,10]*(1-eepars[,3])*100)/sum(eepars[,10]))
# output tubular measure:
res1 = z-gp1$zhat
(HET1 = 100 - 100*max(0, 1-(mean(res1^2)/var(z))))
quantile(hets0,.75)
plot(hets0,eepars[,1],xlim=c(0,100),pch=20,t='b',
		main=paste("Heterogeneity (H0=",round(HET0,1),", H1=",round(HET1,1),")",sep=""),
		axes=F,
		# xlab="1-R^2",
		xlab="Transverse heterogeneity",
		ylab="Spinal location (mm)")
axis(1); axis(2) 
# output (regularized?) spinal het
points(hets1,eepars[,1], t='b', pch=15, col=2)
legend("right",legend=c("input","output"),col=c(1,2),pch=c(20,15))
# 3d volume rendering:
getwd()
# initial input volume
xv=c(0,1,0);rrx=1.1
viewxv(xv,cbind(bxt[,1],bxt[,2],bxt[,3]),cbind(nxt[,1],nxt[,2],nxt[,3]),
		"x3'","x1'","Initial volume",rrx) 
# final Gamma volume
xv=c(0,1,0);rrx=1.1
viewxv(xv,cbind(bgam1[,1],bgam1[,2],bgam1[,3]),cbind(ngam1[,1],ngam1[,2],ngam1[,3]),
		"x3'","x1'","Regularized volume",rrx) 
dev.off()

# FIG: gradient mapping and weighted gradients v
pdf("demo_gradient_mapping.pdf",# height=5,width=12)
quartz(height=5,width=12)
par(mfrow=c(1,3),font=2,font.lab=2,font.axis=2)
cx.main = 2.0
cx.ax = 1.8
cx.lab = 1.5
# suv+g transverse slice (Case 1 sl 6)
slc.t = 7
im.int.suv = make.interp.im(roi.zT,slc.t)
im.s = make.interp.im(s.mask,slc.t)
im.g = make.interp.im(v.mask.g,slc.t)
make.m.img(im.int.suv,im.g,im.s,m=paste("(T",slc.t,") raw gradients",sep=""))
# suv+v transverse slice (Case 1 sl 6)
im.ap = make.interp.im(v.mask.a,slc.t)
make.m.img(im.int.suv,im.ap,im.s,m=paste("(T",slc.t,") v values",sep=""))
# v-mapping coronal/sagittal slice	
slc.s = 15
im.int.suv = make.interp.slice(t(roi.zT[slc.s,,]))
im.s = make.interp.slice(t(s.mask[slc.s,,]))
im.ap = make.interp.slice(t(v.mask.a[slc.s,,]))
make.m.img(im.int.suv,im.ap,im.s,m=paste("(S",slc.s,") v values",sep=""))
dev.off()
	
