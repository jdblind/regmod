#! /usr/bin/Rscript --vanilla
while( dev.next()>1 ){ dev.off() }
rm(list=ls(pos=.GlobalEnv), pos=.GlobalEnv)

library(splines)
library(gplots)
library(ic.infer)
library(Matrix)
library(plot3D)
library(inline)
library(fields)
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
source("funs/mapping_tools.R")

#------------------------------------------------------------------------------------
# (NB: cf. original case.plot() fun in ../opt2015)
case.plot <- function(slc.t,dopdf=FALSE,cutoff=cutoff){
# ptid: patient ID (global variable)
	# sample transverse slice
	if(!length(slc.t)){slc.t=round(dim(roi.zT)[3]/2,0)}

	if(dopdf){
		pdf(file=paste(ptid,"_dashboard_",slc.t,".pdf",sep=""),
			width=10,height=7)
	}
	###
	par(pch=20,lwd=2,font=2,font.lab=2,font.axis=2,font.main=2)
	cx.main = 2.0
	cx.ax = 1.5
	cx.lab = 1.3
	layout(
		matrix(c(1,1,2,2,3,3,
				  1,1,2,2,3,3,
				  1,1,2,2,3,3,
				  4,4,5,5,6,6,
				  4,4,5,5,6,6,
  				  4,4,5,5,6,6),
			nc=6,byr=T)
	)
	# --- data and fits
		
	# observed uptake
	im.int.suv = make.interp.im(roi.zT,sl=slc.t)
	image(im.int.suv,col=cg,axes=F,
		main=paste(ptid," (T, N=",nrow(xt),")",sep=""),
		cex.axis=cx.ax,cex.lab=cx.lab,cex.main=cx.main)

	v.im = cphg.a # any of cphg.ap, cphg.g, cphg.a, cphg.p
	v.im = cphg.e # any of cphg.ap, cphg.g, cphg.a, cphg.p
	qs = map.quantiles(v.im)
	ni = length(qs)
	qmap = make.qmap(v.im,qs=qs)
	map = make.interp.im(qmap,slc.t)
	mapcols = make.mapcols(v.im,qs,a=.5)
	image(im.int.suv,col=cg,axes=F,
		main="iso-quantile mapping",
		cex.axis=cx.ax,cex.lab=cx.lab,cex.main=cx.main)
	image(map,col=mapcols,add=T)
	image(im.int.suv,col=gray(c(0:255)/255,alpha=.95),add=T)
	image(map,col=mapcols,add=T)
	#
	pie(qs,col=mapcols,labels=round(qs,2))

	# --- 3D rendering
	tit = paste("cut-off =",round(cutoff,1),sep=" ")
	chut = render.meshing(cutoff,tit=tit,
		cex.axis=cx.ax,cex.lab=cx.lab,cex.main=cx.main)
	###
	if(dopdf){
		dev.off()
	}
}

# ---------------------------------------------------------------------------------
cp = .75
dopdf = T

itgt = 1
for(itgt in c(1:2)){
	ptid = c("Case_1","Case_2")[itgt]
	dfile = paste("data/",ptid,"/_doit_large_VOI.Rdata",sep="")
	load(dfile)
	cutoff = as.numeric(quantile(z,cp,na.rm=TRUE))

	slc.t = which(sort(unique(xx[,3]))==0)
	slc.t = 14
	
	source("utils/impro.R")
	source("funs/S_rfuns.R")
	source("scripts/script_summaries.R")
	source("scripts/script_summaries_gradients.R")
	
	case.plot(slc.t, dopdf=dopdf, cutoff=cutoff)
	
	v = c(cphg.g)
	(P90 = quantile(v,.9,na.rm=TRUE))
	summary(v)
}
