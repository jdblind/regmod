source("funs/mapping_tools.R")

make.m.img <- function(im.int.suv,im.x,im.s,m="",a=.5){
# In: script_summaries_gradients.R
# requires source("tmi_figs/mapping_tools.R") be run in main script
	mcols1 = c(rgb(0,0,0,alpha=.05),rgb(1,0,0,alpha=a))
	mcols2 = c(rgb(0,0,0,alpha=.05),rgb(0,0,1,alpha=a))
	image(im.int.suv,col=cg,axes=F,
		main=m,
		cex.axis=cx.ax,cex.lab=cx.lab,cex.main=cx.main)
	image(im.x,col=mcols1,add=T)
	image(im.int.suv,col=gray(c(0:255)/255,alpha=.5),add=T)
	image(im.s,col=mcols1,add=T)
	image(im.x,col=mcols2,add=T)
}

rr = c(range(roi.zT),range(roi.z0),range(roi.z1))
cg = gray(c(0:255)/255)
slc.c = round(length(unique(xx[,1]))/2,0) # coronal slice
slc.s = round(length(unique(xx[,2]))/2,0) # sagittal slice
slc.t = round(length(unique(xx[,3]))/2,0) # transverse slice
slcs = round(c(slc.c,slc.s,slc.t),0)
dims = dim(roi.zT)

# Final het measure:
res1 = z-gp1$zhat
HET = 100 - 100*max(0, 1-(mean(res1^2)/var(z)))

# phase classification image
cph.abs = as.out*roi.pmask # aterm * phase_term raster
cph = as.out*roi.p1 # aterm * phase_term raster

# gradient of profile function - use negative of value - div diffs
ip = order(gp1$phases)
gp = gp1$phases[ip]
gh = gp1$ghat[ip]
grp = diff(gh)/diff(gp)
grp = c(grp,NA)
g.im = rasterize.voi(cbind(-grp,c(1+0*z),xx[ip,]),def=NA)
gp.im = rasterize.voi(cbind(gp,c(1+0*z),xx[ip,]),def=NA)
as.out = rasterize.voi(cbind(aas[ip],c(1+0*z),xx[ip,]))
as.out = as.out/max(as.out)
cphg = as.out*gp.im*g.im # aterm * grad_term raster

# weighted-gradient (v-values)
cphg.g = g.im 
cphg.p = gp.im*g.im 
cphg.a = as.out/max(as.out)*g.im 
cphg.ap = as.out/max(as.out)*gp.im*g.im
#
s.mask = make.qmask(roi.zT,q=.90)
v.mask.g = make.qmask(cphg.g,q=.90)
v.mask.p = make.qmask(cphg.p,q=.90)
v.mask.a = make.qmask(cphg.a,q=.90)
v.mask.ap = make.qmask(cphg.ap,q=.90)

# gradient of profile function - use negative of value - exact
gp = gp1$phases
gp.e = rasterize.voi(cbind(gp,c(1+0*z),xx),def=NA)
gh.e = rasterize.voi(cbind(gp1$ghat,c(1+0*z),xx),def=NA)
grp.e = -gp^(alpha-2)*exp(-gp)*(alpha-1-gp)
g.im.e = rasterize.voi(cbind(grp.e,c(1+0*z),xx),def=NA)
as.out = rasterize.voi(cbind(aas,c(1+0*z),xx))
cphg.e = as.out*g.im.e # aterm * grad_term raster # checked: matches div diffs version cphg

# midpoints and voxel labelling
hb = hist(c(cphg),plot=F)
hmids = hb$mids
hbn = hmids[hmids<0]
nln = 5
nlh = 8
if(length(hbn)){
	mids = seq(min(hbn),max(hbn),l= nln)
} else {
	mids = NULL
}
hbp = hmids[hmids>=0]
mids = c(mids, seq(min(hbp),max(hbp),l=nlh))
nl = length(mids)
