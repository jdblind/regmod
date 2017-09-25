# Includes:
# - rastering
# - (more conventional) image displaying
# - rotating
# - interpolating
# - image information metrics (entropy, variogram, some Laplacian, ...)
#

#------------------------------------------------------------------------------------------ GENERIC FUNS

as.real <- function(x){
# In: impro.R
# as.real being deprecated, this is just a wrapper :)
	as.double(x)
}

rad.to.deg <- function(th){
# In: impro.R
	return(th*180/pi) 
}

is.2D <- function(voi){
# In: impro.R
	return(length(dim(voi))==2)
}

is.3D <- function(voi){
# In: impro.R
	return(length(dim(voi))==3)
}

rebin <- function(x,xg,dd=diff(xg)[1]){
# In: impro.R
# Rebins x data according to REGULAR xg grid of width dd
	shft = abs(min(xg/dd))
	return((round(x/dd+shft)-shft)*dd)
}

fill.roi <- function(x_in){
# In: impro.R
# Input x_in is in list format.
# Output is also in list format, only filled with 0's... 
	ras = rasterize.voi(x_in)
	return(list.raster(ras))
}

#------------------------------------------------------------------------------------------ BASIC IMAGING

r.image <- function(img,col=gray(c(0:255)/255),...){
# In: impro.R
# Adapts R's basic image() function to perform conventional raster display. 
# Input img is a raster in matrix format:
#
# (0,0)
#	+-------------------> (j, x)
#	|
#	|		O	O
#	|	  	  U
#	|		\___/
#   |         
#	v
# (i,y)
#
# Usage example: r.image(img, main="toto")
#
	imr = t(img)[,c(nrow(img):1)]
	image(imr,axes=F,useRaster=T,col=col,...)
}	

r.points <- function(pts,...){
# In: impro.R
# Adapts R's basic points() function to add to r.image (using conventional raster display). 
# Input pts is a list of 2-column pixel coordinates.
# Usage example: 
#	r.image(imr, main="toto")
#	r.points(sbox,pch=15,col='red',cex=.8)
#
	if(nrow(pts)>1){
		points(pts[,c(2,1)],...)		
	} else {
		points(pts[c(2,1)],...)
	}
}	

draw.box <- function(lx,rx,ly,ry,nx,ny){
# In: impro.R
	out=rbind(	
		cbind(c(lx:rx),rep(ly,(rx-lx+1))),
		cbind(c(lx:rx),rep(ry,(rx-lx+1))),
		cbind(rep(c(lx),(ry-ly+1)),c(ly:ry)),
		cbind(rep(c(rx),(ry-ly+1)),c(ly:ry)))
	outl=cbind((c(0:nx)/nx)[out[,1]],(c(0:ny)/ny)[out[,2]])		
	return(outl)		
}

#------------------------------------------------------------------------------------------ RASTER FUNS

vec.to.raster <- function(x, nc, nr, dim=1){
# In: impro.R
# 
# Re-arranges a collection of (intensity) values into a raster format, i.e. either a 
# 2D or 3D matrix defining a VOI.
# Arguments:
# 	x: 	 a collection of values (vector or sequence)
#   nc:	 number of columns
#   nr:	 number of rows
#   dim: alternatively(?), use "dim" - as done in array()
# The 3rd dimension is worked out automatically based on nc and nr.
# Value:
#   returns the 2D or 3D volume.
#
# Cf. demo.rasterization.from.list() for demo.
#
	if(length(dim)>1){
		nr=dim[1]
		nc=dim[2]
	} 
	nf=ifelse(length(dim)==3,dim[3],length(x)/(nc*nr)) # nb of frames
	if(nc*nr==length(x)){ # 2D raster
		voi=matrix(x,ncol=nc,byrow=TRUE) # image of whole slice
	} else {
		if(floor(nf)==ceiling(nf)){
			# voi=array(x,dim=c(nr,nc,nf))  
			## problem with this line above: sort.roi() sorts ROIs according to 3rd col first,
			## then 4th col and then 5th col. This means the output is in such a way that
			## element need to be arranged into a matrix BY ROW... this option is not 
			## available for array()
			voi=array(0,dim=c(nr,nc,nf))
			for(kz in 1:nf){
				inds=c(1:(nc*nr))+(nc*nr)*(kz-1)
				voi[,,kz]=matrix(x[inds],ncol=nc,byrow=TRUE) # image of whole slice				
			}
		} else {
			stop("in vec.to.raster: non-integer number of frames requested...")
		}
	}
	return(voi)
}

raster.to.list <- function(rin,gx,gy,gz){
# In: impro.R
# Rearranges raster in list format (eg Amide form (uptake, w, gx, gy, gz)), with dummy weights w=1.
# Test:
#   rin=array(c(1:36),c(4,3,3))
#	gx=c(1:5)[-c(3,4)]
#	gy=c(11:15)[-2]
#	gz=c(21:24)[-3]	
#   raster.to.list(rin,x,y,z)
#
	nx=length(gx)
	ny=length(gy)
	nz=length(gz)
	if(prod(dim(rin))!=(nx*ny*nz)){
		stop("In raster.to.list: grid legnths do not match raster size, 
				i.e. prod(dim(rin))!=(nx*ny*nz)...")
	}
	liste = cbind(rep(x,ny*nz), rep(rep(y,each=nx),nz), rep(z,each=nx*ny)) # full list
	lout = NULL
	for(k in 1:nz){
		combo = cbind(c(t(rin[,,k])), rep(1,nx*ny), rep(x,ny), rep(y, each=nx), rep(z[k],nx*ny))
		lout = rbind(lout, combo)
	}
	return(lout)
}

rasterize.voi <- function(voi,def=0){
# In: impro.R
# Input:
#   voi is in list-format (u,w,x,y,z) as in AMIDE
# Value:
#   returns a raster form of the voi
# The output raster is filled with 0's where values 
# are not found in voi.
	xg = sort(unique(voi[,3]))
	yg = sort(unique(voi[,4]))
	zg = sort(unique(voi[,5]))
	NR = length(yg)
	NC = length(xg)
	NS = length(zg)
	ras = array(def,dim=c(NR,NC,NS))
	for(v in 1:nrow(voi)){
		i = which(yg==voi[v,4])
		j = which(xg==voi[v,3])
		k = which(zg==voi[v,5])
		ras[i,j,k] = voi[v,1]
	}
	return(ras)
}

rasterize.voi.loose <- function(voi,xl=NULL,yl=NULL,zl=NULL){
# In: impro.R
# Essentially the same as rasterize.voi, except it produces a 
# raster spanning a wider hypergrid specified by xl, yl,zl.
# Input:
#   voi is in list-format (u,w,x,y,z) as in AMIDE
# Value:
#   returns a raster form of the voi
# The output raster is filled with 0's where values 
# are not found in voi.
	xg = sort(unique(voi[,3]))
	yg = sort(unique(voi[,4]))
	zg = sort(unique(voi[,5]))
	dx = diff(xg)[1]
	dy = diff(yg)[1]
	dz = diff(zg)[1]
	# expand grid first
	if(length(xl)==2){
		xa = rev(seq(min(xg),xl[1],by=sign(xl[1]-min(xg))*dx)[-1])
		xb = seq(max(xg),xl[2],by=sign(xl[2]-max(xg))*dx)[-1]
		xg = c(xa,xg,xb)
	}
	if(length(yl)==2){
		ya = rev(seq(min(yg),yl[1],by=sign(yl[1]-min(yg))*dy)[-1])
		yb = seq(max(yg),yl[2],by=sign(yl[2]-max(yg))*dy)[-1]
		yg = c(ya,yg,yb)
	}
	if(length(zl)==2){
		za = rev(seq(min(zg),zl[1],by=sign(zl[1]-min(zg))*dz)[-1])
		zb = seq(max(zg),zl[2],by=sign(zl[2]-max(zg))*dz)[-1]
		zg = c(za,zg,zb)
	}
	# now carry on
	NR = length(yg)
	NC = length(xg)
	NS = length(zg)
	ras = array(0,dim=c(NR,NC,NS))
	for(v in 1:nrow(voi)){
		i = which(yg==voi[v,4])
		j = which(xg==voi[v,3])
		k = which(zg==voi[v,5])
		ras[i,j,k] = voi[v,1]
	}
	return(ras)
}

rasterize <- function(v,xx){
# In: impro.R
	return(rasterize.voi(cbind(v,0*v+1,xx)))
}

list.raster <- function(ras,xg,yg,zg){
# In: impro.R
# Transforms input raster ras into a list-format (u,w,x,y,z)
# as per AMiDE format (here w=1)
	NR = length(yg)
	NC = length(xg)
	NS = length(zg)
	u = c(ras)
	x = rep(xg,NR*NS)
	y = rep(rep(yg,each=NC),NS)
	z = rep(zg,each=NR*NC)
	xx = cbind(u,1+0*u,x,y,z)
	yy = subset(xx,xx[,1]>0)
	return(yy)
}

demo.rasterization.from.list <- function(){
# In: impro.R
# following Amide format: (suv, weight, x, y, z) 
# i.e. (suv, w, col, row, slice)...
	A = cbind(rnorm(12),rnorm(12),rep(c(1:3),4),rep(c(11:12),each=3),rep(c(8:9),each=6))
	A = rbind(A,c(rnorm(2),2,11,10))
	print("Initial region A (in list form) should have only 1 element at (1,2) in 3rd slice:")
	print(A)
	## So a rasterized A should display as:
	## , , 1
	##
	##           [,1]       [,2]     [,3]
	## [1,] 0.3130696 -0.1469504 1.073888
	## [2,] 0.5504397 -2.5140672 1.120984
	##
	## , , 2
	##
	##            [,1]      [,2]       [,3]
	## [1,]  2.5812700 0.3164413 -0.6493965
	## [2,] -0.4884614 0.6906856 -0.9952144
	##
	## , , 3
	##
	##   	 [,1]       [,2] [,3]
	## [1,]    0 -0.3265859    0
	## [2,]    0  0.0000000    0
	#
	print("Filling A out to form B using B = fill.roi(A):")
	B = fill.roi(A)
	print(B)
	print("Rasterizing B the WRONG WAY ROUND (rows and cols are switched):")
	print(array(B[,1],dim=c(2,3,3))) 		# not correct: rows and cols are switched
	print("Rasterizing B correctly is done using vec.to.raster(B[,1],nr=2,nc=3):")
	ok=vec.to.raster(B[,1],nr=2,nc=3) 	# correct: this is how we created A initially
	print(ok)
}

#------------------------------------------------------------------------------------------ INTERPOLATION

interpol.im <- function(im,xg,xd){
# In: impro.R	
# 1D interpolation of grid intensities...
# This is used following successive 1D image shifts.
# It could probably be replaced with the bilinear interpolation.
#	
	imo=im
	for(i in 1:ncol(im)){
		imo[,i]=approx(xg,im[,i],c(xg+xd),rule=2)$y
	}
	return(imo)
}

interp.raster <- function(rin,x,y,z){
# In: impro.R
# Interpolates input raster rin so that output has regular grids in (x,y,z) plan.
# For the moment, proceeds with linear interpolation in each direction (x,y,z) successively.
# Will upgrade to bilinear interpolation or some other type (kriging, NN, etc) later.
# Arguments:
#	rin: input raster (array format), where dim(rin) are the effective lengthes of the current grids x (columns), 
#		 y (rows) and z (slices).
# Value:
# 	rout: raster (array format) with regular grids where intermediate values from input x, y and z are interpolated.
# Test:
#   rin=array(c(1:36),c(4,3,3))
#	x=c(1:5)[-c(3,4)]
#	y=c(11:15)[-2]
#	z=c(21:24)[-3]
#	...
# 
	# just in case... and to avoid R bug in comparing floats: we blow the grid up by a factor 10
	x = round(sort(unique(x))*10)
	y = round(sort(unique(y))*10) 
	z = round(sort(unique(z))*10)
	dx = min(unique(diff(x)))
	dy = min(unique(diff(y)))
	dz = min(unique(diff(z)))
	dimin = dim(rin)
	# prepare output
	gx = seq(x[1],x[length(x)],by=dx)
	gy = seq(y[1],y[length(y)],by=dy)
	gz = seq(z[1],z[length(z)],by=dz)
	gdim = c(length(gy),length(gx),length(gz))
	rout1 = array(0,dim=c(dimin[1],gdim[2],dimin[3]))
	rout2 = array(0,dim=c(gdim[1],gdim[2],dimin[3]))
	rout3 = array(0,dim=c(gdim[1],gdim[2],gdim[3]))
	# interpolate raster...
	# interpolate new columns first 
	for(k in 1:dimin[3]){
		for(i in 1:dimin[1]){
			ligne = rin[i,,k]
			rout1[i,,k] = approx(x=x,y=ligne,xout=gx)$y
		}
	}
	# interpolate new rows then
	for(k in 1:dimin[3]){
		for(j in 1:gdim[2]){
			colonne = rout1[,j,k]
			rout2[,j,k] = approx(x=y,y=colonne,xout=gy)$y
		}
	}
	# and interpolate new slices last
	for(i in 1:gdim[1]){
		for(j in 1:gdim[2]){
			coupe = rout2[i,j,]
			rout3[i,j,] = approx(x=z,y=coupe,xout=gz)$y
		}
	}
	return(list(rout=rout3,gx=gx/10,gy=gy/10,gz=gz/10))
}

#------------------------------------------------------------------------------------------ RESIZING

shrink <- function(tim,xt,yt,gx=NULL,gy=NULL,meth="nearest",doc=TRUE,ratio=TRUE){
# In: impro.R
# Resizing (downsampling) of 2D images (transverse slices).
# It performs 2 univariate approximations using R's approx() function.
# This is therefoer NOT an optimal 2D spatial interpolation procedure.
# NB: this version does not require package adimpro. 
# Arguments:
# 	tim: 2D image - like a raster slice 
#	xt, yt: number of desired output x- and y-axis bins 
#			(resp. j and i array coordinates)
#			Note that if ratio=TRUE, these may be adjusted.
# 			NB: xt = dim(image)[2], yt = dim(image)[1] in the input; 
#				but this is switched in call to shrink.image()
#				since the latter expects the following:
#				xt = dim(image)[1], yt = dim(image)[2]
#	gx,gy: 	x- and y-grids required for interpolation.
#	meth: 	used to be used in version calling adimpro... Now obsolete.
#	doc: 	used to be used in version calling adimpro... Now obsolete.
#	ratio: 	used to be used in version calling adimpro... Now obsolete.
# Value:
#	output 2D array.
#
	if(is.null(gx)){gx = c(1:ncol(tim))}
	if(is.null(gy)){gy = c(1:nrow(tim))}
	tam = matrix(0,nr=nrow(tim),nc=xt)
	tom = matrix(0,nr=yt,nc=xt)
	gxi = seq(min(gx),max(gx),l=ncol(tim))
	gyi = seq(min(gy),max(gy),l=nrow(tam))
	xgi = seq(min(gxi),max(gxi),l=xt)
	ygi = seq(min(gyi),max(gyi),l=yt)
	for(i in 1:nrow(tim)){
		tam[i,] = approx(gxi,tim[i,],xo=xgi,rule=2)$y
	}
	for(j in 1:ncol(tam)){
		tom[,j] = approx(gyi,tam[,j],xo=ygi,rule=2)$y
	}
	return(tom)
}

shrink.adimpro <- function(ras,xt,yt,meth="nearest",doc=TRUE,ratio=TRUE){
# In: impro.R
# Resizing (downsampling) of 2D images (transverse slices).
# Requires package adimpro:
# 	http://cran.r-project.org/web/packages/adimpro/adimpro.pdf
# 	http://www.jstatsoft.org/v19/i01/paper
# Arguments:
# 	ras: 2D image - like a raster slice 
#	xt, yt: number of desired output x- and y-axis bins 
#			(resp. j and i array coordinates)
#			Note that if ratio=TRUE, these may be adjusted.
# 			NB: xt = dim(image)[2], yt = dim(image)[1] in the input; 
#				but this is switched in call to shrink.image()
#				since the latter expects the following:
#				xt = dim(image)[1], yt = dim(image)[2]
#	meth: 	adimpro::shrink.image()'s argument for interpolation method.
#		  	Can be either "nearest" (fastest), "median" or "mean). 		
#			"median" is supposed to give best results.
#	doc: 	adimpro::shrink.image()'s argument 'compress'... 
#			Whether or not to copmress the data (in intermediate objects
#			which are of class "adimpro")
#	ratio: 	adimpro::shrink.image()'s argument 'ratio' - whether to keep 
#			input image aspect ratio or not...
# Value:
#	output 2D array.
#
	require(adimpro)
	img = make.image(ras,compress=doc)
	omg = shrink.image(img,xt=yt,yt=xt,met=meth,compress=doc,ratio=ratio)
	return(extract.image(omg))
}

shrink.raster <- function(ras,xt,yt,zt,meth="nearest",doc=TRUE,ratio=FALSE){
# In: impro.R
# Resizing (downsampling) of 3D raster images (bunch of transverse slices).
# NB: requires package adimpro via shrink().
# Transverse slices are resized using shrink(), then slices are interpolated
# to desired output slice number 'zt' by linear interpolation along each 
# (i,j) coordinate.
# NB: xt = dim(image)[2], yt = dim(image)[1]
# NB: does not operate z-axis interpolation if zt >= dim(ras)[3].
#
	# (1) shrink transverse slices
	nz = dim(ras)[3]
	omg = array(0,dim=c(yt,xt,nz))
	for(k in 1:nz){
		omg[,,k] = shrink(ras[,,k],xt=xt,yt=yt,met=meth,doc=doc,ratio=ratio)
	}
	# (2) interpolate voxel lines in z-direction
	if(zt<nz){
		out = array(0,dim=c(yt,xt,zt))
		for(j in 1:xt){
			for(i in 1:yt){
				out[i,j,] = approx(c(1:nz),omg[i,j,],xo=seq(1,nz,l=zt),rule=2)$y
			}
		}
	} else {
		out = omg
	}
	return(out)
}

demo.shrink.raster <- function(){
# In: impro.R
# Testing...
	xr = matrix(scan(file="./data/ABRAMS/tumor.tsv",skip=8),nc=5,byrow=TRUE)
	xr = sort.roi(xr)
	# (1) rasterize.voi...
	ras = rasterize.voi(xr)
	dim(ras)
	gs = gray(c(0:255)/255)
	par(mfrow=c(2,4))
	for(k in 9:12){
		image(ras[,,k],col=gs,ax=F)
	}
	# (2) list.raster...
	xg = sort(unique(xr[,3]))
	yg = sort(unique(xr[,4]))
	zg = sort(unique(xr[,5]))
	lvoi = list.raster(ras,xg,yg,zg)
	nrow(xr)
	nrow(lvoi)
	# (3) rasterize.voi.loose...
	dx = diff(xg)[1]
	dy = diff(yg)[1]
	dz = diff(zg)[1]
	xl = c(min(xg)-dx,max(xg)+10*dx)
	yl = c(min(yg)-10*dy,max(yg)+2*dy)
	zl = NULL
	lras = rasterize.voi.loose(xr,xl,yl,zl)
	for(k in 9:12){
		image(lras[,,k],col=gs,ax=F)
	}
	
	quartz()
	par(mfcol=c(2,4))
	print(dim(ras))
	xl = 12
	yl = 14
	for(k in 9:12){
		image(ras[,,k],col=gs,ax=F)
		ro = shrink(ras[,,k],xt=xl,yt=yl)
		image(ro,col=gs,ax=F)
	}
	
	xl = 12
	yl = 14
	zl = 33
	tst = shrink(ras[,,1],xt=xl,yt=yl,ratio=FALSE)
	xl = dim(tst)[2]
	yl = dim(tst)[1]
	rbs = array(0,dim=c(yl,xl,zl))
	print(paste("Initial number of voxels:",prod(dim(ras))))
	print(paste("Final number of voxels:",prod(dim(rbs))))
	for(k in 1:dim(ras)[3]){
		rbs[,,k] = shrink(ras[,,k],xt=xl,yt=yl,ratio=FALSE)
	}
	
	quartz()
	par(mfcol=c(2,4))
	for(k in 9:12){
		image(ras[,,k],col=gs,ax=F)
		image(rbs[,,k],col=gs,ax=F)
	}
	
	rcs = shrink.raster(ras,xt=xl,yt=yl,zt=15)
	
	quartz()
	par(mfcol=c(3,5))
	for(k in 1:(dim(rcs)[3])){
		image(rcs[,,k],col=gs,ax=F,main=paste(k))
	}
}


#------------------------------------------------------------------------------------------ EXTRA

raster.contigs <- function(im){
# In: impro.R
# Identifies contiguous regions...
# Counts number of contiguous regions in a rasterized mask. 
# Input raster im may be 2D or 3D. 
# Values
# 	coordinates: 	coordinates of contiguous points; all regions are represented and 
# 					separated by a row of 0's
#	list: 			list of voxel IDs (as in ordered entries in the 2D or 3D array im)
#	number:		number of distinct contiguous regions found
#
# Tests
# 	# test matrix with 3 distinct contiguous regions
# 	A = array(0,dim=c(8,5,4))
# 	A[1,4:5,1] = A[1,4:5,2] = 1 		# 1st region
# 	A[1,1:2,4] = 1						# 2nd region
# 	A[3:4,2:5,1] = A[7,1:2,1] = 1		# 3rd region
# 	A[3:5,2:5,2] = A[7,1:2,2] = 1
# 	A[2,2,2] = 1
# 	A[7,1:2,3] = A[4:6,3:5,3] = 1
# 	A[6,2:4,3] = A[8,4:5,3] = 1
# 	A[3,2:3,3] = 1
# 	A[4,3:4,3] = 0
# 	A[6:8,2:5,4] = 1
# 	A[3,3:4,4] = 1
# 	A[4,2:5,4] = 1
# 	A
# 	# build raster of voxel IDs
# 	N=prod(dim(A))
# 	A.voxid = array(c(1:N),dim=dim(A))
#	# 	1er test... 2D
# 	tim=A[,,1]
# 	tim.id=A.voxid[,,1]
# 	raster.2D.contigs(tim)$coordinates
# 	raster.contigs(tim)$coordinates
# 	raster.2D.contigs(tim)$coordinates-raster.contigs(tim)$coordinates
# 	tom=tim
# 	tom[7,1]=0
# 	tom[1,]=0
# 	raster.2D.contigs(tom)$coordinates
# 	raster.contigs(tom)$coordinates
# 	raster.2D.contigs(tom)$coordinates-raster.contigs(tom)$coordinates
# 	tam=tim
# 	tam[3,4]=tam[4,3]=0
# 	tam
# 	raster.2D.contigs(tam)$coordinates
# 	raster.contigs(tam)$coordinates
# 	raster.2D.contigs(tam)$coordinates-raster.contigs(tam)$coordinates
# 	# 	2e test... 3D
# 	raster.3D.contigs(A)
# 	raster.contigs(A)
# 	raster.3D.contigs(A)$coordinates-raster.contigs(A)$coordinates
#
	dimim = dim(im)
	im.id = array(c(1:prod(dim(im))),dim=dimim)
	cands = which(im==1)
	cands.array = arrayInd(cands,.dim=dimim)
	contig.list = NULL
	contig.coords = NULL
	while(length(cands)){
		liste = cands[1]
		lind = 1 # index for swiping liste
		while(lind<=length(liste)){
			# identify target from candidates
			tgt.ind = which(cands==liste[lind])
			tgt = cands[tgt.ind]
			# get its coordinates and remove this target from candidates
			if(length(cands.array)>length(dimim)){
				tgt.coord = cands.array[tgt.ind,]	
				cands.array = cands.array[-tgt.ind,]
			} else {
				tgt.coord = cands.array
				cands.array = NULL
			}			
			cands = cands[-tgt.ind]
			# find its neighbours in raster
			imin = max(1,tgt.coord[1]-1)
			imax = min(dimim[1],tgt.coord[1]+1)
			jmin = max(1,tgt.coord[2]-1)
			jmax = min(dimim[2],tgt.coord[2]+1)
			if(length(dimim)==3){
				kmin = max(1,tgt.coord[3]-1)
				kmax = min(dimim[3],tgt.coord[3]+1)				
				neighbs = im.id[imin:imax,jmin:jmax,kmin:kmax][which(im[imin:imax,jmin:jmax,kmin:kmax]==1)]
			} else {
				neighbs = im.id[imin:imax,jmin:jmax][which(im[imin:imax,jmin:jmax]==1)]			
			}
			new.neighbs = neighbs[!is.element(neighbs,liste)]
			if(length(new.neighbs)){
				# add any new neighbours to contiguity list
				liste = c(liste, new.neighbs)
			}
			lind = lind+1 # update index for swiping liste				
		}		
		contig.list = c(contig.list, liste, 0)
		contig.coords = rbind(contig.coords, arrayInd(liste,.dim=dimim), c(rep(0,length(dimim))))
	}
	nb = sum(contig.list==0) # number of contiguous regions found
	return(list(coordinates=contig.coords, list=contig.list, number=nb))	
}

#------------------------------------------------------------------------------------------ polar rotation

im.rot.pol <- function(imin,phi){
# In: impro.R
# wrapper function!
	phi=-phi # conventions, conventions...
	im.rot.pol.bwd.interp(imin,phi)
}
#
im.rot.pol.fwd <- function(imin,phi){
# In: impro.R
# Matrix rotation implemented in polar representation: rotates input image imin 
# (in raster/matrix form) by angle phi (in rads).
# Returns rotated matrix object (forward rotation).
# Example usage:
#	image(t(imin),col=gray(c(0:255)/255),axes=F,main="orig",useRaster=T)
#	imout=im.rot.pol(imin,phi=0)
#	image(t(imout),col=gray(c(0:255)/255),axes=F,main="0",useRaster=T)
#
	imout=0*imin
	Ni=nrow(imin)
	Nj=ncol(imin)
	xss=seq(1:Nj)
	yss=seq(1:Ni)
	mx=Nj/2  #mean(xss)  
	my=Ni/2	 #mean(yss)
	for(i in 1:Ni){
		for(j in 1:Nj){
			# raster to cartesian (with (0,0) = image center)
			x=j-mx
			y=-i+my
			# cartesian to polar
			r=sqrt(x^2+y^2)
			th=atan2(y,x)
			if((x==0)&(y==0)){ # center
				imout[i,j]=imin[i,j]
			} else {
				if(x==0){
					if(y<0){th=pi*3/2}else{th=pi/2}
				}
				# throwing values from source to destination...
				th=th+phi
				# polar (via cartesian) to raster
				x=round(r*cos(th)+mx)
				y=round(-r*sin(th)+my)
				# check bounds here...
				if((x<1)|(x>Nj)|(y<1)|(y>Ni)){
					val=0
					x=min(max(x,1),Nj)
					y=min(max(y,1),Ni)	
				}else{val=imin[i,j]}
				imout[y,x]=val
			}
		}
	}
	return(imout)
}
###
im.rot.pol.bwd <- function(imin,phi){
# In: impro.R
# Matrix rotation implemented in polar representation: rotates input image imin 
# (in raster/matrix form) by angle phi (in rads).
# Returns rotated matrix object (backward rotation).
# Example usage:
#	image(t(imin),col=gray(c(0:255)/255),axes=F,main="orig",useRaster=T)
#	imout=im.rot.pol(imin,phi=0)
#	image(t(imout),col=gray(c(0:255)/255),axes=F,main="0",useRaster=T)
#
	imout=0*imin
	Ni=nrow(imin)
	Nj=ncol(imin)
	xss=seq(1:Nj)
	yss=seq(1:Ni)
	mx=Nj/2  #mean(xss)  
	my=Ni/2	 #mean(yss)
	for(i in 1:Ni){
		for(j in 1:Nj){
			# raster to cartesian (with (0,0) = image center)
			x=j-mx
			y=-i+my
			# cartesian to polar
			r=sqrt(x^2+y^2)
			th=atan2(y,x)		
			if((x==0)&(y==0)){ # center
				imout[i,j]=imin[i,j]
			} else {
				if(x==0){
					if(y<0){th=pi*3/2}else{th=pi/2}
				}
				# lifting values from source to destination...
				th=th-phi
				# polar (via cartesian) to raster
				x=round(r*cos(th)+mx)
				y=round(-r*sin(th)+my)
				# check bounds here...
				if((x<1)|(x>Nj)|(y<1)|(y>Ni)){val=0}
				else{val=imin[y,x]}
				imout[i,j]=val
			}
		}
	}
	return(imout)
}
###
im.rot.pol.bwd.interp <- function(imin,phi){
# In: impro.R
# Matrix rotation implemented in polar representation: rotates input image imin 
# (in raster/matrix form) by angle phi (in rads).
# Returns rotated matrix object (backward rotation). WITH BILINEAR INTERPOLATION.
# Example usage:
#	image(t(imin),col=gray(c(0:255)/255),axes=F,main="orig",useRaster=T)
#	imout=im.rot.pol(imin,phi=0)
#	image(t(imout),col=gray(c(0:255)/255),axes=F,main="0",useRaster=T)
#
	mimin=min(c(imin))
	Mimin=max(c(imin))
	imout=0*imin
	Ni=nrow(imin)
	Nj=ncol(imin)
	xss=seq(1:Nj)
	yss=seq(1:Ni)
	mx=Nj/2  #mean(xss)  
	my=Ni/2	 #mean(yss)
	for(i in 1:Ni){
		for(j in 1:Nj){
			# raster to cartesian (with (0,0) = image center)
			x=j-mx
			y=-i+my
			# cartesian to polar
			r=sqrt(x^2+y^2)
			th=atan2(y,x)
			if((x==0)&(y==0)){ # center
				imout[i,j]=imin[i,j]
			} else {
				if(x==0){
					if(y<0){th=pi*3/2}else{th=pi/2}
				}
				# lifting values from source to destination...
				th=th-phi
				# polar (via cartesian) to raster
				xt=r*cos(th)+mx
				yt=-r*sin(th)+my
				xt=min(max(xt,1),Nj)
				yt=min(max(yt,1),Ni)	
				# bilinear interpolation
				minx=floor(xt)
				miny=floor(yt)
				maxx=ceiling(xt)
				maxy=ceiling(yt)
				# check bounds here...
				if((minx<1)|(miny<1)|(maxx>Nj)|(maxy>Ni)){val=0}
				else{				
					dx=xt-minx
					dy=yt-miny
					# pixel neighborhood 
					tl=imin[miny,minx]
					tr=imin[maxy,minx]
					bl=imin[miny,maxx]
					br=imin[maxy,maxx]
					# 1st interpolation pass 
					iv=(1-dx)*tl+dx*tr
					iiv=(1-dx)*bl+dx*br
					# 2nd pass
					# fv=round((1-dy)*iv+dy*iiv)	
					fv=((1-dy)*iv+dy*iiv)
					# final interpolated intensity
					val=max(mimin,min(Mimin,fv))
				}
				imout[i,j]=val
			}
		}
	}
	return(imout)
}

#------------------------------------------------------------------------------------------ info quant

lap.het <- function(im){
# In: impro.R
# Some sort of cumulated image Laplacian to score edge dynamics (hence heterogeneity)...
# This version is reliable for 2D rasters im. For 3D rasters, it returns the mean of all
# frame scores.
# 
	nr=nrow(im)
	nc=ncol(im)
	lh=NULL
	if(length(dim(im))==3){ # then it's a 3D raster
		nf=dim(im)[3]
	} else { # 2D raster
		nf=1
		im=array(im,dim=c(dim(im)[1:2],nf))
	}
	for(k in 1:nf){
		lhk=0
		for(i in 2:(nr-1)){
			for(j in 2:(nc-1)){
				lhk = lhk + 6*im[i,j,k] - sum(im[i,c(j-1,j+1),k]) - sum(im[c(i-1,i+1),j,k])
			}
		}
		lh=c(lh,lhk)
	}
	return(mean(lh))
}

vario <- function(im){
# In: impro.R
# Implements a simple variogram on direct voxel pairs...
# This version is reliable for 2D rasters im. For 3D rasters, it returns the mean of all
# frame scores.
# 
	nr=nrow(im)
	nc=ncol(im)
	vo=NULL
	if(length(dim(im))==3){ # then it's a 3D raster
		nf=dim(im)[3]
	} else { # 2D raster
		nf=1
		im=array(im,dim=c(dim(im)[1:2],nf))
	}
	for(k in 1:nf){
		vok=0
		for(i in 1:nc){ # row-wise pairs
			x=im[,i,k]; 
			vo=vo+sum((x[-1]-x[-nr])^2)
		}
		for(i in 1:nr){ # col-wise pairs
			x=im[i,,k]; 
			vo=vo+sum((x[-1]-x[-nc])^2)
		}
		vok=vok/(nc*(nr-1)+nr*(nc-1))
		vo=c(vo,vok)
	}
	return(mean(vo))
}

al.entropy <- function(x,h=0,met=1){
# In: impro.R
# Ahmad-Lin 1D entropy estimator. h is the bandwidth for the cross-validated KDE. 
	n=length(x)
	if(h==0){ 
		xx=sort(x)
		sx=(xx[3*n/4]-xx[n/4])/1.3490  # not sd(x)
		if(met==1){	h=2.345*sx*n^(-1/5) }
		else{ h=1.664*sx*n^(-1/5) }
	}
	ent=0
	for(i in 1:n){
		xi=x[i]
		kdei=sum(exp(-(xi-x[-i])^2/(2*h^2)))/(sqrt(2*pi)*(n-1)*h)
		ent=ent+log(kdei)
	}
	return(-ent/n)
}

al.entropy.2d <- function(x,hi=0,hj=0,met=1){
# In: impro.R
# Ahmad-Lin 2D entropy estimator. h is the bandwidth for the cross-validated KDE. 
	nr=nrow(x)
	nc=ncol(x)
	if(hi==0|hj==0){ 
		xx=sort(c(x))
		sx=(xx[3*n/4]-xx[n/4])/1.3490  # not sd(x)
		if(met==1){	h=2.345*sx*n^(-1/5) }
		else{ h=1.664*sx*n^(-1/5) }
		hi=hj=h
	}
	ent=0
	for(i in 1:nc){
		for(j in 1:nr){
			xk=x[i,j]
			# 2D KDE = product of 2 univariate KDEs
			kdekj=sum(exp(-(xk-x[-i,j])^2/(2*hj^2)))/(sqrt(2*pi)*(n-1)*hj)
			kdeki=sum(exp(-(xk-x[i,-j])^2/(2*hi^2)))/(sqrt(2*pi)*(n-1)*hi)
			ent=ent+log(kdekj*kdeki)
		}
	}
	return(-ent/(nr*nc))
}
