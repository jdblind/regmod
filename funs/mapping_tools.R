make.qmask <- function(im,q=.90){
# input image im must be in raster form
# e.g. im = cphg or im = roi.zT
	vim = im
	v.mask = 0*vim
	v.mask[vim> quantile(c(vim),q,na.rm=TRUE)] = 1
	return(v.mask)
}

make.qqmask <- function(im,q1=.90,q2=1,lab=1){
# input image im must be in raster form
# q1 = min-quantile, q2 = max-quantile
# mask indicates values within (q1,q2)
# e.g. im = cphg or im = roi.zT
	vim = im
	v.mask = 0*vim
	v.mask[vim> quantile(c(vim),q1,na.rm=TRUE)] = lab
	v.mask[vim> quantile(c(vim),q2,na.rm=TRUE)] = 0
	return(v.mask)
}

make.qmap <- function(im,qs=seq(.7,1,by=.1)){
	imf = 0*im
	for(ll in 1:(length(qs)-1)){
		mll = make.qqmask(im,q1=qs[ll],q2=qs[(ll+1)],lab=ll)
		imf = imf+mll
	}
	return(imf)
}

make.interp.im <- function(vmask,sl=round(c(dim(vmask)[3]/2),0),iby=.25){
# Interpolates transverse slice sl of raster vmask.
# Defaults:
#	- sl corresponds to the mid-height slice
#	- interp to the 1/4 pixel i.e. iby=.25
# Requires fields library.
	require(fields)
	im = t(vmask[,,sl])
	nr = nrow(im)
	nc = ncol(im)		
	im.xyz = list(y=seq(1,nc,by=1),x=seq(1,nr,by=1),z=im)
	im.int = interp.surface.grid(im.xyz,
		list(y=seq(1,nc,by=iby),x=seq(1,nr,by=iby)))
	return(im.int)
}

make.interp.slice <- function(im,iby=.25){
# Interpolates transverse/coronal/sagittal slice im
# Defaults:
#	- interp to the 1/4 pixel i.e. iby=.25
# Requires fields library.
	require(fields)
	nr = nrow(im)
	nc = ncol(im)		
	im.xyz = list(y=seq(1,nc,by=1),x=seq(1,nr,by=1),z=im)
	im.int = interp.surface.grid(im.xyz,
		list(y=seq(1,nc,by=iby),x=seq(1,nr,by=iby)))
	return(im.int)
}

map.quantiles <- function(im,Q1=.1,Q2=.6){
	# qs = seq(.5,1,by=.05)
	# or:
	# looking for max of neg quantiles
	qsearch = seq(Q1,1,by=.1)
	qpmin = rev(which(quantile(im,qsearch,na.rm=T)<0))[1]
	qs = seq(0,qsearch[max(qpmin,1,na.rm=T)],by=.1) # neg values
	qs = c(qs,seq(qsearch[max(qpmin,1,na.rm=T)+1],Q2,by=.25))
	qs = c(qs,seq(Q2,1,by=.05))
	return(unique(qs))	
}

make.mapcols <- function(im,qs,a=.5){
# qs = output of map.quantiles()
	ni = length(qs) # nb of quantiles used for isomasks
	if(ni==3){
		mapcols = c(rgb(0,0,0,alpha=.05), # 0 = black
				rgb(0,0,1,alpha=a),  # 1 = blue
				rgb(0,1,0,alpha=a),  # 2 = green
				rgb(1,0,0,alpha=a))  # 3 = red
	} else {# or, if using ni>3 labels:
		qvals = quantile(c(im),qs,na.rm=TRUE)
		# mapcols = c(0,rev(heat.colors(ni-1,alpha=a)))
		mapcols = c(rev(rainbow(ni,alpha=a)))
		# if any neg quantiles paint 1st color pink
		if(length(which(qvals<0))){ 
			mapcols[1] = col2hex("pink")		
		}
	}
	return(mapcols)	
}

make.mapcols.2pies <- function(im){
# palette control: nlc cold steps, nlhhot steps
	qim = quantile(im,qs,na.rm=T)
	nlc = sum(qim<0)
	nlh = length(qs)-nlc
	steps.cold <- (c("blue4", "blue2", "cyan"))
	pal.cold <- color.palette(steps.cold, c(2,1), space="rgb")
	steps.hot <- rev(c("red", "orange", "yellow", "beige"))
	pal.hot <- color.palette(steps.hot, c(1,2,2), space="rgb")
	mapcols = c(pal.cold(nlc),pal.hot(nlh))
	return(mapcols)
}
