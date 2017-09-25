as.real <- function(x){
## aliasing: this function was deprcated in R 3.0.1
	return(as.double(x))
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
get.tac.ns = function(x_in){
##
## Used when dealing with TACs (Time Activity Curves).
## Get:
## - number of time points ntp contained in x_in
## - number of voxels n contained in x_in
## Returns:
##   c(n,ntp)
## Note:
## dim(x_in) = [n x ntp, 5]
## Columns of x_in: uptake, weight, x, y, z (Amide format)
##
  # Pick one voxel (arbitrary)
  i=1
  voxdata <- as.matrix(subset(x_in,x_in[,3]==x_in[i,3]&x_in[,4]==x_in[i,4]&x_in[,5]==x_in[i,5]))
  # count number of time points
  ntp=nrow(voxdata)
  # count number of voxels in ROI
  n=nrow(x_in)/ntp
  return(c(n,ntp))
}
getns = function(x_in){
## Alias function for previous code
  get.tac.ns(x_in)
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
is.roi.tac.sorted = function(x_in){
#
# Tests whether input 5-col ROI table x_in (Amide format), with multiple time frames,
# has voxels ordered the same way from all time frames.
# Returns a boolean.
#
  nntp=getns(x_in)
  n=nntp[1]; ntp=nntp[2]
  if(ntp<2){warning("Only one time frame in is.roi.sorted"); return(F)}
  i=1    # Pick one voxel (arbitrary)
  tx=x_in[i,3]
  ty=x_in[i,4]
  tz=x_in[i,5]
  indx=which((x_in[,3]==tx)&(x_in[,4]==ty)&(x_in[,5]==tz))
  tst=indx[2:ntp]-indx[1:(ntp-1)]
  test = prod(tst[2:(ntp-1)]==tst[1:(ntp-2)])
  if(!test){print("!! ROI is not ordered !!")}
  return( as.logical(test) )
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
sort.roi <- function(xin,c1=5,c2=4,c3=3){
##
## Sorts input 5-column ROI table x_in according to 3 columns,
## sorting by column c1 first, then c2, then c3.
## Default: c1=5,c2=4,c3=3.
## Returns re-ordered ROI.
## Specifying values for c1,c2,c3 allows this function to be applied to both Amide and Alice ROI's.
##
  yy=NULL
  for(k in sort(unique(xin[,c1]))){
    xk=subset(xin,xin[,c1]==k)
    if(length(c(xk))>5){xk=xk[order(xk[,c2]),]}
    for(j in unique(xk[,c2])){
      xjk=subset(xk,xk[,c2]==j)
      if(length(c(xk))>5){xjk=xjk[order(xjk[,c3]),]}
      if(length(xjk)){yy=rbind(yy,unlist(xjk))}
    }
  }
  return(yy)
}
sortROIxyz <- function(x_in,ord1=3,ord2=0,ord3=0){
## Alias function for previous code. Should become obsolete ASAP.
  sort.roi(x_in,c1=ord1,c2=ord2,c3=ord3)
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
unique.roi <- function(x,rm.dups=FALSE){
##
## If rm.dups, removes multiple entries for same voxel. Otherwise, affects the mean uptake 
## value to that voxel.
## SHOULD BE RUN AFTER SORTING ROI w.r.t. Z.
## which sorts ROI but also according to weight and as last instance, z...
##
  if(rm.dups){
  	return(x[as.logical(1-duplicated(x[,c(3:5)])),])
  } else {
  	if(sum(duplicated(x[,c(3:5)]))){
  		# then there are duplicates
	  	y=x[as.logical(1-duplicated(x[,c(3:5)])),]
	  	for(i in 1:nrow(y)){
	  		inds=find.v.element(y[i,],x)
	  		y[i,1]=mean(x[inds,1])
	  	}
	  	return(y)
  	} else {
  		return(x)
  	}
  }
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
fast.unique.roi <- function(xx){
#
# Input: binned and sorted ROI xx (Amide format).
#
  	if(sum(duplicated(xx[,c(3:5)]))){
		zz=NULL; 
		gz=unique(xx[,5])
		for(k in gz){
			sz=subset(xx,xx[,5]==k)
			gy=unique(sz[,4])
			for(j in gy){
				sy=subset(sz,sz[,4]==j)
				if(length(sy)>5){
					gx=unique(sy[,3])
					for(i in gx){
						sx=subset(sy,sy[,3]==i)
						if(length(sx)>5){
							zz=rbind(zz,c(mean(sx[,1]),mean(sx[,2]),sx[1,3:5]))
						} else {
							zz=rbind(zz,sx)
						}
					}
				}
			}
		}
	} else {
		zz=xx
	}
	return(zz)
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
find.v.element <- function(vox,roi){
##
## Returns the collection of indices at which vox was found within the input ROI.
## A value of 0 indicates vox was not found within the ROI.
## vox is a 5-vector (uptake, weight, x1, x2, x3), same format as roi.
## The function does not care (test) whether uptakes match.
##
  ind=0
  if(prod(is.element(vox[3:5],roi[,3:5]))){
	for(i in 1:nrow(roi)){
		found=((roi[i,3]==vox[3]) & (roi[i,4]==vox[4]) & (roi[i,5]==vox[5]))
		if(found){ind=c(ind,i)}
	}
  }
  if(sum(ind)>0){ind=ind[2:length(ind)]}
  return(ind)
}
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
is.v.element <- function(vox,roi){
##
## Returns a logical indicating whether vox is in the input ROI.
## Uses find.v.element(vox,roi).
##
  return(as.logical(prod(find.v.element(vox,roi)>0)))
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
round.roi<-function(x,roundto=1,roundzto=roundto){
# Rounds up the 3D coordinate values to 'roundto' digits.
	y=x
	y[,3]=round(x[,3],roundto)
	y[,4]=round(x[,4],roundto)
	y[,5]=round(x[,5],roundzto)	
	return(y)
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
prep.roi <- function(x,doround=FALSE,roundto=1,roundzto=roundto){
##
## Preps Amide ROI x, i.e. rounds, sorts and makes unique.
##
	xo=x
	if(doround){xo=round.roi(x,roundto,roundzto)}
	xo=sort.roi(xo)
	return(unique.roi(xo))
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
is.grid.regular<-function(g){
# Is input grid a regular one?	
	if(length(unique(g))>1){
		return((sum(g[-1]-g[-length(g)])/(g[2]-g[1]))==(length(g)-1))
	} else { 
		return(length(unique(g))==1)
	}
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
is.roi.regular<-function(roi){
# Checks whether all roi grids are regular 	
	su<-function(v){return(sort(unique(v)))}
	gxx=su(roi[,3]); gxy=su(roi[,4]); gxz=su(roi[,5])
	return((is.grid.regular(gxx) & is.grid.regular(gxy) & is.grid.regular(gxz)))
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
bin.roi<-function(roi){
# Bins roi onto .5-thin grids for x and y, and a unit grid fo z.
	bin.v<-function(v){
	# Bins vector v onto a .5-thin grid
		s=sign(v)
		w=floor(abs(v))
		p=(abs(round(v,1))-w)*10
		return(s*(w+.5*(p>=3)+.5*(p>=8)))
	}
	x=roi
	x[,3]=bin.v(x[,3])
	x[,4]=bin.v(x[,4])
	x[,5]=round(x[,5])
	return(x)
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
pair.rois<-function(x,y,rm.na=TRUE){
# Input ROIs x and y must be prep'd first, ie sorted and uniqued.
	z=x[,c(3,4,5,1,2)]
	for(i in 1:nrow(x)){
		indi=find.v.element(x[i,],y)
		if(indi){
			z[i,5]=y[indi,1]
		} else {
			z[i,5]=NA
		}
	}
	if(rm.na){z=z[(!is.na(z[,5])),]}
	return(z)
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
safe.pair.rois<-function(x,y,rto=1,rzto=rto,rm.na=TRUE){
# Same as pair.rois, but tries to make sure that the x and y grids "match"...
	su<-function(v){return(sort(unique(v)))}
	gxx=su(x[,3]); gxy=su(x[,4]); gxz=su(x[,5])
	gyx=su(y[,3]); gyy=su(y[,4]); gyz=su(y[,5])	
	# make sure that (gxz,gxy,gxz) and (gyz,gyy,gyz) match
	if(!(sum(is.element(gxx,gyx))&sum(is.element(gxy,gyy))&sum(is.element(gxz,gyz)))){
		x=bin.roi(x); y=bin.roi(y)
		gxx=su(x[,3]); gxy=su(x[,4]); gxz=su(x[,5])
		gyx=su(y[,3]); gyy=su(y[,4]); gyz=su(y[,5])			
	}
	if(!sum(is.element(gxx,gyx))){
		if(length(gyx)>1){dx=(gyx[2]-gyx[1])/2}else{dx=.5}
		if(sum(is.element((gxx+dx),gyx))>sum(is.element((gxx-dx),gyx))){
			x[,3]=x[,3]+(dx)
		} else {
			x[,3]=x[,3]-(dx)
		}
	}		
	if(!sum(is.element(gxy,gyy))){
		if(length(gyy)>1){dy=(gyy[2]-gyy[1])/2}else{dy=.5}
		if(sum(is.element((gxy+dy),gyy))>sum(is.element((gxy-dy),gyy))){
			x[,4]=x[,4]+(dy)
		} else {
			x[,4]=x[,4]-(dy)
		}
	}		
	if(!sum(is.element(gxz,gyz))){
		if(length(gyz)>1){dz=(gyz[2]-gyz[1])/2}else{dz=1}
		if(sum(is.element((gxz+dz),gyz))>sum(is.element((gxz-dz),gyz))){
			x[,5]=x[,5]+(dz)
		} else {
			x[,5]=x[,5]-(dz)
		}
	}		
	return(pair.rois(x,y,rm.na))	
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
diff.roi <- function(x,y){
##
## Removes y from x, 2 Amide ROIs. Assumes ROIs have been prep'd.
## Example: test diff.roi
## xs=sort.roi(xx); ys=sort.roi(yy)
## xsu=unique.roi(xs); yxs=unique.roi(ys)
## xs[,5]=xs[,5]+1
## inds=(sample(c(1:nrow(xsu)),15))
## x=xsu[inds,]; y=x#yy[inds,]
## cbind(x,y)
## diff(abs(x),abs(y))
##
  inds=duplicated(rbind(x[,c(3:5)],y[,c(3:5)]),fromLast=TRUE)
  return(x[as.logical(1-inds[1:nrow(x)]),])
}

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
image.roi <- function(x_in,slicenb=0,framenb=1,defcol=0,colpal=gray((0:255)/256),main=NULL,...){ 
##
## Plots input ROI x_in (Amide format, read from txt file).
## Input ROI may be box-shaped or not.
## - framenb indicates time frame number to be plotted, =1 by default.
## - slicenb=0 (default) indicates to plot all slices of ROI.
##   If only one specific slice of interest, assign target slice nb to slicenb.
## - defcol = default colour for filled-in voxels (default 0=K)
## - colpal = colour palette
## Output ROI xout should be in Amide format, i.e. s.t.
##      xout[,1] = uptake value
##      xout[,2] = weight
##      xout[,3] = column number
##      xout[,4] = row number
##      xout[,5] = slice number
##  To test:
##    - open data/testplot.xif in Amide
##    - compare with call image.roi(read.table("data/testplotbox.tsv"),1)
##    - or with call image.roi(read.table("data/testplotcyl.tsv"),1)
##    - or with call image.roi(read.table("data/testplotell.tsv"),1)
##    - or with call image.roi(read.table("data/testplotbox.tsv"),22,colpal=topo.colors(256))
##
  nntp=getns(x_in)               ## check if TAC...
  if(nntp[2]==1){ framenb=1 }
  x = x_in[((framenb-1)*nntp[1]+(1:nntp[1])),]
  # prepare gridx
  ggridx = sort(unique(x[,3]))
  gridx = ggridx/(ggridx[2]-ggridx[1])
  gridx = gridx-gridx[1]
  # prepare gridy
  ggridy = sort(unique(x[,4]))
  gridy = ggridy/(ggridy[2]-ggridy[1])
  gridy = gridy-gridy[1]
  # prepare gridz
  ggridz = sort(unique(x[,5]))
  gridz = ggridz/(ggridz[2]-ggridz[1])
  gridz = gridz-gridz[1]
  # nn's
  nx1 = length(gridx)
  nx2 = length(gridy)
  nk  = length(gridz)
  # prepare plottable ROI and plot
  xout=NULL
  slices=as.real(names(table(x[,5])))
  cols=as.real(names(table(x[,3])))
  rows=as.real(names(table(x[,4])))
  if(!prod(slicenb)){
    vslices=slices
  } else {
    vslices=slices[slicenb]
  }
  for(k in vslices) {
    xpz = x[(x[,5]==k),]
    # test ROI shape
    if((as.numeric(nx1)*as.numeric(nx2)*as.numeric(nk))==nrow(x)){ ## then it is a box
      xpo = matrix(sort.roi(xpz,3,4,5)[,1],ncol=nx2,byrow=T)
    } else {                   ## fill the box around the ellipse-ROI with 0's
      # If xpz is not a box-roi, deal with it
      xpo = matrix(0,nx1,nx2)
      for(i in 1:nrow(xpz)){
        xpo[(cols==xpz[i,3]),(rows==xpz[i,4])]=xpz[i,1]
      }
    }
    if(!defcol){    
      if(!prod(slicenb)){
      	  if(!length(main)){main=as.character(which(vslices==k))}
          image(ggridx,ggridy,xpo,zlim=range(xpz[,1]),
            axes="T",main=main,useRaster=TRUE,col=colpal)
      } else {
       	  if(!length(main)){main=as.character(slicenb[which(vslices==k)])}
          image(ggridx,ggridy,xpo,zlim=range(xpz[,1]),axes="T",col=colpal,main=main,useRaster=TRUE,...)
      }
    } else {
      if(!length(main)){main=as.character(which(vslices==k))}
      image(xpo,zlim=range(xpz[,1]),axes="T",
          main=main,
          col=topo.colors(256))  #col=terrain.colors(100))
    }
    xout=rbind(xout,xpo)
  }
  return(xout)	
}
