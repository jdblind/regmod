# FUNCTIONS FOR SARCOMA SYNTHESIS

##################################################
# Polar Transform of x-data given a 2-d mean and 
# radial scaling function sg[nh,2]
##################################################

polar <- function(x,m,sg,thcorr=2) { 
	if(length(x)<=2){x = matrix(x,nc=2)}
	xn = cbind(x[,1]-m[1],x[,2]-m[2])
	th = atan2(xn[,2],xn[,1])
	th = ifelse(xn[,2]<0,th+thcorr*pi,th)
	r = sqrt(xn[,1]^2+xn[,2]^2) 
	u = r/(approx(sg,xo=th,rule=2)$y)
	out = matrix(NA,nr=length(u),nc=4)
	out[,c(1,4)]=th; out[,2]=u; out[,3]=r
	return(out)
} 

polarE <- function(x,m,s,thcorr=2) { 
	out = eigen(s) 
	gam = out$vectors
	lam = out$values 
	lam = ifelse(lam<0,0,lam)
	xi  = atan2(gam[2,1],gam[1,1]) 
	xi = ifelse(gam[2,1]<0,xi+2*pi,xi) 
	xii = 0-xi
	ee = sqrt( 1- sqrt(lam[1]*lam[2])/( .5*(sqrt(lam[1])+sqrt(lam[2])))^2 )
	xn = cbind(x[,1]-m[1],x[,2]-m[2]) 
	th = atan2(xn[,2],xn[,1])
	th = ifelse(xn[,2]<0,th+thcorr*pi,th)
	r = sqrt(xn[,1]^2+xn[,2]^2)
	# E Radius
	re = (1-ee)*(1+ee)/sqrt((1-ee)^2*cos(th+xii)^2+(1+ee)^2*sin(th+xii)^2)
	u = r/re
	out = matrix(NA,nr=length(u),nc=3)
	out[,1]=th; out[,2]=u; out[,3]=re
	return(out)
}

##################################################
# Unimodal Smoother
##################################################

unismooth <- function(x,y,w) {
	n = length(x)
	n0 = ifelse(n>200,n,max(n,3))
	nk = max(n0,length(x)/5) 
	nk = min(100,nk)
	out = smooth.spline(x,y,w,nknots=nk) 
	zhat = out$y
	xv = out$x 
	ns = length(xv)
	zmax = max(zhat) 
	xmax = xv[order(-zhat)[1]] 
	zu = zhat*0 
	zm = zu 
	if(((xmax>xv[1])&(xmax<xv[ns]))) { 
		zu[xv<=xmax]=isoreg(zhat[xv<=xmax])$yf 
		zu[xv>xmax]=0-isoreg(0-zhat[xv>xmax])$yf 	
		zm[xv<=xmax]=max(zu) 
		zm[xv>xmax]=zu[xv>xmax]
	}
	if((xmax<=xv[1])){ 
		zu[xv>=xmax] = 0-isoreg(-zhat[xv>=xmax])$yf
		zm=zu				     
	}
	if((xmax>=xv[ns])){
		zu[xv<=xmax]=isoreg(zhat[xv<=xmax])$yf   
		zm=median(zu)+0*zu
	}
	out = matrix(NA,nr=length(xv),nc=2)
	out[,1]=xv; out[,2]=zu
	return(out)
}

#######################################################
# Data Threshold Values   
#######################################################

thresholds <- function(x,y,pcut) {
	par(mfrow=c(2,2))
	# SELECT REGION OF DATA
	# CUT-OFF LOWER AND UPPER EXTREMES (LOWEST VALUE IS AN ARTIFACT OF PAS)
	
	#Create value (2-Pass)
	my=min(y) 
	yl=y
	yl=yl[y>my] 
	ny=length(yl)
	yll=sort(yl)[max(.01*ny,1)] 
	ylu=sort(yl)[min(.99*ny,ny)] 
	yv=yl[(yl>yll)&(yl<ylu)]
	ylow=yll+.5*mad(yv) 
	yhigh=ylu-.5*mad(yv)
	yv=yl[(yl>ylow)&(yl<yhigh)]
	md=mad(yv)
	my=my+.001*md 
	
	yl=y
	yl=yl[y>my] 
	ny=length(yl)
	yll=sort(yl)[max(.01*ny,1)] 
	ylu=sort(yl)[min(.99*ny,ny)] 
	yv=yl[(yl>yll)&(yl<ylu)]
	ylow=yll+.5*mad(yv) 
	yhigh=ylu-.5*mad(yv)
	
	# Data Selection
	yv=yl[(yl>ylow)&(yl<yhigh)] 
	yy=y[y>my]
	
	ny=length(yv)
	y1=sort(yv)[1:(.25*ny)] 
	pp=(c(0:(ny-1))+.5)/(ny)
	plot(sort(yv),qnorm(pp),pch=".")
	value=sort(yv)[pcut*ny]	# CUT OUT pcut% of  VOXELS
	thres=value+.2*mad(yv) 
	abline(v=thres,col=c(5))
	abline(v=value,col=c(6))
	par(new=T)
	hist(yv,xlab=" ",ylab="",main="",axes=F)
	thresn=ylow+.1*mad(yv)
	
	pcent=100*length(yy[yy>thres])/length(yy)
	pcentv=100*length(yv[yv>thres])/length(yv)
	
	hist(yy)
	abline(v=thres,col=c(5))
	abline(v=value,col=c(6))
	return(c(thres,pcut,round(c(100-pcentv,100-pcent,length(y),length(yy),length(yv)))) )
}

#########################################################
# AMIDE ROIs and ANALYSIS
#########################################################

pasbins <- function(z,nz,nv,n){
	zb=sort(z)[c(1,c(1:(nz-1))*nv,n)] 
	if(sum(duplicated(zb))){
		# must preserve length nz, so unique(zb) is forbidden...
		zb = seq(min(z),max(z),len=nz+1)
	}
	zb[1]=zb[1]-max(.1e-9,.001*(zb[2]-zb[1])) 	# pas bins
	return(zb)
}	

pas <- function(x,y,thres){
# Principl Axis Sampling (projection), feeds into arotx().
# Used ot be top part of arotx(). Its output is useful in 
# other places too.
	#Create value (2-Pass)
	my=min(y) 
	yl=y
	yl=yl[y>my] 
	ny=length(yl)
	yll=sort(yl)[max(.01*ny,1)] 
	ylu=sort(yl)[min(.99*ny,ny)] 
	yv=yl[(yl>yll)&(yl<ylu)]
	ylow=yll+.5*mad(yv) 
	yhigh=ylu-.5*mad(yv)
	yv=yl[(yl>ylow)&(yl<yhigh)]
	md=mad(yv) 
	my=my+.001*md 
	#
	yl=y
	yl=yl[y>my]
	ny=length(yl) 
	yll=sort(yl)[max(.01*ny,1)] 
	ylu=sort(yl)[min(.99*ny,ny)] 
	yv=yl[(yl>yll)&(yl<ylu)]
	ylow=yll+.5*mad(yv) 
	yhigh=ylu-.5*mad(yv)
	#
	mi=min(y[((y>ylow)&(y<yhigh))]) 
	ma=max(y[((y>ylow)&(y<yhigh))])
	# Data Selection for Mean and Variance Calculations
	yv=yl[(yl>ylow)&(yl<yhigh)]
	xs=x[((y>ylow)&(y<yhigh)),]  
	zs=ifelse(y[((y>ylow)&(y<yhigh))]>thres,(y[((y>ylow)&(y<yhigh))]-mi)/(ma-mi),.1e-9) 
	ns=length(zs)
	mm = c(weighted.mean(xs[,1],zs),weighted.mean(xs[,2],zs),weighted.mean(xs[,3],zs))
	# Mean and Covariance
	xs=cbind(xs[,1]-mm[1],xs[,2]-mm[2],xs[,3]-mm[3]) 
	s=matrix(rep(0,9),ncol=3)
	s[1,1]=sum(xs[,1]*xs[,1]*zs)  
	s[1,2]=sum(xs[,1]*xs[,2]*zs)
	s[2,1]=s[1,2]
	s[1,3]=sum(xs[,1]*xs[,3]*zs)
	s[3,1]=s[1,3]
	s[2,2]=sum(xs[,2]*xs[,2]*zs)
	s[2,3]=sum(xs[,2]*xs[,3]*zs)
	s[3,2]=s[2,3]
	s[3,3]=sum(xs[,3]*xs[,3]*zs)
	s=s/sum(zs) 
	s=s+diag(rep(1,3))*0.01*mean(diag(s)) 
	gam=eigen(s)$vectors[,c(3,2,1)]
	for(j in 1:3) { 
		u=gam[j,]/max(abs(gam[j,]))
		if(min(u)<=-.99) { gam[j,]=0-gam[j,] }
	}
	# rows for amide rotation
	mu = c(mm[1],mm[2],mm[3])
	return(list(xs=xs,ys=y[((y>ylow)&(y<yhigh))],zs=zs,
				mu=mu,gam=gam,S=s,s=sqrt(eigen(s)$values)[c(3,2,1)]))
}

arotx <- function(x,y,thres) {
# Merely reshapes output of pas() for subsequent use...
# CUT-OFF LOWER AND UPPER EXTREMES (LOWEST VALUE IS AN ARTIFACT OF PAS)
	pas.out = pas(x,y,thres)
	S = pas.out$S
	s = pas.out$s
	mu = pas.out$mu
	gam = pas.out$gam
	amide = round(cbind(mu,t(gam),sqrt(rev(eigen(S)$values))),3)
	output = list(amide=amide,c=mu,s=s,Gxi=gam,S=pas.out$S)
	return(output)				# Manual entries for AMIDE rotation
}


#######################################################
# Given (PAS) transformed Data x[,3] and using x[,3] as 
# the axial direction compute z,m,vol,sg for the tumor 
# boundary
#######################################################

msg<-function(xt,y,thres,alpha,dd,nb,nres=100,doplot=TRUE){
	out=ccore(xt,y,thres,dd,nb,doplot)
	z=out[,1] 
	mk=out[,3:4]
	sk=out[,5:8]
	outn=gsg(xt,y,z,mk,sk,thres,alpha,nres,doplot)
	e <- new.env()
	e$ccore=out
	e$nb=nb
	e$z=outn$z 
	e$vol=outn$vol 
	e$m=outn$m 
	e$th=outn$th 
	e$sg=outn$sg
	# add-ons:
	e$z.pasbins=z 
	e$z.midpoints=outn$z
	e$gsg.ee=outn$ee
	e$gsg.sg0=outn$sg0	
	e$S=outn$s
	e$lams=outn$lams		
	e$siv=outn$siv
	e$res=outn$res
	e$u=outn$u
	e$re=outn$re
	e$uu=outn$uu # actually = v[siv]
	return(e)
}

bpts <-function(z,m,th,sg) {
# Creates Boundary Points
	nz=length(z)
	nt=length(th)
	xxvs=NULL
	for(j in 1:nz) { # PAS contours
		r=sg[j,] 
		xn1=r*cos(th) 
		xn2=r*sin(th)
		xn= t(cbind(xn1,xn2))
		xvs = xn + m[j,]
		xxvs=cbind(xxvs,c(t(xvs)))
	}
	bxt=NULL 
	for(j in 1:nz) { 
		bxt=rbind(bxt,cbind(matrix(xxvs[,j],ncol=2),rep(z[j],nt))) 
	}
	return(bxt)
}

npts <-function(z,m,th,sg) {
# Normals to Boundary Points
	nz=length(z)
	nt=length(th)
	# a is the derivative (differences) of surface pts wrt th (nt index)
	xxvs=NULL
	for(j in 1:nz) { # PAS contours
		r=sg[j,] 
		dr=r*0 
		dr[1]=(r[2]-r[nt]) 
		dr[2:(nt-1)]=r[3:nt]-r[1:(nt-2)]
		dr[nt]=r[1]-r[nt-1]
		xn1=0-r*sin(th)+dr*cos(th) 
		xn2=r*cos(th)+dr*sin(th)
		xn= t(cbind(xn1,xn2))
		xvs = xn 
		xxvs=cbind(xxvs,c(t(xvs)))
	}
	a=NULL 
	for(j in 1:nz) { 
		a=rbind(a,cbind(matrix(xxvs[,j],ncol=2),rep(0,nt))) 
	}
	# b is the derivative (differences) of surface pts wrt z (nz index)
	dzsg=sg*0
	for(i in 1:nt) { 
		dzsg[,i]=c(0,sg[c(3:nz),i]-sg[c(1:(nz-2)),i],0)
	}
	dzm=m*0 
	dzm[,1]=c(0,m[c(3:nz),1]-m[c(1:(nz-2)),1],0)
	dzm[,2]=c(0,m[c(3:nz),2]-m[c(1:(nz-2)),2],0)
	dz = c(z[2]-z[1],z[c(3:nz)]-z[c(1:(nz-2))],z[nz]-z[nz-1])
	xxvs=NULL
	for(j in 1:nz) { # PAS contours
		dr=dzsg[j,] 
		xn1=dr*cos(th) 
		xn2=dr*sin(th)
		xn= t(cbind(xn1,xn2))
		xvs = xn+dzm[j,]
		xxvs=cbind(xxvs,c(t(xvs)))
	}
	b=NULL 
	for(j in 1:nz) { 
		b=rbind(b,cbind(matrix(xxvs[,j],ncol=2),rep(dz[j],nt))) 
	}
	nxt=NULL 
	nn=nrow(a) 
	for(i in 1:nn) { 
		nxt=rbind(nxt,axb(c(a[i,]),c(b[i,]))) 
	}
	return(nxt)	
}

vis <- function(xv,nxt) {
# Visibility of Points
	nxv=xv/sqrt(sum(xv*xv)+.1e-9)
	isee=NULL 
	n=nrow(nxt)
	for(i in 1:n) { 
		yv=c(nxt[i,]) 
		nyv=yv/sqrt(sum(yv*yv)+.1e-9)
		isee=c(isee,ifelse(sum(nxv*nyv)>0,1,0))  # angle less than 90 dgrees
	}
	return(isee)
}

#######################################################
# APPROXIMATE THE CENTRAL CORE AND ELLIPTICAL SHAPE OF 
# THE TUMOR USING THE PAS AXIS
#######################################################

ccore <- function(x,y,thres,dd,nb,doplot=T) { 
	# Average of 100 voxels per PAS slice (dd is PAS resolution)
	n=length(y) 
	nv=max(100,n/nb) 
	nb=round(max(1,n/nv)) 
	z=pasbins(x[,3],nb,nv,n)
		
	# SELECT REGION OF DATA
	# CUT-OFF LOWER AND UPPER EXTREMES (LOWEST VALUE IS AN ARTIFACT OF PAS)
	# Create value (2-Pass)
	my=min(y) 
	yl=y
	yl=yl[y>my] 
	ny=length(yl) 
	yll=sort(yl)[max(.01*ny,1)] 
	ylu=sort(yl)[min(.99*ny,ny)] 
	yv=yl[(yl>yll)&(yl<ylu)]
	ylow=yll+.5*mad(yv) 
	yhigh=ylu-.5*mad(yv)
	ylow=yll
	yhigh=ylu
	yv=yl[(yl>ylow)&(yl<yhigh)]
	md=mad(yv) 
	my=my+.001*md 
	
	# Data Selection for Mean and Variance Calculations
	yv=yl[(yl>ylow)&(yl<yhigh)] 
	yhigh=max(y)
	yy=y  
	yo=min(y[((y>ylow)&(y<yhigh)&(y>thres))]) 
	mo=max(y[((y>ylow)&(y<yhigh)&(y>thres))])
	
	mk=NULL; lamk=NULL; sk=NULL; kv=NULL; wv=NULL
	xik=NULL; ek=NULL; cork=NULL; zk=NULL; wk=NULL
	
	y2=thres
	y3=thres+.125*(yhigh-thres)
	y4=thres+.25*(yhigh-thres)
	
	for(k in 1:nb) { 
		yk = yy[(x[,3]<=z[k+1])&(x[,3]>z[k])] 
		no = length(yk) 
		tho = sort(yk)[max(no/3,no-600)]
		thm = sort(yk)[no]
		xs=x[(x[,3]<=z[k+1])&(x[,3]>z[k]),c(1:2)]
		xs=xs[((yk>ylow)&(yk<yhigh)),c(1:2)]
		zs=ifelse(yk[((yk>ylow)&(yk<yhigh))]>tho,((yk[((yk>ylow)&(yk<yhigh))]-tho+.1e-5)/(thm-tho+.1e-5)) ,0)
		yk=yk[((yk>ylow)&(yk<yhigh))]
		  
		mm = c(weighted.mean(xs[,1],zs),weighted.mean(xs[,2],zs)) ### not used

		# Center and Covariance
		out=center(xs,yk,dd)
		m=out$xc
		s=out$s
		mk=rbind(mk,m)
		
		# Shape Analysis
		s=s+diag(rep(1,2))*0.01*mean(diag(s)) 
		sk=rbind(sk,c(s)) 
		corr= s[1,2]/sqrt(s[1,1]*s[2,2]) 
		cork=c(cork,corr)
		out=eigen(s) 
		gam=out$vectors
		gam[,1]=gam[,1]/sign(gam[,1])[order(abs(gam[,1]))[2]]
		gam[,2]=gam[,2]/sign(gam[,2])[order(abs(gam[,2]))[2]]
		lam=out$values
		lam=ifelse(lam>0,lam,0)
		lamk=rbind(lamk,lam) 
		xi=atan2(gam[2,1],gam[1,1]) 
		xi=ifelse(gam[2,1]<0,xi+2*pi,xi)
		xik=c(xik,xi)
		ee=sqrt( 1 - sqrt(lam[1]*lam[2])/(.5*(sqrt(lam[1])+sqrt(lam[2])))^2 )
		ek=c(ek,ee)
		zk=c(zk,.5*(z[k]+z[k+1]))
		wk=c(wk,sum(zs))
		# NICE PLOT
		if(doplot){
			if(k<(-1000)){
		 		plot(xs[,1],xs[,2],pch="*",xlim=range(c(x[,1:2])),
		 				ylim=range(c(x[,1:2])),
						xlab=" ",ylab=" ",axes=F,col=c(8))
				points(xs[((yk>y2)&(yk<y3)),1],xs[((yk>y2)&(yk<y3)),2],pch="*",col=c(4))
				points(xs[(yk>y3)&(yk<y4),1],xs[(yk>y3)&(yk<y4),2],pch="*",col=c(2))
				points(xs[yk>y4,1],xs[yk>y4,2],pch="*",col=c(7)) 		
				abline(h=m[2]);abline(v=m[1])
				par(new=TRUE)
				oo=polarE(xs,c(m[1],m[2]),s)
				plot(oo[,2],yk,pch=".",axes="F",xlab=" ",
						ylab=" ",main=as.character(round(zk[k],2)))
				axis(2)
			}
		}
	}
	leftzk = min(zk[cumsum(wk)>max(wk)]) 
	rightzk=max(rev(zk)[cumsum(rev(wk))>max(wk)])
	dfs=max(2,.2*length(zk))
	m1=smooth.spline(zk,smooth(mk[,1]),wk,df=dfs)$y ; m10=mk[,1] ; mk[,1]=m1 
	m2=smooth.spline(zk,smooth(mk[,2]),wk,df=dfs)$y ; m20=mk[,2] ; mk[,2]=m2
	s1=smooth.spline(zk,smooth(sqrt(sk[,1])),wk,df=dfs)$y
	s2=smooth.spline(zk,smooth(sqrt(sk[,4])),wk,df=dfs)$y
	rho=sin(smooth.spline(zk,smooth(asin(cork)),wk,df=dfs)$y)
	sk[,1]=s1*s1 
	sk[,2]=s1*s2*rho 
	sk[,3]=sk[,2]
	sk[,4]=s2*s2
	
	w=z*0 
	w=approx(zk,wk,xo=z,rule=2)$y
	out=cbind(z,w,
			approx(zk,mk[,1],xo=z,rule=2)$y,approx(zk,mk[,2],xo=z,rule=2)$y,
			approx(zk,sk[,1],xo=z,rule=2)$y,approx(zk,sk[,2],xo=z,rule=2)$y,
			approx(zk,sk[,3],xo=z,rule=2)$y,approx(zk,sk[,4],xo=z,rule=2)$y,
			approx(zk,cork,xo=z,rule=2)$y,approx(zk,ek,xo=z,rule=2)$y)
	return(out)
} # end ccore

####################################################################################################
####################################################################################################

center <- function(xv,yv,dd) { 
#
# FIND LOCAL CENTER OF THE DATA
#
	# Peel 40% of Layers off the REGION
	xn=xv 
	kp=0 
	while(nrow(xn)>(5+length(chull(xn)))) { 
		xn=xn[-chull(xn),] 
		kp=kp+1 
	}
	xn=xv 
	yn=yv; 
	p=NULL 
	kv=0
	while(kv<round(max(.4*kp,min(kp/2,5)))) { 
		p=c(p,chull(xn)) 
		pv=chull(xn) 
		xn=xn[-pv,] 
		yn=yn[-pv];kv=kv+1
	}
	xouter = xn[chull(xn),] 
	# outer layer of points - center location must be inside this

	# Guess from Quadratic Fit to intensity
	wn=yn*0+1/length(yn)
	x1=xn[,1] ; x2=xn[,2] ; x11=x1*x1 ; x22=x2*x2 ; x12=x1*x2 ; 
	bmat=matrix(rep(0,4),ncol=2)
	b=rep(0,2)
	out=lm(yn~x1+x2+x11+x12+x22,weights=wn)
	ynh=out$fitted.values 
	coef=summary(out)$coef[,1]
	b[1]=0-coef[2] 
	b[2]=0-coef[3] 
	bmat[1,1]=2*coef[4]
	bmat[1,2]=coef[5]
	bmat[2,1]=bmat[1,2]
	bmat[2,2]=2*coef[6]
	e=eigen(bmat) 
	l=e$val;gam=e$vec

	if((l[1]*l[2]<0)) { 
		if(abs(l[1])<=abs(l[2])) { l[1]=.01*l[2] }
        if(abs(l[2])<=abs(l[1])) { l[2]=.01*l[1] } 
    }
	
	if(l[1]>0&l[2]>0) { xc0=xn[order(ynh)[1],]    }		#Valley
	if(l[1]<0&l[2]<0) { xc0=xn[order(-ynh)[1],]   }		#Hill
	s = gam[,1]%*%t(gam[,1])/abs(l[1])+gam[,2]%*%t(gam[,2])/abs(l[2])
	ss=s[1,1]+s[2,2]
	ss0=var(xv[,1])+var(xv[,2]^2)
	s=s*ss0/ss

	# DECIDE BETWEEN 4 POINTS
	bflag=1
	xcm=solve(bmat,b);nn=nrow(xn) ; if(min(chull(rbind(xcm,xn)))<2) { xcm=xc0; bflag=2 }	# Model Point
	xcyL1=xn[order(-yn)[1],]					# y-max 
	xn=xv ;yn=yv; p=NULL ; kv=0
	while(kv<round(max(.8*kp,min(kp,3)))) { p=c(p,chull(xn)) ; pv=chull(xn) ; xn=xn[-pv,] ; yn=yn[-pv];kv=kv+1}
	xcyL=xn[order(-yn)[1],]						# y-max (inner core max)
	xcyS=xn[order(yn)[1],]						# y-low (inner core min)
	ms=rbind(xcm,xcyL1,xcyL,xcyS)
	nxn=nrow(xn);ms=rbind(ms,xn[sample(c(1:nxn),min(70,nxn)),])     # random sample of inner shell points

	s=s+.000001*mean(diag(s))
	ua=.5*log(s[1,1]/s[2,2]) ; ua=max(log(1/3),ua) ; ua=min(log(3),ua) 
	pua=(ua-log(1/3))/(log(3)-log(1/3)); pua=max(.001,pua) ; pua=min(pua,.999)
	ua=log(pua/(1-pua))
	ur= s[1,2]/(.000000001+sqrt(s[1,1]*s[2,2])); ur=max(-.8,ur);ur=min(.8,ur)
	pur=(ur-(0-0.8))/(0.8-(0-0.8)); pur=max(.001,pur) ; pur=min(pur,.999)
	ur=log(pur/(1-pur))

	s=sv(ua,ur)

	rss=rep(1.e+14,nrow(ms)) ; sn=s;sns=matrix(c(s),nrow=1); wyt=ifelse(yv>0,yv,.1e-8) ; wyt=wyt/max(wyt)
	for(j in bflag:nrow(ms)) { 
                        sn[1,1]=sum(wyt*(xv[,1]-ms[j,1])^2)/sum(wyt);sn[2,2]=sum(wyt*(xv[,2]-ms[j,2])^2)/sum(wyt)
			sn[1,2]=sum(wyt*(xv[,2]-ms[j,2])*(xv[,1]-ms[j,1]))/sum(wyt);sn[2,1]=sn[1,2]
			ss0=mean(diag(sn)) ; sn=sn/ss0;sns=rbind(sns,c(sn))
			rss[j]=rssvf2(ms[j,1],ms[j,2],sn,xv,yv)
			}
	xc=c(ms[order(rss)[1],]); rsso=min(rss,na.rm=T); xcn=xc ;sn=matrix(sns[order(rss)[1],],ncol=2)

	# Now do steepest descent starting at xc and but restricted to the convex hull of the inner layers
	ibetter=1	# PASS 1
	for(j in 1:2){
		if(ibetter>0){
			out = descent(xc,sn,xv,yv,dd,xouter)
			ibetter=0
			if(out$rssn<rsso){
				xc=out$xcn
				rsso=out$rssn
				s=out$s
				ibetter=1
			}
		}
	}
	
	e <- new.env()
	e$xc=xc ; e$s=s
	return(e)
}

# ---------------------------------------------------------------
descent <- function(xo,s,x,y,dd,x30){
#
# STEEPEST DESCENT DIRECTION
#
	s=s+.000001*mean(diag(s))*diag(rep(1,2))
	ua=.5*log(s[1,1]/s[2,2]) ; ua=max(log(1/3),ua) ; ua=min(log(3),ua) 
	pua=(ua-log(1/3))/(log(3)-log(1/3)); pua=max(.001,pua) ; pua=min(pua,.999)
	ua=log(pua/(1-pua)) ; dua=.1 
	ur= s[1,2]/(.000000001+sqrt(s[1,1]*s[2,2])); ur=max(-.8,ur);ur=min(.8,ur)
	pur=(ur-(0-0.8))/(0.8-(0-0.8)); pur=max(.001,pur) ; pur=min(pur,.999)
	ur=log(pur/(1-pur)) ; dur=.1 

	s=sv(ua,ur)

	# Calculate descent direction (+/- dd from xo)
	p1 = (rssvf2((xo[1]+dd),xo[2],s,x,y)-rssvf2((xo[1]-dd),xo[2],s,x,y))/(2*dd)
	p2 = (rssvf2(xo[1],(xo[2]+dd),s,x,y)-rssvf2(xo[1],(xo[2]-dd),s,x,y))/(2*dd)
	su=sv(ua+dua,ur) ; sl=sv(ua-dua,ur)
	pu=(rssvf2(xo[1],xo[2],su,x,y)-rssvf2(xo[1],xo[2],sl,x,y))/(2*dua)
	su=sv(ua,ur+dur) ; sl=sv(ua,ur-dur)
	pr=(rssvf2(xo[1],xo[2],su,x,y)-rssvf2(xo[1],xo[2],sl,x,y))/(2*dur)
	p = 0.00001*c(1,1,1,1)-c(p1,p2,pu,pr)/sqrt(p1*p1+p2*p2+pu*pu+pr*pr+.1e-20) # Negative Gradient
	p=p/sqrt(sum(p*p))
	# Scale to achieve a step-size  of dd in components 1 and 2
	fac = 1*dd/sqrt(p[1]*p[1]+p[2]*p[2])

	# Optimize descent step
	rss=NULL ; xs=NULL ; rssb=1.e+30 ; xb=xo ; uab=ua;urb=ur; iflag=0 ; bflag=1
	for(j in 1:10) { xx = xo+(j-1)*fac*p[1:2] 
			uav=ua+(j-1)*fac*p[3] 
			urv=ur+(j-1)*fac*p[4] 
			s=sv(uav,urv)
			if(min(chull(rbind(xx,x30)))>1)  { rssv=rssvf2(xx[1],xx[2],s,x,y) }
			if(min(chull(rbind(xx,x30)))<=1) { rssv=rssvf2(xb[1],xb[2],s,x,y) ; xx=xb}
			if(rssv<=rssb) { xb=xx ; uab=uav ; urb=urv ; rssb=rssv}
			xs=rbind(xs,c(xx)) ; rss=c(rss,rssv)
			}

	e <- new.env()
	e$xcn=c(xs[order(rss)[1],]) ; e$s=sv(uab,urb) ; e$rssn=rssb
e
}

# ---------------------------------------------------------------
sv <- function(ua,ur) {
#
# SHAPE MATRIX GIVEN ua and ur
#
	ua=max(-20,ua);ua=min(20,ua) 
	x = log(1/3) + (exp(ua)/(1+exp(ua)))*(log(3)-log(1/3)) ; a=exp(2*x) 
	s11=a/(1+a) ; s22=1/(1+a); s1=sqrt(s11) ; s2=sqrt(s22)

	ur=max(-20,ur);ur=min(20,ur)
	rl=-.8 ; rh=.8 ; r =  rl+ (exp(ur)/(1+exp(ur)))*(rh-rl) 

	s=matrix(rep(0,4),ncol=2)
	s[1,1]=s11 ; s[1,2]=s1*s2*r ; s[2,1]=s[1,2] ; s[2,2]=s22
	
	s=s+diag(rep(mean(diag(s))*.001,2))
	return(s)
}

# ---------------------------------------------------------------
rssvf2<-function(m1,m2,s,x,y){ 
# Unimodal Model fit based on given center and covariance matrix
	s=s+diag(rep(1,2))*0.01*mean(diag(s))
	oo=polarE(x,c(m1,m2),s); 
	u=oo[,2]; 
	wt=1+u*0
	out=unismooth2(u,y)
	yh=approx(out[,1],out[,2],xo=u,rule=2,ties="ordered")$y
	rssv=sum(wt*(y-yh)^2) 
	return(rssv)
}

# ---------------------------------------------------------------
unismooth2 <- function(u,y) {
#
# Unimodal Smooth
#
	yy=smooth(y[order(u)]) 
	uu=sort(u)
	um=min((min(uu)+max(uu))/2,uu[order(-yy)[1]])	# MODE <= (min(u)+max(u))/2
	u1=uu[uu<=um] 
	my1=mean(yy[uu<=um])
	u2=uu[uu>=um] 
	my2=mean(yy[uu>=um])
	if((length(u1)<=3)|(um<=min(uu))){  
		yh1=u1*0+my1 
	}
	if((um>min(uu))&(length(u1)>3))  { 
		u1=uu[uu<=um]  
		y1=yy[uu<=um][order(u1)]
		u1=sort(u1) 
		sy1=smooth(y1)  
		yh1=isoreg(u1,sy1)$yf   
	}	
	if((length(u2)<=3)|(um>=max(uu))){  
		yh2=u2*0+my2 
	}
	if((um<max(uu))&(length(u2)>3))  { 
		u2=uu[uu>=um] 
		y2=yy[uu>=um][order(u2)]
		u2=sort(u2)
		sy2=0-smooth(y2);
		yh2=0-isoreg(u2,sy2)$yf  
	}
	ux=c(u1,u2)
	yx=c(yh1,yh2)
	nx=length(as.real(names(table(ux)))) 
	yh=approx(ux,yx,xo=u,rule=2,ties="ordered")$y
	if(nx>10){ 
		out=smooth.spline(ux,yx,all.knots=FALSE,nknots=max(2,nx/2))
		yh=approx(out$x,out$y,xo=u,rule=2,ties="ordered")$y 
	}
	return(cbind(sort(u),yh[order(u)]))
}


#######################################################
# sg by slice  
#######################################################

# ---------------------------------------------------------------
bc <- function(th,u,yh,y,thres,nres,dfv,igraph,doplot=TRUE) {
#
# BOUNDARY CONTOUR in u=r/re parameter keep uvol and voxel count
#
	if(min(yh)>thres) {  # constant sg(th) phi will be small
		uo=max(u)
		o = cbind(seq(0,2*pi,length=nres),rep(uo,nres))
		dfvv=1  
	}
	if(max(yh)<thres) { # constant sg(th) phi will be large
		uo=min(u)
		o = cbind(seq(0,2*pi,length=nres),rep(uo,nres)); 
		dfvv=1 
	}
	if(max(yh)>thres&min(yh)<thres) { 
		um=u[order(-yh)[1]]
	    uo=min(u[(u>um)&(yh<(thres))]) 
	    n=length(yh)
	    # select points on the boundary
	    au=abs(u-uo)
	    uvs= sort(au)[1:min(800,max(1,.7*n))]
	    cut=max(uvs)
	    na=length(uvs)
	    dfvv=max(1,min(1.5*na,3*dfv))
		s= sqrt(var(y[au<cut]))
		xth=th[au<=cut] 
		wt = 1/(1+.5*((y[au<=cut]-thres)/s)^2)
		o = sbc(th[au<=cut],u[au<=cut],y[au<=cut],thres,dfvv,nres) 
	}
	uvol=sum(o[,2]^2*pi/nres) 
	volv=length(u[u<=uo])
	if(doplot){
		if((igraph<=-201)&(max(yh)>thres&min(yh)<thres)) {
			plot((u[au<=cut]*cos(xth)),(u[au<=cut]*sin(xth)),
					ylim=c(-max(u[au<=cut]),max(u[au<=cut])),
					xlim=c(-max(u[au<=cut]),max(u[au<=cut])),
					xlab=" ",ylab=" ",pch="*",axes=F,main="HI")
			lines(o[,2]*cos(o[,1]),o[,2]*sin(o[,1]),col=c(4))
		}
	}
	e <- new.env()
	e$volu=uvol
	e$volc=volv
	e$o=o
	return(e)
}

# ---------------------------------------------------------------
sbc <- 	function(th,u,y,thres,dfvv,nres) {
#
# smooth boundary
#
	na=length(th) 
	ub=mean(u) 
	if(na>20){ 
		e=fitth(u,y,thres)
		ub=e$uv
	}
	uth=ub
	o = cbind(seq(0,2*pi,length=nres),rep(uth,nres))
	# divide 0,2i into nb bins with an average of at least 20 points per bin
	nbin=round(na/20)
	if(nbin>5) {
		u=u[order(th)] 
		y=y[order(th)] 
		th=sort(th) 
		thj=NULL ; uthj=NULL ; wthj =NULL ; nvj=NULL
		for(j in 1:nbin) { 
			thl=(j-1)*2*pi/nbin 
			thh=(j)*2*pi/nbin 
	        thj = c(thj ,(thl+thh)/2)
			nv=length(th[th<=thh&th>=thl]) 
			val=ub 
			wt=0
			if(nv>3) { 
				e=fitth(u[th<=thh&th>=thl],y[th<=thh&th>=thl],thres)
				val=e$uv
				wt=e$seuv^2
				wt=1/max(.1e-5,wt) 
			}
			uthj=c(uthj,val) ; wthj=c(wthj,wt) ; nvj=c(nvj,nv)
		}
		xna=cbind(thj,uthj,wthj)
		xna=na.omit(xna)
		thj=xna[,1]
		uthj=xna[,2]
		wthj=xna[,3]
		uvals = smooth(rep(uthj,3))
		wvals=rep(wthj,3)
		out=approx(smooth.spline(c(thj-2*pi,thj,thj+2*pi),uvals,wvals,df=dfvv),
						xout=seq(0,2*pi,length=nres))
		o = cbind(out$x,out$y)  
	}
	return(o)
}

# ---------------------------------------------------------------
dgammp<- function(tau,phi,beta,cv,u,alpha) {
# Add the penalty term to gamma
	yhatp=c(beta*dgamma(tau+phi*u,alpha),sqrt(beta)*cv)
	return(yhatp)
}

# ---------------------------------------------------------------
gsg <- function(x,y,z,mk,sk,thres,alpha,nres=100,doplot=TRUE) {
	n = length(y) 
	mo = max(y) 
	lo = min(y[y>thres]) 
	myv = sqrt(var(sort(y)[c(1:(.99*n))]))
	dfbc = 8				# resolution for sg and the df for bound countour smoother
	nb = length(z)-1
	low = lo+myv*.001 
	pv = 0.1
	beta0s=NULL;betas=NULL;taus=NULL ; phis=NULL
	sebetas=NULL;setaus=NULL;sephis=NULL
	wts=NULL;uu=NULL;yyuu=NULL;yyw=NULL;yyhat=NULL 
	rssa=NULL;sg0=NULL;sgU=NULL;vols=NULL;zvs=NULL;ms=NULL
	res=NULL;lams=ss=ees=NULL
	
	nbb = 5
	y2=thres
	y3=thres+.125*(mo-thres)
	y4=thres+.3*(mo-thres)
	nthres=length(y[y>thres])
	if(nthres>4){
		syt=sort(y[y>thres]) 
		y3=syt[.33*nthres] 
		y4=syt[.66*nthres]
	}
	k33=round(nb*.33) 
	k66=round(nb*.66) 
	indx=c(k33,k66)
	
	for(k in 1:nb) {
		fac=0
		beta0=0
		beta=1
		tau=qgamma(.1,alpha)
		phi=1.
		wwt=0
		setau=sephi=sebeta=1000
		zv=(z[k+1]+z[k])/2
		yu=y[((x[,3]<=z[k+1])&(x[,3]>=z[k]))] 
		xs=x[((x[,3]<=z[k+1])&(x[,3]>=z[k])),c(1:2)]
		m=rep(0,2) 
		s=rep(0,4)
		m[1]=approx(z,mk[,1],xo=zv,rule=2)$y 
		m[2]=approx(z,mk[,2],xo=zv,rule=2)$y 
		s[1]=approx(z,sk[,1],xo=zv,rule=2)$y 
		s[2]=approx(z,sk[,2],xo=zv,rule=2)$y 
		s[3]=s[2]
		s[4]=approx(z,sk[,4],xo=zv,rule=2)$y 
		out = polarE(xs,m,matrix(s,ncol=2)) 
		re = c(out[,3])
		u = c(out[,2])
		th = c(out[,1])
		out = eigen(matrix(s,ncol=2))
		gam = out$vectors
		lam = out$values
		lam = ifelse(lam<0,0,lam)
		xi=atan2(gam[2,1],gam[1,1])
		xi=ifelse(gam[2,1]<0,xi+2*pi,xi)
		ee=sqrt(1-sqrt(lam[1]*lam[2])/(.5*(sqrt(lam[1])+sqrt(lam[2])))^2)
		xii = 0-xi
		ms=rbind(ms,m)
		ss=rbind(ss,s)
		lams=rbind(lams,lam)
		thv=seq(0,2*pi,length=nres) 
		sgv=1+0*thv 
		voln=0 
		ivol=0
		ees = c(ees, ee)
		
		if(length(yu[yu>low])>20) { # Enough points to process
			ivol = 1
			wwt = length(yu[yu>low]) 
			fac = max(yu)
			if(abs(fac)<0.1e-20) {fac=0.1e-20}
			yu = yu/fac
			cv = sqrt(pv*var(yu)*length(yu))
			yp = c(yu,0)
			yu = yu*fac 
			yp = c(yu,0) 
			fac=1
			cc=fitaa(u,yu,alpha)
			sbeta=cc[1]
			tau0=cc[2]
			phi0=cc[3]
			sebeta=sbeta/2
			setau=tau0/2 
			sephi=phi0/2
			tau=tau0
			beta=sbeta
			phi=phi0
			result=try(out <- nls(yp~dgammp(tau,phi,beta,cv,u,alpha),
					control=list(warnOnly =TRUE),
			       	start=list(beta=sbeta,tau = tau0, phi = phi0),
					lower=c(beta=sbeta/2,tau0/1.25,phi0/1.25),
					upper=c(beta=sbeta*2,tau0*1.25,phi0*1.25),
			        trace = FALSE,algorithm = "port"),silent=TRUE)
			if(!inherits(result,"try-error")){
				tau=summary(out)$coef[2,1]
				phi=summary(out)$coef[3,1]
				beta=summary(out)$coef[1,1]
				setau=summary(out)$coef[2,2]
				sephi=summary(out)$coef[3,2]
				sebeta=summary(out)$coef[1,2]
			}
			
			v = tau+phi*u
			siv = v<qgamma(.999,alpha)
			uu = c(uu,v[siv]) 
			yyuu = c(yyuu,(yu[siv]-beta0)/beta)
			yyhat = c(yyhat,
				c(fac*(beta0+beta*dgamma(abs(tau+phi*u),alpha))[siv]))
			yyw = c(yyw,(yu[siv]*0+(beta*fac)^2))
			res = c(res,re[siv])
			k33=round(nb*.33) 
			k66=round(nb*.66) 
			
			if((k==k33)|(k==k66)|(k<0)) {
				yyu=fac*yu
				if(doplot){
					plot(xs[,1],xs[,2],pch="*",
							xlim=range(c(x[,1:2])),ylim=range(c(x[,1:2])),
							main=as.character(round(zv)),
							xlab=" ",ylab=" ",axes=F,col=c(8))
					points(xs[((yyu>y2)&(yyu<=y3)),1],
							xs[((yyu>y2)&(yyu<=y3)),2],pch="*",col=c(4))
					points(xs[((yyu>y3)&(yyu<=y4)),1],
							xs[((yyu>y3)&(yyu<=y4)),2],pch="*",col=c(2))
					points(xs[yyu>y4,1],xs[yyu>y4,2],pch="*",col=c(7)) 
					abline(h=m[2])
					abline(v=m[1])
				}
			}
			igraph=0
			if((k==k33)|(k==k66)) {
				igraph=0
			}
			out = bc(th,u,c(beta0+beta*dgamma(tau+phi*u,alpha)),yu,
						(thres/fac),nres,dfbc,igraph,doplot)
			volu = out$volu
			volc = out$volc
			thv = out$o[,1]
			sgv = out$o[,2]
		}
		sgU=rbind(sgU,sgv)
		sgo=(1-ee)*(1+ee)/sqrt((1-ee)^2*cos(thv+xii)^2+(1+ee)^2*sin(thv+xii)^2)
		sg0=rbind(sg0,sgo)
		zvs=c(zvs,zv)
		r=sgv*sgo*ivol
		voln=sum((r)^2*pi/nres)
		vols=c(vols,voln)
		
		if((k==k33)|(k==k66)|(k<0)) {
			fval=mean(r)/mean(sgo)
			if(doplot){
				lines(m[1]+r*cos(thv),m[2]+r*sin(thv),col=c(5),lwd=c(2))
				plot(u,yu,ylim=range(y),pch=".",axes="F",xlab=" ",ylab=" ",main=" ")
				axis(1)
				lines(unismooth2(u,yu),lwd=c(2),col=c(3))
				abline(h=thres/fac,col=c(2))	
				axis(2)
			}
		}
	}
	
	# Smooth VOLUME
	dev = vols-as.vector(smooth(vols)) 
	dev = (dev-median(dev))/sqrt(var(dev-median(dev))) 
	wtd = 1/(1+dev^2)
	svols = smooth.spline(zvs,vols,wtd,df=max(4,.25*length(zvs)))$y
	low = max(min(vols)+1,mean(vols)-2*sqrt(var(vols))) 
	svols = ifelse(svols>low, svols,low)
	high = min(max(vols), mean(vols)+2*sqrt(var(vols))) 
	svols = ifelse(svols<high,svols,high)
	
	# Smooth sg across slices
	sgUU=sg0*sgU
	for(i in 1:nres) {
		cc=smooth.spline(zvs,sgUU[,i],svols,df=max(5,.5*length(zvs)))$y 
		cc=ifelse(cc>min(sgUU[svols>.05*mean(svols),i]),cc,min(sgUU[svols>.05*mean(svols),i]))
		cc=ifelse(cc<max(sgUU[svols>.05*mean(svols),i]),cc,max(sgUU[svols>.05*mean(svols),i]))
		sgUU[,i]=cc 
	}
	
	# normalizing phases wrt partial volume here... ###
	sg=sgUU 
	for(j in 1:nb) { 
		r=sgUU[j,] 
		volr=sum((r)^2*pi/nres) 
		r=r*sqrt(svols[j]/volr)
		sg[j,]=r 
	}
	# output 
	e <- new.env()
	e$z=zvs
	e$m=ms
	e$vol=svols 
	e$th=thv  
	e$sg=sg
	# add-ons:
	e$ee=ees
	e$sg0=sg0	
	e$s=ss
	e$lams=lams	
	e$siv=siv
	e$res=res
	e$u=u
	e$re=re
	e$uu=uu # actually = v[siv]
	return(e)
} # end of gsg()
		  
# ---------------------------------------------------------------
axb <- function(a,b) { 
	ab = c( a[2]*b[3]-a[3]*b[2] , a[3]*b[1]-a[1]*b[3] , a[1]*b[2]-a[2]*b[1] )
	return(ab)
}

# ---------------------------------------------------------------
viewxv <- function(xv,bxt,nxt,nylab,nxlab,nme,rr,...) {
# In: S_rfuns.R
# 3D meshing of input volume (PAS-transformed VOI)
# Inputs:
# xv = axis of representation (transverse, coronal or sagittal)
# bxt = boundary points (surface)
# nxt = normal points
# rr = range
	#
	# QR decomposition is used to determine visibility?
	qrx = qr(matrix(xv,ncol=1))
	qrx$qr = qrx$qr * sign(qrx$qr[order(-abs(qrx$qr))[1]])
	isee = vis(xv,nxt)*0+1
	rr1 = rr*range(c(t(qr.qty(qrx,t(bxt)))[,2]))
	rr2 = rr*range(c(t(qr.qty(qrx,t(bxt)))[,3]))
	plot(t(qr.qty(qrx,t(bxt)))[,2:3],type="n",xlab=nxlab,ylab=nylab,axes="F",xlim=rr1,ylim=rr2,main=nme,...)
	#
	# slice binning
	sl = as.real(names(table(bxt[,3])))
	nl = length(sl) 
	nld = min(25,ceiling(nl*.75))
	u = approx(sl,c(1:nl),xout=seq(sl[1],sl[nl],length=nld),rule=2,ties="ordered")$y  
	fu = floor(u)
	cu = ceiling(u)
	isl = ifelse(abs(u-fu)<abs(cu-u),fu,cu)
	#
	# slice-by-slice "lining" - what direction?
	for(i in c(1:nld)){ 
		ii=isl[i] 
		iv=isee[abs(bxt[,3]-sl[ii])<.001]
		dline(iv,t(qr.qty(qrx,t(bxt[abs(bxt[,3]-sl[ii])<.001,])))[,2:3],8) 
	}
	#
	nt = ceiling(nrow(bxt)/nl) 
	ntd = min(nt,25)
	for(t in c(1:ntd)){ 
		tt=round(1+(t-1)*(nt-1)/(ntd-1))		
		bxv=NULL 
		iseeb=NULL
		for(i in 1:nl) { 
			bxv=rbind(bxv,(bxt[abs(bxt[,3]-sl[i])<.001,])[tt,])
			iseeb=c(iseeb,isee[abs(bxt[,3]-sl[i])<.001][tt])   
		}
		dline(iseeb,t(qr.qty(qrx,t(bxv)))[,2:3],8) 
	}
	#
	isee=vis(xv,nxt)
	sl = as.real(names(table(bxt[,3])))
	nl = length(sl) 
	nld = min(25,ceiling(nl*.75))
	u = approx(sl,c(1:nl),xout=seq(sl[1],sl[nl],length=nld),rule=2,ties="ordered")$y
	fu = floor(u)
	cu = ceiling(u)
	isl = ifelse(abs(u-fu)<abs(cu-u),fu,cu)
	for(i in c(1:nld)) {
		ii=isl[i] 
		iv=isee[abs(bxt[,3]-sl[ii])<.001]
		dline(iv,t(qr.qty(qrx,t(bxt[abs(bxt[,3]-sl[ii])<.001,])))[,2:3],1) 
	}
	#
	nt=ceiling(nrow(bxt)/nl)
	ntd=min(nt,25)
	for(t in c(1:ntd)) { 
		tt=round(1+(t-1)*(nt-1)/(ntd-1))	
		bxv=NULL 
		iseeb=NULL
		for(i in 1:nl){
			bxv=rbind(bxv,(bxt[abs(bxt[,3]-sl[i])<.001,])[tt,])
			iseeb=c(iseeb,isee[abs(bxt[,3]-sl[i])<.001][tt])   
		}
		dline(iseeb,t(qr.qty(qrx,t(bxv)))[,2:3],1) 
	}
}

# ---------------------------------------------------------------
viewxv2 <- function(xv,bxt,bxo,nxo,nylab,nxlab,nme,rr) {
	qrx=qr(matrix(xv,ncol=1))
	qrx$qr=qrx$qr*sign(qrx$qr[order(-abs(qrx$qr))[1]])
	rr1=rr*range(c(t(qr.qty(qrx,t(bxo)))[,2]))
	rr2=rr*range(c(t(qr.qty(qrx,t(bxo)))[,3]))
	plot(t(qr.qty(qrx,t(bxo)))[,2:3],type="n",
		xlab=nxlab,ylab=nylab,axes="F",xlim=rr1,ylim=rr2,main=nme )
	sl= as.real(names(table(bxt[,3])))
	nl=length(sl) 
	nld=min(25,ceiling(nl*.75))
	isee=vis(xv,nxo)*0+1
	u=approx(sl,c(1:nl),xout=seq(sl[1],sl[nl],length=nld),rule=2,ties="ordered")$y
	fu=floor(u)
	cu=ceiling(u)
	isl=ifelse(abs(u-fu)<abs(cu-u),fu,cu)
	for(i in c(1:nld)) { 
		ii=isl[i] 
		iv=isee[abs(bxt[,3]-sl[ii])<.001]
		dline(iv,t(qr.qty(qrx,t(bxo[abs(bxt[,3]-sl[ii])<.001,])))[,2:3],8) 
	}
	nt=ceiling(nrow(bxt)/nl) 
	ntd=min(nt,25)
	for(t in c(1:ntd)) {
		tt=round(1+(t-1)*(nt-1)/(ntd-1))
		bxv=NULL 
		iseeb=NULL
		for(i in 1:nl) {
			bxv=rbind(bxv,(bxo[abs(bxt[,3]-sl[i])<.001,])[tt,]) 
			iseeb=c(iseeb,isee[abs(bxt[,3]-sl[i])<.001][tt])   
		}
		dline(iseeb,t(qr.qty(qrx,t(bxv)))[,2:3],8) 
	}
	isee=vis(xv,nxo)
	sl= as.real(names(table(bxt[,3])))
	nl=length(sl) 
	nld=min(25,ceiling(nl*.75))
	u=approx(sl,c(1:nl),xout=seq(sl[1],sl[nl],length=nld),rule=2,ties="ordered")$y
	fu=floor(u);cu=ceiling(u);isl=ifelse(abs(u-fu)<abs(cu-u),fu,cu)
	for(i in c(1:nld)) { 
		ii=isl[i] 
		iv=isee[abs(bxt[,3]-sl[ii])<.001]
		dline(iv,t(qr.qty(qrx,t(bxo[abs(bxt[,3]-sl[ii])<.001,])))[,2:3],1) 
	}
	nt=ceiling(nrow(bxt)/nl) 
	ntd=min(nt,25)
	for(t in c(1:ntd)) { 
		tt=round(1+(t-1)*(nt-1)/(ntd-1))
		bxv=NULL
		iseeb=NULL
		for(i in 1:nl) { 
			bxv=rbind(bxv,(bxo[abs(bxt[,3]-sl[i])<.001,])[tt,]) 
			iseeb=c(iseeb,isee[abs(bxt[,3]-sl[i])<.001][tt])   
		}
		dline(iseeb,t(qr.qty(qrx,t(bxv)))[,2:3],1) 
	}
}

# connect contiguous visible arcs on a curve
# ---------------------------------------------------------------
dline <- function(isee,u,color) {
nn=nrow(u) ; is=0 ;ie=0 ; ivis=0  
 while(is<nn) {
	while((is<nn)&(ivis<1)) { is=is+1 ; ivis=isee[is] ; ie=is }
	while((ie<nn)&(ivis>0)) { ie=ie+1 ; ivis=isee[ie]         } 
	ns=ie-is+1 ; if(ns>1) { u1=u[is:ie,1];u2=u[is:ie,2];vis=isee[is:ie]
	lines(u1[vis[1:ns]>0],u2[vis[1:ns]>0],col=c(color)) }
	is=ie
	}
}


# ---------------------------------------------------------------
viewxvp <- function(xv,bxt,nme,rr) {
	qrx=qr(xv)
	qrx$qr=qrx$qr*sign(qrx$qr[order(-abs(qrx$qr))[1]])
	
	plot(t(qr.qty(qrx,t(bxt)))[,2:3],type="n",pch=".",xlab=" ",ylab=" ",axes="F",xlim=range(rr*1.55*c(bxt)),ylim=range(1.55*c(bxt)),main=nme )
	nn=nrow(bxt)/3
	lines(t(qr.qty(qrx,t(bxt)))[c(1:nn),2:3],lwd=c(2))
	text(1.2*t(qr.qty(qrx,t(bxt)))[nn,2],1.2*t(qr.qty(qrx,t(bxt)))[nn,3],"y'")
	
	lines(t(qr.qty(qrx,t(bxt)))[c(nn+c(1:nn)),2:3],lwd=c(2))
	
	text(1.5*t(qr.qty(qrx,t(bxt)))[2*nn,2],1.5*t(qr.qty(qrx,t(bxt)))[2*nn,3],"x'")
	
	lines(t(qr.qty(qrx,t(bxt)))[c(2*nn+c(1:nn)),2:3],lwd=c(2))
	
	text(1.2*t(qr.qty(qrx,t(bxt)))[3*nn,2],1.2*t(qr.qty(qrx,t(bxt)))[3*nn,3],"z'")
}

# ---------------------------------------------------------------

roid <- function(x,y,z,m,sg) { 
#
# ROI and ZUY Data for Heterogeneity Analysis and Synthesis
# Assumes PAS data. 
# Returns a subset of input VOI xt of all voxels with phase u<1.
# Arguments
#	x:  3D PAS coordinates of initial input VOI
#	y:  observed voxel uptake 
#	z:  msg()$z, spine-grid midpoints (actually, output of gsg())
#	m:  msg()$m, 2D spine coordinates (slice-per-slice basis)
#	sg: msg()$sg, initial boundary adjustment for VOI x
# Value (output fields)
#	zuy: triplet (h,u,y) s.t. phase u<1; h = spine-axis values; y = uptake values
#	roi: PAS-domain, (x1,x2) spine-adjusted, 3D coordinates of output VOI
#	r:	 radial coordinates (for convenience)
#	th:	 corresponding angular coordinates (for convenience)
#	u:	 voxel phases (testing: roid()$u==roid()$zuy[,2])
#
	nb=length(z) 
	n=length(y) 
	nv=round(max(1,n/nb)) 
	zb=pasbins(x[,3],nb,nv,n)
	roid=NULL ; zuyd=NULL ; rs=us=ths=NULL
	alli = c(1:nrow(x))
	icut = NULL
	
	for(k in 1:nb){
		ii = ((x[,3]<=zb[k+1])&(x[,3]>zb[k]))
		xs=x[ii,] 
		ys=y[ii]
		
		p.out=polar(xs[,1:2],m[k,],sg[k,])
		u=p.out[,2] 
		r=p.out[,3] 
		th=p.out[,4]
		
		# Adjust for Spine Center
		xs[,1]=xs[,1]-m[k,1]
		xs[,2]=xs[,2]-m[k,2]
		
		tempi = alli[ii]
		icut = c(icut,tempi[which(u<1)])
 		roid = rbind(roid,cbind(xs,ys)[u<1,])
 		ths = c(ths,th[u<1])
		rs = c(rs,r[u<1])
		us = c(us,u[u<1])
		zuyd = rbind(zuyd,cbind(xs[,3],u,ys)[u<1,])	
    }
	e <- new.env()
	e$zuy=zuyd
	e$roi=roid
	e$r=rs
	e$th=ths
	e$u=us
	e$icut=icut
	return(e)
}


# ---------------------------------------------------------------
roid2 <- function(x,y1,y2,e,z,m,sg) { 
# assumes pas data
	nb = length(z) 
	n = length(y1) 
	nv = round(max(1,n/nb)) 
	zb = pasbins(x[,3],nb,nv,n)
	roid=NULL ; zuyd=NULL
	for(k in 1:nb) { 		
		xs = x[((x[,3]<=zb[k+1])&(x[,3]>zb[k])),] 
        y1s = y1[((x[,3]<=zb[k+1])&(x[,3]>zb[k]))]
		y2s = y2[((x[,3]<=zb[k+1])&(x[,3]>zb[k]))]
		u = polar(xs[,1:2],m[k,],sg[k,])[,2]
		# Adjust for Spine Center
		xs[,1] = xs[,1]-m[k,1]
		xs[,2] = xs[,2]-m[k,2]
 		roid = rbind(roid,cbind(xs,y1s,y2s)[(((1-e)<=u)&(u<1)),])
		zuyd = rbind(zuyd,cbind(xs[,3],u,y1s,y2s)[(((1-e)<=u)&(u<1)),])
	}
	e <- new.env()
	e$zuy=zuyd
	e$roi=roid
	return(e)
}

# ---------------------------------------------------------------
izuyd <- function(zuy,nz,nu) {
#
# Image zuy data
#
	z=zuy[,1] ; u=zuy[,2] ; y=zuy[,3]
	n=length(y) ; nv=round(max(1,n/nz)) 
	zb=pasbins(z,nz,nv,n)
	nv=round(max(1,n/nu)) 
	ub=sort(u[c(1,c(1:(nu-1))*nv,n)])
	iy=matrix(rep(NA,nz*nu),ncol=nz)
	for(j in 1:nz) { 
		for(i in 1:nu) { 
			iy[i,j] = median(y[(u<=ub[i+1])&(u>=ub[i])&(z<=zb[j+1])&(z>=zb[j])]) 
		}
	}
	return(iy)
}

# ---------------------------------------------------------------
fitaa <- function(u,y,alpha) { 
#
# FIT GAMMA MODEL USING  APPROXIMATE QUANTILE TRANSFORM
#
# return tau,phi,beta y=yu 
	yp=y[y>0] 
	zs=u[y>0] 
	yp=yp[order(zs)] 
	zs=sort(zs)	# Sample Quantiles
	fits=NULL 
	po=NULL
	for(k in 1:100) { 
		po=c(po,.9*(k)/100)
		fp=cumsum(yp) 
		fp=(fp/max(fp)) 
		np=length(fp) 
		fp=(fp*np+.5)/(np+2) 
		fp=po[k]+(1-po[k])*fp
		zp=qgamma(fp,alpha) 	# Theoretical Quantiles
		tau=qgamma(po[k],alpha)
		out=lm((zp-tau)~0+zs) 
		coef=summary(out)$coef[,1]
		phi=coef[1]
		yh=dgamma(tau+phi*u,alpha)
		beta=sum(yh*y)/sum(yh^2)
		fits=c(fits,sum( (y-beta*yh)^2) ) 
	}
	pov=po[order(fits)[1]] 
	ww=exp(0-.5*fits/(.1e-9+min(fits)))
	pom= mean(po*ww)/mean(ww)
	fp=cumsum(yp) 
	fp=(fp/max(fp)) 
	np=length(fp) 
	fp=(fp*np+.5)/(np+2) 
	fp=pov+(1-pov)*fp
	zp=qgamma(fp,alpha) 					# Theoretical Quantiles
	tau=qgamma(pom,alpha)
	out=lm((zp-tau)~0+zs) 
	coef=summary(out)$coef[,1]
	phi=coef[1]
	yh=dgamma(tau+phi*u,alpha)
	beta=sum(yh*y)/sum(yh^2)
	fits=c(fits,sum( (y-beta*yh)^2) ) 
	        
	return(c(beta,tau,phi,sum( (y-beta*yh)^2)))
}

# ---------------------------------------------------------------
fitg <- function(zuy,nz,nu) { 
#
# FIT Gaussian MODEL
#
	z=zuy[,1]   # h-range values
	uv=zuy[,2]  # voxel phases
	y=zuy[,3]   # uptakes
	pv=.1
	n=length(y) 
	nv=round(max(1,n/nz)) 
	zb=pasbins(z,nz,nv,n) 			# discretized z-grid
	zbb=(zb[-1]+zb[-(nz+1)])/2					# mid-points
	nv=round(max(1,n/nu)) 
	ub=sort(uv)[c(1,c(1:(nu-1))*nv,n)]   		# discretized phase-grid
	ubb=(ub[-1]+ub[-(nu+1)])/2					# mid-points
	
	# pass 1
	betak=NULL; tauk=NULL; phik=NULL; uk=NULL
	sebetak=NULL; setauk=NULL; sephik=NULL; rssk=NULL; nk=NULL				
	ghat=matrix(rep(0,(nu+1)*nz),ncol=nz) 
	mhat=ghat 
	yhat=y
	for(k in 1:nz) {
		yu=y[(z<=zb[k+1])&(z>zb[k])]
		u=uv[(z<=zb[k+1])&(z>zb[k])]
		fac=mean(yu) 
		if(abs(fac)<0.1e-20) {
			fac=0.1e-20
		}
		yu=yu/fac					# thus mean(yu)=1
		yu=yu[order(u)] 
		u=sort(u) 
		out.us=unismooth2(u,yu) 		# smoothing yielding (quasi-)unimodal profile
		ux=seq(min(u),max(u),length=100) 
		yx=approx(u,yu,xo=ux)$y  	# interpolation on ux-grid
		wx=ifelse(yx>0,yx,.1e-7) 	# weight by uptake
		# profile mode estimate:
		tau0 = min((min(u)+max(u))/2,u[order(-out.us[,2])][1])
		# weighted (inverse) scale estimate:
		phi0=1/(1.e-5+sqrt(sum((ux-tau0)*(ux-tau0)*wx)/(.1e-7+sum(wx))))
		# Gaussian profile for "standardized" phase:
		fu=exp(-.5*((u-tau0)*phi0)^2) 
		sbeta=max(0,mean(fu*yu)/mean(fu*fu))
		# rescale regular points taken from spline-profile:
		mhat[,k]=fac*approx(unismooth2(u,yu),xo=ub,rule=2)$y 
		sebeta=abs(sbeta)/2 
		setau=abs(tau0)/2 
		sephi=abs(phi0)/2
		tau=tau0
		beta=sbeta
		phi=phi0
		# Gaussian fit on standardized phase:
		result=try(out <- nls(yu ~ beta*exp(0-.5*( (u-tau)*phi )^2),
				control=list(warnOnly =TRUE),
			    start=list(beta=sbeta,tau = tau0, phi = phi0),
				lower=c(beta=sbeta/2,tau=max(0,tau0-2/phi0),phi=phi0/3),
				upper=c(beta=sbeta*2,tau=tau0+2/phi0,phi=phi0*3),
			    trace = FALSE,algorithm = "port"),silent=TRUE)
		if(!inherits(result,"try-error")){
			tau=summary(out)$coef[2,1]
			phi=summary(out)$coef[3,1];
			beta=summary(out)$coef[1,1]
			setau=summary(out)$coef[2,2]
			sephi=summary(out)$coef[3,2]
			sebeta=summary(out)$coef[1,2]
		}
		v=(ub-tau)*phi
		ghat[,k]=fac*beta*exp(0-.5*v*v)
		v=(u-tau)*phi 
		yuh=fac*beta*exp(0-.5*v*v)    # "final" intermediate Gaussian profile
		yhat[(z<=zb[k+1])&(z>zb[k])]=yuh 
		yu=fac*yu
		betak=c(betak,beta*fac)
		tauk=c(tauk,tau)
		phik=c(phik,phi)
		sebetak=c(sebetak,sebeta*fac)
		setauk=c(setauk,setau)
		sephik=c(sephik,sephi)
		res=c(yu-yuh) 
		rssk=c(rssk,1-(mean(res*res)/var(yu))) 
		nk=c(nk,length(yu))
	}
	## Smooth parameters and recompute
	pars=cbind(betak,tauk,phik) 
	ses=cbind(sebetak,setauk,sephik)
	spars=pars
	for(k in 1:3) { 
		if(1){
			we=1/ifelse(ses[,k]>.001,ses[,k],100) 
		} else {
			we=ifelse(rssk>.001,rssk,0)  
		}
		we=we*we 
		we=we/median(we)
		op=smooth.spline(zbb,(pars[,k]),we,df=min(length(zbb)/2))
		op=smooth.spline(zbb,(pars[,k]),we,df=max(op$df*.7,2))
		spars[,k]=approx(op,xo=zbb,rule=2)$y 
		wmn=sum(spars[,k]*we)/sum(we)
		wvr=sqrt(sum(((spars[,k]-wmn)^2)*we)/sum(we))
		low=max(.01*wmn,wmn-2.5*wvr) 
		high=wmn+2.5*wvr
		if(k==1) { low=.001*wmn }
		spars[,k]=ifelse(spars[,k]>low,spars[,k],low) 
		spars[,k]=ifelse(spars[,k]<high,spars[,k],high)
	}
	
	# pass 2	
	betak=NULL;tauk=NULL ; phik=NULL ; wwk=NULL
	q2k=NULL;q1k=NULL;q3k=NULL
	ghat=matrix(rep(0,(nu+1)*nz),ncol=nz) 
	ynhat=y
	rssk=NULL; nk=NULL 
	vv=y*0
	wf=y*0 
	yn=y
	uk=0*y
	phasem=ghat
	betain=ghat
	for(k in 1:nz) { # yu is NOT SCALED HERE
		yu=y[(z<=zb[k+1])&(z>zb[k])] 
		u=uv[(z<=zb[k+1])&(z>zb[k])]
		fac=1
		tau=spars[k,2]
		beta=spars[k,1]
		phi=spars[k,3]
		result=try(out <- nls(yu ~ beta*exp(0-.5*( (u-tau)*phi )^2),
				control=list(warnOnly =TRUE),
			    start=list(beta=spars[k,1],tau = spars[k,2],phi = spars[k,3]),
				lower=c(beta=spars[k,1]*.75,tau=spars[k,2]*.75,phi=spars[k,3]*.75),
				upper=c(beta=spars[k,1]*1.2,tau=spars[k,2]*1.25,phi=spars[k,3]*1.25),
			    trace = FALSE,algorithm = "port"),silent=TRUE)
		if(!inherits(result,"try-error")){
			tau=summary(out)$coef[2,1]
			phi=summary(out)$coef[3,1]
			beta=summary(out)$coef[1,1] 
		}
		v=(ub-tau)*phi 
		ghat[,k]=fac*beta*exp(0-.5*v*v)
		phasem[,k]=(ub-tau)*phi 
		betain[,k]=beta
		v=(u-tau)*phi 
		yuh=fac*beta*exp(0-.5*v*v)
		ynhat[(z<=zb[k+1])&(z>zb[k])]=yuh/beta 
		yu=fac*yu
		vv[(z<=zb[k+1])&(z>zb[k])]=v 
		wf[(z<=zb[k+1])&(z>zb[k])]=beta^2
		yn[(z<=zb[k+1])&(z>zb[k])]=yu/beta
		q2k=c(q2k,median(v))
		q1k=c(q1k,sort(v)[.25*length(yu)])
		q3k=c(q3k,sort(v)[.75*length(yu)])
		uk[(z<=zb[k+1])&(z>zb[k])]=u
		betak=c(betak,beta)
		tauk=c(tauk,tau)
		phik=c(phik,phi)
		res=c(yu-yuh) 
		rssk=c(rssk,max(0,1-(mean(res*res)/var(yu)))) 
		nk=c(nk,length(yu))
		wwk=c(wwk,sum(yu))
	}
	
	e <- new.env()
	sdk=1/phik 
	muk=0-tauk/sdk # spread and core-phase
	# e$pars:
	#	1: mid-points of z-axis grid
	#	2: number of points evaluated within for each slice
	#	3: RSS of Gaussian fit per slice
	#	4: median-phase
	#	5: Q1-phase
	# 	6: Q3-phase
	#	7: uptake model amplitude (per slice)
	#	8: core phase (mu = tau*phi where tau is initial shift)
	#	9: phase spread (sd = 1/phi)
	#  10: weights (uptake-based)
	e$pars=cbind(zbb,nk,rssk,q2k,q1k,q3k,betak,muk,sdk,wwk)
	e$se=ses 
	e$fit=cbind(vv,wf,yn,yhat) 
	e$im=ghat 
	e$ims=mhat 
	e$phase=phasem 
	e$poten=betain
	### new outputs...
	e$nz = nz  # nb of slices
	e$nk = nk  # nb of points per slice
	e$uk = uk
	e$tauk=tauk
	e$bk=phik
	e$ak=betak
	e$muk=muk
	e$sdk=sdk
	e$vox.phase=vv
	return(e)
}

# ---------------------------------------------------------------
fitth <- function(u,yu,thres) { 
#
# FIT Gaussian MODEL and use it to compute an estimated u for a 
# given threshold and its se (uv,seuv returned)
#
	yu=yu[order(u)] 
	u=sort(u) 
	out=unismooth2(u,yu) 
	n=length(u)
	ux=seq(min(u),max(u),length=100) 
	yx=approx(u,yu,xo=ux,ties="ordered")$y
	tau0 = min((min(u)+max(u))/2,u[order(-out[,2])][1])
	wx=ifelse(yx>0,yx,.1e-7)
	phi0=1/(1.e-5+sqrt(sum((ux-tau0)*(ux-tau0)*wx)/(.1e-7+sum(wx))))
	fu=exp(-.5*((u-tau0)*phi0)^2)
	sbeta=max(0.1e-8,mean(fu*yu)/mean(fu*fu))
	sebeta=abs(sbeta)/2 
	setau=abs(tau0)/2 
	sephi=abs(phi0)/2
	tau=tau0
	beta=sbeta
	phi=phi0
	cov=diag(rep(1,3))
	frac=max(.1e-7,thres/beta)
	frac=min(1,frac) 
	frac=1/frac
	uv=sqrt(abs(2*log(frac)))/max(.1e-7,phi)+tau  
	yhat=beta*exp(0-.5*( (uv-tau)*phi )^2) 
	vuv=n*(max(u)-min(u))
	if(yx[100]>thres) { uv=max(u) ; vuv=n*(max(u)-min(u)) }
	if(max(yx)<thres) { uv=min(u) ; vuv=n*(max(u)-min(u)) }
	result=try(out <- nls(yu ~beta*exp(0-.5*( (u-tau)*phi )^2),
				control=list(warnOnly =TRUE),
	        	start=list(beta=sbeta,tau = tau0, phi = phi0),
				lower=c(beta=sbeta/2,tau=max(0,tau0-2/phi0),phi=phi0/3),
				upper=c(beta=sbeta*2,tau=tau0+2/phi0,phi=phi0*3),
	            trace = FALSE,algorithm = "port"),silent=TRUE)
	if(!inherits(result,"try-error")){
		if(abs(out$convergence)<.1e-8) { 
			tau=summary(out)$coef[2,1]
			phi=summary(out)$coef[3,1]
			beta=summary(out)$coef[1,1]
			setau=summary(out)$coef[2,2]
			sephi=summary(out)$coef[3,2]
			sebeta=summary(out)$coef[1,2]
			cov=(summary(out)$cov)*summary(out)$sigma^2
		}
		frac=max(.1e-7,thres/beta)
		frac=min(1,frac) 
		frac=1/frac
		uv=sqrt(abs(2*log(frac)))/max(.1e-7,phi)+tau 
		uv=min(uv,max(u))
		uv=max(min(u),uv) 
		yhat=beta*exp(0-.5*( (uv-tau)*phi )^2)
		guv = c(1/(sqrt(2*log(frac))*phi*beta),-sqrt(2*log(frac))/phi^2,1)
		vuv = matrix(guv,nrow=1)%*%cov%*%matrix(guv,ncol=1)
		if(abs(out$convergence)>.1e-8)  { 
			vuv=n*(max(u)-min(u)) 
		}
		if(is.na(vuv)) {  #covariance might not be right
			vuv=n*(max(u)-min(u)) 
		} 		
		if(yhat>1.01*thres) { 
			uv=max(u)
			vuv=n*(max(u)-min(u)) 
		}
		if(yhat<0.99*thres) { 
			uv=min(u)
			vuv=n*(max(u)-min(u)) 
		}
	}
	e <- new.env()
	e$uv=uv 
	e$seuv=sqrt(vuv) 
	return(e)
}

as.real <- function(x){
# as.real() became obsolete...
	return(as.double(x))
}

euler.angles <- function(gam,method=0){
# gam is a rotation matrix (eg eigenvectors from covariance matrix)	
# Decomposes Euler angles wrt gam = Rz(xi3) Ry(xi2) Rx(xi1)
# Returns c(xi1,xi2,xi3)
# NB: within extractroi() code we can find:
# 	xi  = atan2(gam[2,1],gam[1,1])
# 	xi = ifelse(gam[2,1]<0,xi+2*pi,xi)
# 
	if(method>0){
		xi2 = -1/sin(gam[3,1])
		if(abs(gam[3,1])!=1){
			xi1 = atan2(gam[3,2]/cos(xi2),gam[3,3]/cos(xi2))
			xi3 = atan2(gam[2,1]/cos(xi2),gam[1,1]/cos(xi2))
		} else {	
			# then abs(theta) = pi/2, must fix xi3
			xi3 = 0
			if(gam[3,1]==1){
				xi1 = xi3 + atan2(gam[1,2],gam[1,3])
			} else {
				xi1 = - xi3 + atan2(-gam[1,2],-gam[1,3])
			}
		}
	} else {
		# code taken from init()
		# Get xi-representation of G (input gam) 
		G = gam
		sg = sign(apply(G,2,sum)) 
		for(j in 1:3) { 
			G[,j]=G[,j]/sg[j] 
		}
		e1=G[,1] 
		xi1=atan2(e1[2],e1[1]) 
		xi2=acos(e1[3])
		q=qr.Q(qr(e1),complete=T)[,2:3]
		xi3=acos(sum(q[,1]*G[,2]))	
		xi = c(xi1,xi2,xi3)
		xi.o = xi # old xi (input)
		# Check : pmat(xi1,xi2,xi3)-G
	}
	return(c(xi1,xi2,xi3))
}
