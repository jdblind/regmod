get.phases <- function(theta,z,orph,xinds,alpha){
####
	ax2s = (modcon$axs==2)
	r1 = orph$rphs[xinds,1]
	rphs = orph$rphs[xinds,]
	x1s = orph$x1s[xinds,]
	x2s = orph$x2s[xinds,]
	x3s = orph$x3s[xinds,]
	n = length(z)
	ps = vs = gh = za = zh = numeric(n)
	for(i in 1:n){
		aterm = ifelse(modcon$axs==1,sum(x1s[i,]*theta[all.inds$a]),sum(x2s[i,]*theta[all.inds$a]))
		bterm = ifelse(modcon$bxs==1,sum(x1s[i,]*theta[all.inds$b]),sum(x2s[i,]*theta[all.inds$b]))
		ps[i] = sum(x1s[i,]*theta[all.inds$tau]) + r1[i]/bterm
		gh[i] = ps[i]^(alpha-1)*exp(-ps[i])
		za[i] = z[i]/aterm		
		zh[i] = aterm * ps[i]^(alpha-1)*exp(-ps[i])
	}
	return(list(phases=ps,ghat=gh,z=z,za=za,zhat=zh))
}

plot.profile <- function(gp,titre=NULL,xlimi=NULL){
# gp is the output of get.phases()
	if(is.null(titre))
		titre = "amplitude-adjusted profile"
	plot(gp$phases,gp$za,pch=20,cex=.7,xlim=xlimi,
		xlab="spline-phases",
		ylab="g-uptake (za, using splines)",
		main=titre)
	points(gp$phases,gp$gh,col=4,pch=20,cex=.6)
}

plot.fit <- function(gp,titre=NULL,xlimi=NULL){
# gp is the output of get.phases()
	if(is.null(titre))
		titre = "profile fit"
	plot(gp$phases,gp$z,pch=20,cex=.7,xlim=xlimi,
		xlab="spline-phases",
		ylab="uptake fit",
		main=titre)
	points(gp$phases,gp$zhat,col=4,pch=20,cex=.6)
}

plot.phases <- function(theta,z,orph,xinds,alpha){
	plot.profile(get.phases(theta,z,orph,xinds,alpha))
}

