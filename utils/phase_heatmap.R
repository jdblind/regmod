phase.heatmap <- function(orph,hv,phiv,th,inds){
# Generates a heatmap of mean phase values...
# Arguments
#	 orph: output of eval.x.rph
#	 hv: vector of h-values (bins), also indicative of number of bins
#	 phiv: vector of phi-values (sectors), also indicative of number of sectors
#	 th: estimates
#	 inds: indices as used throughout
# Usage:
# 	hv = seq(a2,b2,l=J2)
# 	phiv = seq(a1,b1,l=J1)
# 	umap = phase.heatmap(orph,hv,phiv,theta1,all.inds)
# 	heatmap.plus(umap,col=gray((0:255)/255),xlab="phi",ylab="h",Rowv=NA,Colv=NA)
#
	rs=orph$rphs[,1]
	phis=orph$rphs[,2]
	hs=orph$rphs[,3]
	x1s=orph$x1s # h-spline
	x2s=orph$x2s # (phi,h)-spline
	x3s=orph$x3s # phi-spline
	NP = length(phiv)
	NH = length(hv)
	umap = matrix(0,nr=NP-1,nc=NH-1)
	un = c(1:nrow(orph$rphs))
	ih = ip = 1	
	for(ih in 1:(NH-1)){
		for(ip in 1:(NP-1)){
			i1 = which((phis>=phiv[ip])&(phis<phiv[ip+1]))
			deux = un[i1]
			i2 = which((hs[deux]>=hv[ih])&(hs[deux]<hv[(ih+1)]))
			trois = deux[i2] 
			us = 0*trois
			for(iu in 1:length(trois)){
				ti = trois[iu]
				us[iu] = sum(x1s[ti,]*th[inds$tau])+rs[ti]/sum(x2s[ti,]*th[inds$b])
			}
			umap[ip,ih] = mean(us)
		}
	}
	umap[is.na(umap)]=0
	return(umap)
}

get.phases <- function(orph,tau,b,alpha){
	rphs = orph$rphs
	x1s = orph$x1s
	x2s = orph$x2s
	x3s = orph$x3s
	n = nrow(rphs)
	vals = numeric(n)
	for(id in 1:n){
		vals[id] = sum(x1s[id,]*tau) + rphs[id,1]/sum(x2s[id,]*b)
	}
	return(vals)
}

phase.raster <- function(phases,xx){
	NR = length(unique(xx[,2]))
	NC = length(unique(xx[,1]))
	NZ = length(unique(xx[,3]))
	roi = sort.roi(cbind(phases,c(1+0*phases),xx))
	voi = vec.to.raster(roi[,1],nc=NC,nr=NR,dim=NZ)
	return(voi)
}

phase.transverse <- function(orph,xx,th,inds,alpha,xinds){
# Generates mean phase info in transverse domain
	us = get.phases(orph,th[inds$tau],th[inds$b],alpha)
	voi = phase.raster(us[xinds],xx)
	return(apply(voi,c(1,2),mean))
}

