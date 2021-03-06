source("funs/devfuns.R")

# Set global cut-off uptake value:
cp = .25
(cutoff = as.numeric(quantile(z,cp)))

phig = seq(phia,phib,l=Jphi)
hg = seq(ab[1],ab[2],l=Jh)

rs0 = matrix(0,nc=Jphi,nr=Jh)
rs1 = matrix(0,nc=Jphi,nr=Jh)
ap = matrix(NA,nc=Jphi,nr=Jh)
bp = matrix(NA,nc=Jphi,nr=Jh)
tp = matrix(NA,nc=Jphi,nr=Jh)
for(j in 1:Jphi){
	for(k in 1:Jh){
		phi = phig[j]
		h = hg[k]
		bnd.0 = bnd.r(theta.init,h,phi)
		bnd.1 = bnd.r(theta.final,h,phi)
		rs0[k,j] = bnd.0$r.tgt
		rs1[k,j] = bnd.1$r.tgt
		ap[k,j] = bnd.1$ap
		bp[k,j] = bnd.1$bp
		tp[k,j] = bnd.1$tp
	}
}

# Jh
zs  = sort(unique(xx[,3]))
ms  = matrix(0,nr=Jh,nc=2)
ths = seq(0,2*pi,len=Jphi)
# sg:
brads <- function(th){
	rsd = matrix(0,nc=Jphi,nr=Jh)
	for(j in 1:Jphi){
		for(k in 1:Jh){
			phi = phig[j]
			h = hg[k]
			bnd.d = bnd.r(th,h,phi)
			rsd[k,j] = bnd.d$r.tgt
		}
	}
	return(rsd)
}
# sg = rsd
# sg = msg.out$sg

## original volume
bxt = bpts(zs,ms,ths,brads(theta.true))
nxt = npts(zs,ms,ths,brads(theta.true))
bgam0 = bpts(zs,ms,ths,brads(theta.init))
ngam0 = npts(zs,ms,ths,brads(theta.init))
bgam1 = bpts(zs,ms,ths,brads(theta.final))
ngam1 = npts(zs,ms,ths,brads(theta.final))

bn0 = prune.xv(bgam0,ngam0)
bn1 = prune.xv(bgam1,ngam1)
bgam0 = bn0$bgam
ngam0 = bn0$ngam
bgam1 = bn1$bgam
ngam1 = bn1$ngam

if(0){
	par(mfrow=c(3,3))
	#
	# initial extractroi() volume
	xv=c(0,1,0);rrx=1.1
	viewxv(xv,cbind(bxt[,1],bxt[,2],bxt[,3]),cbind(nxt[,1],nxt[,2],nxt[,3]),
			"x3'","x1'",paste(ptid,"(PAS), cutoff",round(cutoff,2)),rrx) 
	xv=c(1,0,0);rrx=1.1
	viewxv(xv,cbind(bxt[,1],bxt[,2],bxt[,3]),cbind(nxt[,1],nxt[,2],nxt[,3]),
			"x3'","x2'"," ",rrx)
	xv=c(0,0,1);rrx=1.1
	viewxv(xv,cbind(bxt[,1],bxt[,2],bxt[,3]),cbind(nxt[,1],nxt[,2],nxt[,3]),
			"x2'","x1'"," ",rrx)
	#
	# initial Gamma volume
	xv=c(0,1,0);rrx=1.1
	viewxv(xv,cbind(bgam0[,1],bgam0[,2],bgam0[,3]),cbind(ngam0[,1],ngam0[,2],ngam0[,3]),
			"x3'","x1'",paste(ptid,"in, gam",gam),rrx) 
	xv=c(1,0,0);rrx=1.1
	viewxv(xv,cbind(bgam0[,1],bgam0[,2],bgam0[,3]),cbind(ngam0[,1],ngam0[,2],ngam0[,3]),
			"x3'","x2'"," ",rrx)
	xv=c(0,0,1);rrx=1.1
	viewxv(xv,cbind(bgam0[,1],bgam0[,2],bgam0[,3]),cbind(ngam0[,1],ngam0[,2],ngam0[,3]),
			"x2'","x1'"," ",rrx)
	#
	# final Gamma volume
	xv=c(0,1,0);rrx=1.1
	viewxv(xv,cbind(bgam1[,1],bgam1[,2],bgam1[,3]),cbind(ngam1[,1],ngam1[,2],ngam1[,3]),
			"x3'","x1'",paste(ptid,"out, gam",gam),rrx) 
	xv=c(1,0,0);rrx=1.1
	viewxv(xv,cbind(bgam1[,1],bgam1[,2],bgam1[,3]),cbind(ngam1[,1],ngam1[,2],ngam1[,3]),
			"x3'","x2'"," ",rrx)
	xv=c(0,0,1);rrx=1.1
	viewxv(xv,cbind(bgam1[,1],bgam1[,2],bgam1[,3]),cbind(ngam1[,1],ngam1[,2],ngam1[,3]),
			"x2'","x1'"," ",rrx)
}
