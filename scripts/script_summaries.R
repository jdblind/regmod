# ------------------------------------------------------------------------ profile

th1 = theta1
ps0.r = ps1.r = zh0.r = zh1.r = gh0 = gh1 = numeric(n)
aas = aas.0 = betas = betas.0 = taus = taus.0 = numeric(n)
for(i in 1:n){
	bterm = sum(x2s[i,]*th1[b.inds])
	ps0.r[i] = sum(x1s[i,]*th0[tau.inds]) + rs[i]/sum(x2s[i,]*th0[b.inds])
	ps1.r[i] = sum(x1s[i,]*th1[tau.inds]) + rs[i]/bterm
	gh0[i] = ps0.r[i]^(alpha-1)*exp(-ps0.r[i])
	gh1[i] = ps1.r[i]^(alpha-1)*exp(-ps1.r[i])
	zh0.r[i] = gh0[i]*sum(x2s[i,]*th0[a.inds])
	zh1.r[i] = gh1[i]*sum(x2s[i,]*th1[a.inds])
	betas[i] = bterm
	betas.0[i] = sum(x2s[i,]*th0[b.inds])
	aas.0[i] = sum(x2s[i,]*th0[a.inds])
	aas[i] = sum(x2s[i,]*th1[a.inds])
	taus.0[i] = sum(x1s[i,]*th0[tau.inds])
	taus[i] = sum(x1s[i,]*th1[tau.inds])
}

abT = hab(xt)
rphT = orph
gpT = get.phases(theta.true,z,rphT,xinds,alpha)
gp0 = get.phases(th0,z,orph,xinds,alpha)
gp1 = get.phases(theta1,z,orph,xinds,alpha)
# --- images ---
roi.zT = rasterize.voi(sort.roi(cbind(z,c(1+0*z),xx)),def=NA)
roi.z0 = rasterize.voi(sort.roi(cbind(gp0$zhat,c(1+0*z),xx)))
roi.z1 = rasterize.voi(sort.roi(cbind(gp1$zhat,c(1+0*z),xx)))
roi.pT = rasterize.voi(sort.roi(cbind(gpT$phases,c(1+0*z),xx)))
roi.p1 = rasterize.voi(sort.roi(cbind(gp1$phases,c(1+0*z),xx)))
roi.pTa = rasterize.voi(sort.roi(cbind(-abs(2-gpT$phases),c(1+0*z),xx)))
roi.p1a = rasterize.voi(sort.roi(cbind(-abs(2-gp1$phases),c(1+0*z),xx)))
a.true = matrix(theta.true[a.inds],nc=Jh,nr=Jphi,byr=T)
a.in = matrix(theta.init[a.inds],nc=Jh,nr=Jphi,byr=T)
a.out = matrix(theta1[a.inds],nc=Jh,nr=Jphi,byr=T)
as.in = rasterize.voi(cbind(aas.0,c(1+0*z),xx))
as.out = rasterize.voi(cbind(aas,c(1+0*z),xx))
pmask = 0*gp1$phases
if(0){
	pmask[gp1$phases>10] = .1
	pmask[gp1$phases<6] = .3
	pmask[gp1$phases<4] = .5
	pmask[gp1$phases<3.2] = .85
	pmask[gp1$phases<2.5] = 1
	pmask[gp1$phases<1.3] = .45
} else {
	pmask = -abs(gp1$phases-2)
}
roi.pmask = rasterize.voi(sort.roi(cbind(pmask,c(1+0*z),xx)))
tau.pmask = rasterize.voi(sort.roi(cbind(taus,c(1+0*z),xx)))

# ------------------------------------------------------------------------ quantitative analysis

# --- hets and phases ---
nb = length(eroi$z.midpoints)
n = nrow(xt)
nv = round(max(1,n/nb)) 
zb = eroi$z.pasbins
zz = xt[,3]
res0 = res1 = resT = 0*zz
hets0 = hets1 = hetsT = mp = core.mp = numeric(Jh)
for(k in 1:Jh){
	kinds = which((zz<=zb[k+1])&(zz>zb[k]))
	resT[kinds] = z[kinds]-gpT$zhat[kinds]
	res0[kinds] = z[kinds]-gp0$zhat[kinds]
	res1[kinds] = z[kinds]-gp1$zhat[kinds]
	hetsT[k] = max(0, 1-(mean(resT[kinds]^2)/var(z[kinds])))
	hets0[k] = max(0, 1-(mean(res0[kinds]^2)/var(z[kinds])))
	hets1[k] = max(0, 1-(mean(res1[kinds]^2)/var(z[kinds])))
	mp[k] = median(gp1$phases[kinds])
}
hetsT = 100-100*hetsT
hets0 = 100-100*hets0
hets1 = 100-100*hets1

# --- radial summaries (uptake and phase) ---
Nrbins = 21
rbins = seq(0,max(rs),len=(Nrbins+1))
uptk = phs = matrix(NA,nr=Jh,nc=Nrbins)
for(k in 1:Jh){
	kinds = which((zz<=zb[k+1])&(zz>zb[k]))
	if(length(kinds)){
		xzk = cbind(z[kinds],rs[kinds],gp1$phases[kinds],xt[kinds,])
		for(b in 1:Nrbins){
			ri = which((xzk[,2]>=rbins[b])&(xzk[,2]<rbins[(b+1)]))
			if(length(ri)){
				uptk[k,b] = mean(xzk[ri,1])
				phs[k,b] = mean(xzk[ri,3])
			}
		}
	}
}

# --- phases ---
#
# h-binning...
bins = numeric(n)
for(k in 1:Jh){
	kinds = which((zz<=zb[k+1])&(zz>zb[k]))
	bins[kinds] = k
}
#
# phi-binning...
phis = seq(phia,phib,len=(Jphi+1))
pbins = numeric(n)
for(k in 1:Jphi){
	kinds = which((rphs[,2]<=phis[k+1])&(rphs[,2]>phis[k]))
	pbins[kinds] = k
}
#
# work out mean phases per (phi,h)...
mp2 = numeric(Jh)
phases.out = matrix(0,nc=Jphi,nr=Jh)
for(k in 1:Jh){
	ks = which(bins==k)
	mp2[k] = mean(gp1$phases[ks])
	for(j in 1:Jphi){
		js = which(pbins[ks]==j)
		phases.out[k,j] = mean(c(gp1$phases[ks])[js])
	}
}
