all.inds = get.inds(Jphi,Jh)
a.inds = all.inds$a
tau.inds = all.inds$tau
b.inds = all.inds$b
c.inds = all.inds$c
s.inds = all.inds$s
xi.inds = all.inds$xi
mu1.inds = all.inds$mu1
mu2.inds = all.inds$mu2

# structural parameters
io = init.gpars(eroi,g2g)
a.init = io$a0
b.init = io$b0
tau.init = io$t0
th0 = numeric(length(unlist(all.inds)))
th0[c.inds]   = eroi$c
th0[s.inds]   = eroi$s
th0[xi.inds]  = eroi$xi
th0[mu1.inds] = c(matrix(eroi$mu1,nc=1))
th0[mu2.inds] = c(matrix(eroi$mu2,nc=1))
if(modcon$axs==1){
	th0[a.inds]   = eroi$a*fac
} else {
	th0[a.inds]   = rep(eroi$a,each=Jphi)*fac
}
if(modcon$bxs==1){
	th0[b.inds]   = c(apply(matrix(b.init,nc=12),2,'mean'))
} else {
	th0[b.inds]   = b.init
}
th0[tau.inds] = tau.init

# ---------------------------------------------------------------

nb = length(eroi$z.midpoints)
n = length(z)
nv = round(max(1,n/nb)) 
# zb = sort(xt[,3])[c(1,c(1:(nb-1))*nv,n)] 
# zb[1] = zb[1]-max(.1e-9,.001*(zb[2]-zb[1])) 	# pas bins
zb=pasbins(xt[,3],nb,nv,n)
zb.new = zb 
zb = eroi$z.pasbins
zz = xt[,3]

rs = us = phis = numeric(n)
muk = matrix(0,nr=n,nc=2)
for(k in 1:nb){
	kinds = which((zz<=zb[k+1])&(zz>zb[k]))
	p.out = polar(xt[kinds,1:2],c(0,0),eroi$msg.out$sg[k,],thcorr=2)
	us[kinds] = p.out[,2]
	rs[kinds] = p.out[,3]
	phis[kinds] = p.out[,1]
	muk[kinds,] = matrix(c(eroi$mu1[k],eroi$mu2[k]),nc=2,nr=length(kinds),byr=T)
}
rphs = cbind(rs,phis,xt[,3])

ab = range(rphs[,3])
x1s=matrix(0,nr=n,nc=Jh)
x2s=matrix(0,nr=n,nc=(Jh*Jphi))
x3s=matrix(0,nr=n,nc=Jphi)		
x = unlist(xt)
for(i in 1:n) {
	phi = rphs[i,2]
	h = x[i,3]
	x1 = ecbs(h,ab[1],ab[2],Jh)
	x1s[i,]=x1
	x2 = etpb12(c(phi,h),phia,phib,Jphi,ab[1],ab[2],Jh) 
	x2s[i,]=x2
	x3 = epcbs(phi,phia,phib,Jphi)
	x3s[i,]=x3
}
orph = list(rphs=rphs,x1s=x1s,x2s=x2s,x3s=x3s)
summary(apply(orph$x1s,1,sum))
summary(apply(orph$x2s,1,sum))
summary(apply(orph$x3s,1,sum))

beta0 = eroi$b*g2g$sc.adj
b0 = rep(beta0,each=Jphi)

a0vec = b0vec = tau0vec = numeric(n)
for(k in 1:nb){
	kinds = which((zz<=zb[k+1])&(zz>zb[k]))
	tau0vec[kinds] = eroi$tau[k]
	a0vec[kinds] = c(eroi$a)[k]
	b0vec[kinds] = beta0[k]
}

ps.r = zh.r = zh.ar = numeric(n)
for(i in 1:n){
	bterm=ifelse(modcon$bxs==1,sum(x1s[i,]*th0[b.inds]),sum(x2s[i,]*th0[b.inds]))
	aterm=ifelse(modcon$axs==1,sum(x1s[i,]*th0[a.inds]),sum(x2s[i,]*th0[a.inds]))
	ps.r[i] = sum(x1s[i,]*th0[tau.inds]) + rs[i]/bterm
	zh.r[i] = ps.r[i]^(alpha-1)*exp(-ps.r[i])*a0vec[i]
	zh.ar[i] = ps.r[i]^(alpha-1)*exp(-ps.r[i])*aterm
}
