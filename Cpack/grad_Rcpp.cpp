#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void lamRcpp(NumericVector rs, NumericVector x1s, NumericVector x2s, int adim, int bdim,
	     NumericVector tau, NumericVector a, NumericVector b, int alpha, NumericVector vals)
// ...
{
  int i=0, j=0;
  int n = rs.size();
  int n1 = x1s.size()/n;
  int n2 = x2s.size()/n;
  double t1=0., t2=0., t3=0., u=0.;

  for(i=0; i<(n); i++){
    t1 = 0.;
    t2 = 0.;
    t3 = 0.;
    for(j=0; j<(n1); j++){
        t1 += x1s[(i+j*(n))]*tau[j];
        if (adim==1) {
            t3 += x1s[(i+j*(n))]*a[j];
        }
        if (bdim==1) {
            t2 += x1s[(i+j*(n))]*b[j];
        }
    }
    if ((adim+bdim)>2) {
        for(j=0; j<(n2); j++){
            if (adim==2) {
                t3 += x2s[(i+j*(n))]*a[j];
            }
            if (bdim==2) {
                t2 += x2s[(i+j*(n))]*b[j];
            }
        }
    }
    u = t1+rs[i]/t2;
    vals[i] = t3 * pow(u,((alpha)-1)) * exp(-u);
  }
}

// [[Rcpp::export]]
void dlamRcpp(IntegerVector nghbr, NumericVector lx, int n, NumericVector dvals)
// Evaluates the 3D discrete Laplacian of the uptake model function.
// ...
{
  int i=0;
  for(i=0; i<n; i++){
    dvals[i] = -6*lx[i] + lx[(nghbr[(i*6)]-1)] + lx[(nghbr[(i*6)+1]-1)] 
      + lx[(nghbr[(i*6)+2]-1)] + lx[(nghbr[(i*6)+3]-1)] 
      + lx[(nghbr[(i*6)+4]-1)] + lx[(nghbr[(i*6)+5]-1)];
  }
}

// [[Rcpp::export]]
void blocka(NumericVector rs, NumericVector x1s, NumericVector x2s, int adim, int bdim,
	   NumericVector tha, NumericVector thb, NumericVector tht,
	    double mdx, int alpha, IntegerVector nghbr, NumericVector lvals, NumericVector dvals) 
{
  int i,j;
  int n = rs.size();
  int p = tha.size();
  double eps = .01*(mdx)/8;
  NumericVector bp(p); 
  NumericVector bn(p);
  NumericVector lxp(n);
  NumericVector lxn(n);
  NumericVector dxp(n);
  NumericVector dxn(n);

  // inits
  for(i=0; i<(n); i++){
    lxp[i]=0.;
    lxn[i]=0.;
    dxp[i]=0.;
    dxn[i]=0.;
  }
  for(i=0; i<(n*p); i++){
    lvals[i]=0.;
    dvals[i]=0.;
  }
  
  // block div diff  
  for(i=0; i<(p); i++){    
    for(j=0; j<(p); j++){
      bp[j]=tha[j];
      bn[j]=tha[j];
    }
    bp[i] += eps/2;
    bn[i] -= eps/2;
    lamRcpp(rs,x1s,x2s,adim,bdim,tht,bp,thb,alpha,lxp);
    dlamRcpp(nghbr,lxp,n,dxp);
    lamRcpp(rs,x1s,x2s,adim,bdim,tht,bn,thb,alpha,lxn);
    dlamRcpp(nghbr,lxn,n,dxn);
    for(j=0; j<n; j++){
      lvals[(i*n)+j]=(lxp[j]-lxn[j])/eps;
      dvals[(i*n)+j]=(dxp[j]-dxn[j])/eps;
    }
  }
}

// [[Rcpp::export]]
void blockb(NumericVector rs, NumericVector x1s, NumericVector x2s, int adim, int bdim,
	   NumericVector tha, NumericVector thb, NumericVector tht,
	    double mdx, int alpha, IntegerVector nghbr, NumericVector lvals, NumericVector dvals) 
{
  int i,j;
  int n = rs.size();
  int p = thb.size();
  double eps = .01*(mdx)/8;
  NumericVector bp(p); 
  NumericVector bn(p);
  NumericVector lxp(n);
  NumericVector lxn(n);
  NumericVector dxp(n);
  NumericVector dxn(n);

  // inits
  for(i=0; i<n; i++){
    lxp[i]=0.;
    lxn[i]=0.;
    dxp[i]=0.;
    dxn[i]=0.;
  }
  for(i=0; i<(n*p); i++){
    lvals[i]=0.;
    dvals[i]=0.;
  }
  
  // block div diff  
  for(i=0; i<p; i++){    
    for(j=0; j<p; j++){
      bp[j]=thb[j];
      bn[j]=thb[j];
    }
    bp[i] += eps/2;
    bn[i] -= eps/2;
    lamRcpp(rs,x1s,x2s,adim,bdim,tht,tha,bp,alpha,lxp);
    dlamRcpp(nghbr,lxp,n,dxp);
    lamRcpp(rs,x1s,x2s,adim,bdim,tht,tha,bn,alpha,lxn);
    dlamRcpp(nghbr,lxn,n,dxn);
    for(j=0; j<n; j++){
      lvals[(i*n)+j]=(lxp[j]-lxn[j])/eps;
      dvals[(i*n)+j]=(dxp[j]-dxn[j])/eps;
    }
  }
}

// [[Rcpp::export]]
void blockt(NumericVector rs, NumericVector x1s, NumericVector x2s, int adim, int bdim,
	   NumericVector tha, NumericVector thb, NumericVector tht,
	    double mdx, int alpha, IntegerVector nghbr, NumericVector lvals, NumericVector dvals) 
{
  int i,j;
  int n = rs.size();
  int p = tht.size();
  double eps = .01*(mdx)/8;
  NumericVector bp(p); 
  NumericVector bn(p);
  NumericVector lxp(n);
  NumericVector lxn(n);
  NumericVector dxp(n);
  NumericVector dxn(n);

  // inits
  for(i=0; i<n; i++){
    lxp[i]=0.;
    lxn[i]=0.;
    dxp[i]=0.;
    dxn[i]=0.;
  }
  for(i=0; i<(n*p); i++){
    lvals[i]=0.;
    dvals[i]=0.;
  }
  
  // block div diff  
  for(i=0; i<p; i++){    
    for(j=0; j<p; j++){
      bp[j]=tht[j];
      bn[j]=tht[j];
    }
    bp[i] += eps/2;
    bn[i] -= eps/2;
    lamRcpp(rs,x1s,x2s,adim,bdim,bp,tha,thb,alpha,lxp);
    dlamRcpp(nghbr,lxp,n,dxp);
    lamRcpp(rs,x1s,x2s,adim,bdim,bn,tha,thb,alpha,lxn);
    dlamRcpp(nghbr,lxn,n,dxn);
    for(j=0; j<n; j++){
      lvals[(i*n)+j]=(lxp[j]-lxn[j])/eps;
      dvals[(i*n)+j]=(dxp[j]-dxn[j])/eps;
    }
  }
}
