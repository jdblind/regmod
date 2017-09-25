#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define nelems(x) (sizeof(x)/sizeof((x)[0]))

//------------------------------------------------------
void print_int_array(int *a, int n)
{
  int i=0;
  for(i=0; i<n; i++)
    printf("%d\t",a[i]);
  printf("\n");
}
void print_double_array(double *a, int n)
{
  int i=0;
  for(i=0; i<n; i++)
    printf("%.2f\t",a[i]);
  printf("\n");
}
void print_double_tab(double *a, double *b, double *c, int n)
{
  int i=0;
  for(i=0; i<n; i++)
    printf("%.2f\t%.2f\t%.2f\n",a[i],b[i],c[i]);
  printf("\n");
}

//------------------------------------------------------
void c_which(double *a, int *b, int na, int *nb, double value)
// Input arguments:
//  - a: the vector of values to search
//  - na: length of a
//  - value: the value to search for within a
// Output arguments are b and nb:
//  - b contains the indices at which value was found in a;
//  - nb is the actual length of b (i.e. the number of times
//    found *value was found in a) a value nb=0 means value 
//    was not found in a.
// NB: output indices (i.e. values in b and value of nb) are defined
// so as to start at 0 (C convention) instead of 1 (R convention).
{
  int i, idx=-1;

  for(i=0; i<na; i++){
    if(a[i]==value){
      idx += 1;
      b[idx] = i;
    }
  }

  // now update length of 'which' output for it to be used in wrapper R function
  *nb = idx+1; 
}

//------------------------------------------------------
void c_lookupxyz(int *ijk, double *x1, double *x2, double *x3, int *n,
		 double *gx, double *gy, double *gz, 		 
		 int *out, int *nout)
// C equivalent to R's lookup.xyz(ijk,x,gx,gy,gz)... 
// Input:
//  - ijk is expected to contain index-values in C convention (0:n-1): ijk = c(i,j,k)
//  - x1, x2, x3 stand for x[,1], x[,2], x[,3]
// Assumes the point exists in c(x1,x2,x3), i.e. that it can be found!
// NB: 
//  * ijk[0] should be found in x1
//  * ijk[1] should be found in x2
//  * ijk[2] should be found in x3
// Output:
//  - out contains the indices at which g[ijk] was found in (x1,x2,x3)
//    (where g = c(gy,gx,gz)). 
//    NB: these indices are in C-convention (i.e. start at 0 and not 1)!!!
//  - nout is the number of such occurrences (0 if no occurrence)
{
  int i=0, idx=-1;

  for(i=0; i<*n; i++){
    if(x1[i] == gy[ijk[0]]){
      if(x2[i] == gx[ijk[1]]){
	if(x3[i] == gz[ijk[2]]){
	  idx += 1;
	  out[idx] = i;
	  break;
	}
      }
    }
  }
  *nout = idx+1;
}

//------------------------------------------------------
void c_lookupijk(double *xyz, 
		 double *gx, double *gy, double *gz, 
		 int *nx, int *ny, int *nz,
		 int *out, int *nout)
// C equivalent to R's lookup.ijk(xyz,gx,gy,gz)... 
// Input:
//  - xyz is the 3D value whose position is searched for within grids gx, 
//    gy and gz (searching for xyz[0] in gy, xyz[1] in gx, xyz[2] in gz).
// Output:
//  - out contains the location at which xyz is found in grid system.
//    out must be of length 3.
//    NB: indices are in C-convention (i.e. start at 0 and not 1)!!!
//  - nout is the number of such occurrences (0 if no occurrence)
{
  int i[1]={-1}, j[1]={-1}, k[1]={-1};
  int ni=0, nj=0, nk=0;

  c_which(gy,i,*ny,&ni,xyz[0]);
  c_which(gx,j,*nx,&nj,xyz[1]);
  c_which(gz,k,*nz,&nk,xyz[2]);

  out[0] = i[0];
  out[1] = j[0];
  out[2] = k[0];
  *nout = (int) ni*nj*nk; // 0 or 1...
}

// ------------------------------------------------------------------------------------------
void c_dlap(int *xyztab, 
		int *v1ptab, int *v1ntab, int *v2ptab, int *v2ntab, int *v3ptab, int *v3ntab, 
		double *LX, int *n, double *dval)
{
  double v=0., v1p=0., v1n=0., v2p=0., v2n=0., v3p=0., v3n=0.;
  int i=0;

  for(i=0; i<*n; i++){
    v = LX[(int) fmin(fmax(xyztab[i],0),(*n-1))];
    v1p = LX[(int) fmax(v1ptab[i],0)];
    v1n = LX[(int) fmax(v1ntab[i],0)];
    v2p = LX[(int) fmax(v2ptab[i],0)];
    v2n = LX[(int) fmax(v2ntab[i],0)];
    v3p = LX[(int) fmax(v3ptab[i],0)];
    v3n = LX[(int) fmax(v3ntab[i],0)];
    dval[i] = - 6*v + v1p+v1n+v2p+v2n+v3p+v3n;
  }
}

//------------------------------------------------------

void c_phase(double *rs, int *n, double *x1s, int *n1s, double *x2s, int *n2s, 
	   double *tau, double *b, double *vals)
// Evaluates phase function (i.e. argument of uptake profile function g)...
// ---
// x1s and x2s are expected to be passed as c(x1s) and c(x2s) from R wrapper function.
// *** alpha is of type integer ***
// vals is the output of interest.
{
  int i=0, j=0;
  double t1=0., t2=0.;

  for(i=0; i<(*n); i++){
    t1 = 0;
    t2 = 0;
    for(j=0; j<(*n1s); j++){
      t1 += x1s[(i+j*(*n))]*tau[j];
    }
    for(j=0; j<(*n2s); j++){
      t2 += x2s[(i+j*(*n))]*b[j];
    }
    vals[i] = t1+rs[i]/t2;
  }
}

void c_g(double *rs, int *n, double *x1s, int *n1s, double *x2s, int *n2s, 
	   double *tau, double *b, int *alpha, double *vals)
// Evaluates profile function g(i.e. uptake function lambda without amplitude term a)...
// ---
// x1s and x2s are expected to be passed as c(x1s) and c(x2s) from R wrapper function.
// *** alpha is of type integer ***
// vals is the output of interest.
{
  int i=0, j=0;
  double t1=0., t2=0., u=0.;

  for(i=0; i<(*n); i++){
    t1 = 0;
    t2 = 0;
    for(j=0; j<(*n1s); j++){
      t1 += x1s[(i+j*(*n))]*tau[j];
    }
    for(j=0; j<(*n2s); j++){
      t2 += x2s[(i+j*(*n))]*b[j];
    }
    u = t1+rs[i]/t2;
    vals[i] = pow(u,((*alpha)-1)) * exp(-u);
  }
}

void c_lam(double *rs, int *n, double *x1s, int *n1s, double *x2s, int *n2s, 
	   double *tau, double *a, double *b, int *alpha, double *vals)
// x1s and x2s are expected to be passed as c(x1s) and c(x2s) from R wrapper function.
// *** alpha is of type integer ***
// vals is the output of interest.
// ...
{
  int i=0, j=0;
  double t1=0., t2=0., t3=0., u=0.;

  for(i=0; i<(*n); i++){
    t1 = 0;
    t2 = 0;
    t3 = 0;
    for(j=0; j<(*n1s); j++){
      t1 += x1s[(i+j*(*n))]*tau[j];
    }
    for(j=0; j<(*n2s); j++){
      t2 += x2s[(i+j*(*n))]*b[j];
      t3 += x2s[(i+j*(*n))]*a[j];
    }
    u = t1+rs[i]/t2;
    vals[i] = t3 * pow(u,((*alpha)-1)) * exp(-u);
  }
}

void c_lam_backup(double *rs, int *n, double *x1s, int *n1s, double *x2s, int *n2s, 
	   double *tau, double *a, double *b, int *alpha, double *vals)
// x1s and x2s are expected to be passed as c(x1s) and c(x2s) from R wrapper function.
// *** alpha is of type integer ***
// vals is the output of interest.
// ...
{
  double b1s[*n][*n1s];
  double b2s[*n][*n2s];
  int i=0, j=0;
  double t1=0., t2=0., t3=0., u=0.;

  // reshape x1s and x2s spline bases...
  for(i=0; i<(*n); i++){
    vals[i] = 0.;
    for(j=0; i<(*n1s); j++){
      b1s[i][j] = x1s[(i+j*(*n))];
    }
  }
  for(j=0; i<(*n2s); j++){
    b2s[i][j] = x2s[(i+j*(*n))];
  }

  for(i=0; i<(*n); i++){
    t1 = 0;
    t2 = 0;
    t3 = 0;
    for(j=0; j<(*n1s); j++){
      t1 += b1s[i][j]*tau[j];
    }
    for(j=0; j<(*n2s); j++){
      t2 += b2s[i][j]*b[j];
      t3 += b2s[i][j]*a[j];
    }
    u = t1+rs[i]/t2;
    vals[i] = t3 * pow(u,((*alpha)-1)) * exp(-u);
  }
}

//------------------------------------------------------

void c_dlam_fast(int *nghbr, double *lx, int *n, double *dvals)
// Evaluates the 3D discrete Laplacian of the uptake model function.
// ...
{
  int i=0;
  for(i=0; i<(*n); i++){
    dvals[i] = -6*lx[i] + lx[(nghbr[(i*6)]-1)] + lx[(nghbr[(i*6)+1]-1)] 
      + lx[(nghbr[(i*6)+2]-1)] + lx[(nghbr[(i*6)+3]-1)] 
      + lx[(nghbr[(i*6)+4]-1)] + lx[(nghbr[(i*6)+5]-1)];
  }
}

//------------------------------------------------------

/* void c_grad_block_a(SEXP rs, SEXP x1s, SEXP x2s,  */
/* 		    SEXP tha, SEXP thb, SEXP tht,  */
/* 		    SEXP mdx, SEXP alpha, SEXP nghbr,  */
/* 		    SEXP lvals, SEXP dvals) */
/* // Evaluates the a-block of (divided-differences) function gradient. */
/* // ... */
/* { */
/*   int i,j; */
/*   int n=length(rs); */
/*   int p=length(tha); */
/*   int n1=length(x1s)/n; */
/*   int n2=length(x2s)/n; */
/*   double eps=0.; */
/*   SEXP bp=PROTECT(allocVector(REALSXP,p)); */
/*   SEXP bn=PROTECT(allocVector(REALSXP,p)); */
/*   SEXP lxp=PROTECT(allocVector(REALSXP,n)); */
/*   SEXP lxn=PROTECT(allocVector(REALSXP,n)); */
/*   SEXP dxp=PROTECT(allocVector(REALSXP,n)); */
/*   SEXP dxn=PROTECT(allocVector(REALSXP,n)); */

/*   // inits */
/*   for(i=0; i<(n); i++){ */
/*     lxp[i]=0.; */
/*     lxn[i]=0.; */
/*     dxp[i]=0.; */
/*     dxn[i]=0.; */
/*   } */
/*   for(i=0; i<(n*p); i++){ */
/*     lvals[i]=0.; */
/*     dvals[i]=0.; */
/*   } */
/*   for(i=0; i<(*p); i++){ */
/*     bp[i]=tha[i]; */
/*     bn[i]=tha[i]; */
/*   } */

/*   // block div diff   */
/*   for(i=0; i<(p); i++){ */
/*     eps = .01*(mdx)/8; */
/*     bp[i] += eps/2; */
/*     bn[i] -= eps/2; */
/*     c_lam(rs,&n,x1s,&n1s,x2s,&n2s,tht,bp,thb,&alpha,lxp); */
/*     //c_dlam_fast(nghbr,lxp,n,dxp); */
/*     //c_lam(rs,n,x1s,n1s,x2s,n2s,tht,bn,thb,alpha,lxn); */
/*     //c_dlam_fast(nghbr,lxn,n,dxn); */
/*     for(j=0; j<(*n); j++){ */
/*       lvals[(i*(*n))+j]=(lxp[i]-lxn[i])/eps; */
/*       dvals[(i*(*n))+j]=(dxp[i]-dxn[i])/eps; */
/*     } */
/*   } */

/*   UNPROTECT(6);  */
/* } */
/* void c_grad_block_a_OBS(double *rs, int *n, */
/* 	double *x1s, int *n1s, double *x2s, int *n2s,   */
/* 	double *tha, double *thb, double *tht, int *p,  */
/* 	double *mdx, int *alpha, int *nghbr,  */
/* 	double *lvals, double *dvals) */
/* // Evaluates the a-block of (divided-differences) function gradient. */
/* // ... */
/* { */
/*   int i=0,j=0; */
/*   double eps=0.; */
/*   double bp[(*p)], bn[(*p)]; */
/*   double lxp[(*n)],dxp[(*n)],lxn[(*n)],dxn[(*n)]; */

/*   // inits */
/*   for(i=0; i<(*n); i++){ */
/*     lxp[i]=0; */
/*     lxn[i]=0; */
/*     dxp[i]=0; */
/*     dxn[i]=0; */
/*   } */
/*   for(i=0; i<((*n)*(*p)); i++){ */
/*     lvals[i]=0; */
/*     dvals[i]=0; */
/*   } */
/*   for(i=0; i<(*p); i++){ */
/*     bp[i]=tha[i]; */
/*     bn[i]=tha[i]; */
/*   } */

/*   // block div diff   */
/*   for(i=0; i<(*p); i++){ */
/*     eps = .01*(*mdx)/8; */
/*     bp[i] += eps/2; */
/*     bn[i] -= eps/2; */
/*     c_lam(rs,n,x1s,n1s,x2s,n2s,tht,bp,thb,alpha,lxp); */
/*     c_dlam_fast(nghbr,lxp,n,dxp); */
/*     c_lam(rs,n,x1s,n1s,x2s,n2s,tht,bn,thb,alpha,lxn); */
/*     c_dlam_fast(nghbr,lxn,n,dxn); */
/*     for(j=0; j<(*n); j++){ */
/*       lvals[(i*(*n))+j]=(lxp[i]-lxn[i])/eps; */
/*       dvals[(i*(*n))+j]=(dxp[i]-dxn[i])/eps; */
/*     } */
/*   } */
/* } */

//------------------------------------------------------

void c_dlam(double *x1, double *x2, double *x3, int *n,
	    double *x1h, double *x2h, double *x3h, int *nh,
	    double *gx, double *gy, double *gz, 
	    int *nx, int *ny, int *nz,
	    double *gxh, double *gyh, double *gzh, 
	    int *nxh, int *nyh, int *nzh,
	    double *LX, double *dval,	    
	    int *ijklupo, double *xyzlupo, int *indlupo, int *viso)
	    //	    int (*ijklup)[3], double (*xyzlup)[7], int (*indlup)[7], int (*vis)[7])
// Evaluates the 3D discrete Laplacian of the uptake model function.
// Latest version which incorporates all preliminary loops...
{
  double xyz[3];
  int xyzout[100];
  int ijk[3];
  int ijkm[3];
  int i=0, j=0;
  int nxyzout=0;
  //  double xyzlup[*n][3];
  int xyzv, xyzv1p, xyzv1n, xyzv2p, xyzv2n, xyzv3p, xyzv3n;
  double v, v1p, v1n, v2p, v2n, v3p, v3n;

  int count=0;

  // reshape args
  int ijklup[*n][3];
  int indlup[*n][7];
  int vis[*n][7];
  double xyzlup[*n][7];

  for(i=0; i<*n; i++){
    ijklup[i][0]=0;
    ijklup[i][1]=0;
    ijklup[i][2]=0;
    for(j=0; j<7; j++){
      indlup[i][j]=0;
      vis[i][j]=0;
      xyzlup[i][j]=-1.;
    }
  }

  for(i=0; i<*n; i++){
    v=0.;
    v1p=0.; v1n=0.; 
    v2p=0.; v2n=0.; 
    v3p=0.; v3n=0.; 

    xyz[0] = x1[i];
    xyz[1] = x2[i];
    xyz[2] = x3[i];

    // look-up ijk (position in (i=gyh,j=gxh,k=gzh) grid)
    ijk[0]=-1; ijk[1]=-1; ijk[2]=-1;
    for(j=0; j<*nyh; j++){
      if(gyh[j]==xyz[0]){
	ijk[0]=j;
	break;
      }
    }
    for(j=0; j<*nxh; j++){
      if(gxh[j]==xyz[1]){
	ijk[1]=j;
	break;
      }
    }
    for(j=0; j<*nzh; j++){
      if(gzh[j]==xyz[2]){
	ijk[2]=j;
	break;
      }
    }
    // for checking purposes:
    ijklup[i][0]=ijk[0];
    ijklup[i][1]=ijk[1];
    ijklup[i][2]=ijk[2];

    // look-up xyz's in augmented volume...
    // init
    for(j=0;j<7;j++)
      xyzlup[i][j]=-1.;

    // xyz
    xyzout[0]=-1; 
    c_lookupxyz(ijk,x1h,x2h,x3h,nh,gxh,gyh,gzh,xyzout,&nxyzout); 
    if(nxyzout>0)
      {
	xyzv = (int) fmax(fmin(xyzout[0],*nh-1),0);
	v = LX[xyzv];
	indlup[i][0] = xyzout[0];
	xyzlup[i][0] = v;
	vis[i][0] = xyzv;
      }
    else
      { 
	count+=1; 
      }
    // xyz.v1p
    for(j=0; j<3; j++)
      ijkm[j] = ijk[j]; 
    ijkm[0] += 1;
    xyzout[0]=-1; 
    c_lookupxyz(ijkm,x1h,x2h,x3h,nh,gxh,gyh,gzh,xyzout,&nxyzout); 
    if(nxyzout>0)
      {
		  xyzv1p = (int) fmax(fmin(xyzout[0],*nh-1),0);
		  v1p = LX[xyzv1p];
		  indlup[i][1] = xyzout[0];
		  xyzlup[i][1] = v1p;
		  vis[i][1] = xyzv1p;
      }
    else
      {
	count+=1; 
      }
    // xyz.v1n
    for(j=0; j<3; j++)
      ijkm[j] = ijk[j]; 
    ijkm[0] -= 1; 
    xyzout[0]=-1; 
    c_lookupxyz(ijkm,x1h,x2h,x3h,nh,gxh,gyh,gzh,xyzout,&nxyzout); 
    if(nxyzout>0)
      {
	xyzv1n = (int) fmax(fmin(xyzout[0],*nh-1),0);
	v1n = LX[xyzv1n];
	indlup[i][2] = xyzout[0];
	xyzlup[i][2] = v1n;
	vis[i][2] = xyzv1n;
      }
    else
      {
	count+=1; 
      }
    // xyz.v2p
    for(j=0; j<3; j++)
      ijkm[j] = ijk[j]; 
    ijkm[1] += 1;
    xyzout[0]=-1; 
    c_lookupxyz(ijkm,x1h,x2h,x3h,nh,gxh,gyh,gzh,xyzout,&nxyzout); 
    if(nxyzout>0)
      {
	xyzv2p = (int) fmax(fmin(xyzout[0],*nh-1),0);
	v2p = LX[xyzv2p];
	indlup[i][3] = xyzout[0];
	xyzlup[i][3] = v2p;
	vis[i][3] = xyzv2p;
      }
    else
      {
	count+=1; 
      }
    // xyz.v2n
    for(j=0; j<3; j++)
      ijkm[j] = ijk[j]; 
    ijkm[1] -= 1;
    xyzout[0]=-1; 
    c_lookupxyz(ijkm,x1h,x2h,x3h,nh,gxh,gyh,gzh,xyzout,&nxyzout);
    if(nxyzout>0)
      {
	xyzv2n = (int) fmax(fmin(xyzout[0],*nh-1),0);
	v2n = LX[xyzv2n];
	indlup[i][4] = xyzout[0];
	xyzlup[i][4] = v2n;
	vis[i][4] = xyzv2n;
      }
    else
      {
	count+=1; 
      }
    // xyz.v3p
    for(j=0; j<3; j++)
      ijkm[j] = ijk[j]; 
    ijkm[2] += 1;
    xyzout[0]=-1; 
    c_lookupxyz(ijkm,x1h,x2h,x3h,nh,gxh,gyh,gzh,xyzout,&nxyzout); 
    if(nxyzout>0)
      {
	xyzv3p = (int) fmax(fmin(xyzout[0],*nh-1),0);
	v3p = LX[xyzv3p];
	indlup[i][5] = xyzout[0];
	xyzlup[i][5] = v3p;
	vis[i][5] = xyzv3p;
      }
    else
      {
	count+=1; 
      }
    // xyz.v3n
    for(j=0; j<3; j++)
      ijkm[j] = ijk[j]; 
    ijkm[2] -= 1;
    xyzout[0]=-1; 
    c_lookupxyz(ijkm,x1h,x2h,x3h,nh,gxh,gyh,gzh,xyzout,&nxyzout);
    if(nxyzout>0)
      {
	xyzv3n = (int) fmax(fmin(xyzout[0],*nh-1),0);
	v3n = LX[xyzv3n];
	indlup[i][6] = xyzout[0];
	xyzlup[i][6] = v3n;
	vis[i][6] = xyzv3n;
      }
    else
      {
	count+=1; 
      }

    // eval dlam...
    dval[i] = - 6*v + v1p+v1n+v2p+v2n+v3p+v3n;

    /* if(dval[i]>1000000){ */
    /*   printf("OHOH at i = %d...\n",i); */
    /*   printf("\txyzv = %d, v1p = %d, v1n = %d, v2p = %d, v2n = %d, v3p = %d, v3n = %d\n", */
    /* 	     xyzv,xyzv1p,xyzv1n,xyzv2p,xyzv2n,xyzv3p,xyzv3n); */
    /*   printf("\tv = %.1f, v1p = %.1f, v1n = %.1f, v2p = %.1f, v2n = %.1f, v3p = %.1f, v3n = %.1f\n", */
    /* 	     v,v1p,v1n,v2p,v2n,v3p,v3n); */
    /* } */
  } 

  // reshape output args
  for(i=0; i<*n; i++){
    ijklupo[i*3] = ijklup[i][0];
    ijklupo[i*3+1] = ijklup[i][1];
    ijklupo[i*3+2] = ijklup[i][2];
    for(j=0; j<7; j++){
      indlupo[i*7+j]=indlup[i][j];
      viso[i*7+j]=vis[i][j];
      xyzlupo[i*7+j]=xyzlup[i][j];
    }
  }
  if(count>0)
    printf("\n\n Final count of failed look-ups: %d out of %d\n\n",count,*n*7);
}

// -------------------------------------------------------- 
void c_galam(double *rs, int *n, double *x1s, int *n1s, double *x2s, int *n2s, 
	   double *tau, double *b, int *alpha, double *vals)
// Evaluates gradient of lambda wrt a...
{
  int i=0, j=0;
  double t1=0., t2=0., u=0.;

  for(i=0; i<(*n); i++){
    t1 = 0;
    t2 = 0;
    for(j=0; j<(*n1s); j++){
      t1 += x1s[(i+j*(*n))]*tau[j];
    }
    for(j=0; j<(*n2s); j++){
      t2 += x2s[(i+j*(*n))]*b[j];
    }
    u = t1+rs[i]/t2;
    for(j=0; j<(*n2s); j++){
	  vals[(i+j*(*n))] = x2s[(i+j*(*n))]*pow(u,((*alpha)-1))*exp(-u);
    }
  }
}

// -------------------------------------------------------- 
void c_gblam(double *rs, int *n, double *x1s, int *n1s, double *x2s, int *n2s, 
	   double *tau, double *a, double *b, int *alpha, double *vals)
// Evaluates gradient of lambda wrt b...
{
  int i=0, j=0;
  double t1=0., t2=0., t3=0., u=0., ax=0.;

  for(i=0; i<(*n); i++){
    t1 = 0;
    t2 = 0;
    t3 = 0;
    for(j=0; j<(*n1s); j++){
      t1 += x1s[(i+j*(*n))]*tau[j];
    }
    for(j=0; j<(*n2s); j++){
      t2 += x2s[(i+j*(*n))]*b[j];
      t3 += x2s[(i+j*(*n))]*a[j];
    }
    u = t1+rs[i]/t2;
	ax = ((*alpha)-1)*pow(u,((*alpha)-2))*exp(-u);
	ax -= pow(u,((*alpha)-1))*exp(-u);
	
    for(j=0; j<(*n2s); j++){
	  vals[(i+j*(*n))] = t3 * ax * (-x2s[(i+j*(*n))]) * (rs[i]/(t2*t2));
    }
  }
}

// -------------------------------------------------------- 
void c_gtlam(double *rs, int *n, double *x1s, int *n1s, double *x2s, int *n2s, 
	   double *tau, double *a, double *b, int *alpha, double *vals)
// Evaluates gradient of lambda wrt tau...
{
  int i=0, j=0;
  double t1=0., t2=0., t3=0., u=0., ax=0.;

  for(i=0; i<(*n); i++){
    t1 = 0;
    t2 = 0;
    t3 = 0;
    for(j=0; j<(*n1s); j++){
      t1 += x1s[(i+j*(*n))]*tau[j];
    }
    for(j=0; j<(*n2s); j++){
      t2 += x2s[(i+j*(*n))]*b[j];
      t3 += x2s[(i+j*(*n))]*a[j];
    }
    u = t1+rs[i]/t2;
	ax = ((*alpha)-1)*pow(u,((*alpha)-2))*exp(-u);
	ax -= pow(u,((*alpha)-1))*exp(-u);
	
    for(j=0; j<(*n1s); j++){
	  vals[(i+j*(*n))] = t3 * ax * (x1s[(i+j*(*n))]);
    }
  }
}

/* -------------------------------------------------------- */
/* -------------------------------------------------------- */

int main()
{
  int NX=27, NXH=125;
  double x1[NX], x2[NX], x3[NX];
  double x1h[NXH], x2h[NXH], x3h[NXH];
  double gx[3], gy[3], gz[3];
  double gxh[5], gyh[5], gzh[5];
  double LX[NXH];
  double dval[27];
  //  int ijklup[27][3];
  //double xyzlup[27][7];
  //int indlup[27][7];
  //int vis[27][7];
  int ijklup[27*3];
  double xyzlup[27*7];
  int indlup[27*7];
  int vis[27*7];
  int nx=3,ny=3,nz=3;
  int nxh=5,nyh=5,nzh=5;
  int i,j;

  //  int ijk[3]={2,2,2}; 
  int *iout=malloc(sizeof(*iout));
  //  int nout=0;  

  for(i=0; i<3; i++){
    gx[i]=1.*(i+2);
    gy[i]=1.*(i+2);
    gz[i]=1.*(i+2);
  }
  for(i=0; i<5; i++){
    gxh[i]=1.*(i+1);
    gyh[i]=1.*(i+1);
    gzh[i]=1.*(i+1);
  }
  for(i=0; i<NX; i++){
    x1[i]=gx[(i%3)];
    x2[i]=gy[((int)floor(i/3))%3];
    x3[i]=gz[((int)floor(i/9))%3];
    dval[i]=0.;
    for(j=0; j<7; j++){
      /* xyzlup[i][j]=0.; */
      /* indlup[i][j]=0; */
      xyzlup[i*7+j]=0.;
      indlup[i*7+j]=0;
    }
  } 
  for(i=0; i<NXH; i++){
    x1h[i]=gxh[(i%5)];
    x2h[i]=gyh[((int)floor(i/5))%5];
    x3h[i]=gzh[((int)floor(i/25))%5];
    LX[i]=2.5*(i+1);
  }

  /* printf("..............................................\n"); */
  /* print_double_tab(gx,gy,gz,3); */
  /* printf("..............................................\n"); */
  /* print_double_tab(gxh,gyh,gzh,5); */
  
  /* printf("..............................................\n"); */
  /* print_double_tab(x1h,x2h,x3h,NXH); */
  /* printf("..............................................\n"); */
  /* print_double_tab(x1,x2,x3,NX); */

  /* c_lookupxyz(ijk,x1,x2,x3,&NX,gx,gy,gz,iout,&nout); */
  /*
  print_int_array(iout, nout);
  printf("value found: %.1f %.1f %.1f\n",x1[iout[0]],x2[iout[0]],x3[iout[0]]);
  printf("TOTAL: %d !!!!!!!!!!!!!!!\n\n",nout);
  */

  /* c_lookupxyz(ijk,x1h,x2h,x3h,&NXH,gxh,gyh,gzh,iout,&nout); */
  /*
  print_int_array(iout, nout);
  printf("value found: %.1f %.1f %.1f\n",x1h[iout[0]],x2h[iout[0]],x3h[iout[0]]);
  printf("TOTAL: %d !!!!!!!!!!!!!!!\n\n",nout); 
  */

  printf("..............................................\n");

  printf("dlam before:\n");
  for(i=0; i<5; i++)
    printf("%.2f\t",dval[i]);

  //  xyzlup[1][1]=-1;
  c_dlam(x1,x2,x3,&NX,x1h,x2h,x3h,&NXH,gx,gy,gz,&nx,&ny,&nz,gxh,gyh,gzh,&nxh,&nyh,&nzh,
	 LX,dval,ijklup,xyzlup,indlup,vis);

  printf("..............................................\n");
  printf("\n XYZLUP:\n");
  for(i=0; i<10; i++){
    for(j=0; j<7; j++){
      /* printf("%.2f\t",xyzlup[i][j]); */
      printf("%.2f\t",xyzlup[i*7+j]);
    }
    printf("\n");
  }
  printf("...\n");

  printf("..............................................\n");
  printf("\n VIS:\n");
  for(i=0; i<10; i++){
    for(j=0; j<7; j++){
      /* printf("%d\t",vis[i][j]); */
      printf("%d\t",vis[i*7+j]);
    }
    printf("\n");
  }
  printf("...\n");

  printf("..............................................\n");
  printf("\n INDLUP:\n");
  for(i=0; i<3; i++){
    for(j=0; j<7; j++){
      /* printf("%d\t",indlup[i][j]); */
      printf("%d\t",indlup[i*7+j]);
    }
    printf("\n");
  }
  printf("...\n");

  printf("..............................................\n");
  printf("\n DVAL:\n");
  printf("\nNX = %d\n",NX);
  for(i=0; i<50; i++){
    printf("%.2f\t",dval[i]);
    if((i+1)%10==0)
      printf("\n");
  }
  printf("...\n");

  printf("\n\n dlam size: %lu\n",sizeof(dval)/sizeof(double));
  //  printf("\n\nTEST %1.2f\n", dval[-1]);

  printf("..............................................\n");
  printf("..............................................\n");

  for(i=0; i<3; i++){
    gx[i]=1.*(i+2);
    gy[i]=1.*(i+2);
    gz[i]=1.*(i+2);
  }
  for(i=0; i<5; i++){
    gxh[i]=1.*(i+1);
    gyh[i]=1.*(i+1);
    gzh[i]=1.*(i+1);
  }
  for(i=0; i<NX; i++){
    x1[i]=gx[(i%3)];
    x2[i]=gy[((int)floor(i/3))%3];
    x3[i]=gz[((int)floor(i/9))%3];
    dval[i]=0.;
    for(j=0; j<7; j++){
      /* xyzlup[i][j]=0.; */
      /* indlup[i][j]=0; */
      xyzlup[i*7+j]=0.;
      indlup[i*7+j]=0;
    }
  } 
  for(i=0; i<NXH; i++){
    x1h[i]=gxh[(i%5)];
    x2h[i]=gyh[((int)floor(i/5))%5];
    x3h[i]=gzh[((int)floor(i/25))%5];
    LX[i]=2.5*(i+1);
  }

  print_double_tab(gx,gy,gz,3);
  print_double_tab(gxh,gyh,gzh,5);

  // look-up ijk (position in (i=gyh,j=gxh,k=gzh) grid)
  int ijk[3];
  int xyzout=0;
  double xyz[3];
  int nh=5;
  int nout=0;

  xyz[0]=3.; 
  xyz[1]=3.; 
  xyz[2]=3.; 

  ijk[0]=-1; ijk[1]=-1; ijk[2]=-1;

  for(j=0; j<5; j++){
    if(gyh[j]==xyz[0]){
      ijk[0]=j;
      break;
    }
  }
  for(j=0; j<5; j++){
    if(gxh[j]==xyz[1]){
      ijk[1]=j;
      break;
    }
  }
  for(j=0; j<5; j++){
    if(gzh[j]==xyz[2]){
      ijk[2]=j;
      break;
    }
  }

  print_double_array(xyz,3);
  printf("....................\n");
  print_int_array(ijk,3);
  printf("....................\n");

  /* print_double_array(x1h,3); */
  nh = 5;
  double X1h[5], X2h[5], X3h[5];
  for(i=0; i<nh; i++){
    X1h[i]=i+1;
    X2h[i]=i;
    X3h[i]=i;
  }
  X1h[4]=3.0; 
  X2h[4]=3.0; 
  X3h[4]=3.0; 

  print_double_tab(X1h,X2h,X3h,nh);
  printf("Looking for (%.1f,%.1f,%.1f) in h-region...\n",gyh[ijk[0]],gxh[ijk[1]],gzh[ijk[2]]);
  printf("Found at ijk = (%d,%d,%d)...\n",ijk[0],ijk[1],ijk[2]);
  c_lookupxyz(ijk,X1h,X2h,X3h,&nh,gxh,gyh,gzh,&xyzout,&nout); 
  printf("Value found in h-region at position %d...\n",xyzout);
  printf("Value found = (%.1f,%.1f,%.1f).\n",X1h[xyzout],X2h[xyzout],X3h[xyzout]);

  printf("\n");
  return 0;
}

