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
// Performs lookups and evaluates neighborhing values for each voxel in xx.
// This version does NOT use a hull VOI for evaluation at borders.
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
