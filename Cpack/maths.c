#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void t(double *A, double *tA, int *n, int *p)
// Standard matrix product: given A [n x p], 
// returns tA = t(A) [p x n].
{
  int i=0, j=0;

  for(i=0; i<(*p); i++){
    for(j=0; j<(*n); j++){
      tA[(((i)*(*n))+j)] = A[(((j)*(*p))+i)];
    }
  }
} 

void matmult(double *A, double *B, double *C, int *n, int *p)
// Standard matrix product: given A [p x n] and B [n x p], 
// returns C = A*B [p x p].
{
  int i,j,k;
  int a=0, b=0;
  double temp=0.;

  for(i=0; i<(*p); i++){
    for(j=0; j<(*p); j++){
      temp = 0.;
      for(k=0; k<(*n); k++){
		a = ((i)*(*n)+k);
		b = ((k)*(*p)+j);
        temp += ((A[a])*(B[b]));
      }
      C[((i)*(*p)+j)] = temp;
    }
  }
} 

void c_crossprod(double *X, double *Y, double *Z, int *n, int *p)
// Standard matrix product: given X [n x p] and Y [n x p], 
// returns Z = t(X)*Y [p x p].
{
  double tX[((*n)*(*p))];

  t(X,tX,n,p);
  matmult(tX,Y,Z,n,p);
} 

void c_tcrossprod(double *X, double *Y, double *Z, int *n, int *p)
// Standard matrix product: given X [n x p] and Y [n x p], 
// returns Z = X*t(Y) [p x p].
{
  double tY[((*n)*(*p))];

  t(Y,tY,n,p);
  matmult(X,tY,Z,p,n);
} 

/* -------------------------------------------------------- */
/* -------------------------------------------------------- */

int main()
{
  int n=4, p=3;
  double A[(n*p)], B[(n*p)], C[(p*p)], D[(n*n)], tA[(n*p)];
  int i=0;
  int count=0;

  for(i=0; i<(n*p); i++){
    A[i]=(i+1)*1.0;
	B[i]=(2*i+1)*1.0;
  }
  for(i=0; i<(p*p); i++){
	C[i]=0.0;
  }

  printf("A ....................\n");
  for(i=0; i<(p*n); i++){
	printf("%f\t",A[i]);
    count += 1;
	if(count==p){
      printf("\n");
      count=0;
    }
  }

  printf("B ....................\n");
  for(i=0; i<(p*n); i++){
	printf("%f\t",B[i]);
    count += 1;
	if(count==p){
      printf("\n");
      count=0;
    }
  }

  t(A,tA,&n,&p);

  printf("tA ....................\n");
  for(i=0; i<(p*n); i++){
	printf("%f\t",tA[i]);
    count += 1;
	if(count==n){
      printf("\n");
      count=0;
    }
  }

  matmult(tA,B,C,&n,&p);

  printf("tA*B ....................\n");
  for(i=0; i<(p*p); i++){
	printf("%f\t",C[i]);
    count += 1;
	if(count==p){
      printf("\n");
      count=0;
    }
  }

  c_crossprod(A,B,C,&n,&p);
	printf("C after xprod: \n");
	for(i=0; i<((p)*(p)); i++){
	  printf("%f\t",C[i]);
	  count += 1;
	  if(count==(p)){
	     printf("\n");
	     count=0;
	  }
	}

  c_tcrossprod(A,B,D,&n,&p);
	printf("D after txprod: \n");
	for(i=0; i<((n)*(n)); i++){
	  printf("%f\t",D[i]);
	  count += 1;
	  if(count==(n)){
	     printf("\n");
	     count=0;
	  }
	}
  
  printf("\n");
  return 0;
}

