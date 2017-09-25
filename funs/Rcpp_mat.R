library(Rcpp)
library(inline)
library(RcppArmadillo)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### t()...
transCpp = '
using Eigen::Map;
using Eigen::MatrixXd;
const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
const MatrixXd At(A.transpose());
return wrap(At);
'
ftrans <- cxxfunction(signature(AA="matrix"), transCpp, plugin="RcppEigen")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Element-by-element matrix product; works only for square matrices...
prodCpp = '
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
const MapMatd B(as<MapMatd>(BB));
const MapMatd C(as<MapMatd>(CC));
return wrap(B*C);
'
fprod <- cxxfunction(signature(BB="matrix",CC="matrix"), prodCpp, plugin="RcppEigen")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### crossprod()...
xprodCpp = '
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
const MapMatd B(as<MapMatd>(BB));
const MapMatd C(as<MapMatd>(CC));
return wrap(B.transpose()*C);
'
fxprod <- cxxfunction(signature(BB="matrix",CC="matrix"), xprodCpp, plugin="RcppEigen")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### tcrossprod()...
txprodCpp = '
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
const MapMatd B(as<MapMatd>(BB));
const MapMatd C(as<MapMatd>(CC));
return wrap(B*C.transpose());
'
ftxprod <- cxxfunction(signature(BB="matrix",CC="matrix"), txprodCpp, plugin="RcppEigen")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### equiv of crossprod(A), ie using one matrix argument only...
xuprodCpp = '
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
const int n(A.cols());
MatrixXd AtA(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint()));
return wrap(AtA);
'
fxuprod <- cxxfunction(signature(AA="matrix"), xuprodCpp, plugin="RcppEigen")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### equiv of tcrossprod(A), ie using one matrix argument only...
txuprodCpp = '
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
// equiv of tcrossprod(A), ie using one matrix argument only
const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
const int m(A.rows());
MatrixXd AAt(MatrixXd(m, m).setZero().selfadjointView<Lower>().rankUpdate(A));
return wrap(AAt);
'
ftxuprod <- cxxfunction(signature(AA="matrix"), txuprodCpp, plugin="RcppEigen")


#--------------------------------------------------------------------------------
# test...
# NB: code works for matrices of class numeric() (C's "double")
if(0){ ### not run
	csum <- function(x){summary(c(x))}
	A = matrix(c(1:12)*1.0,nr=4,nc=3,byr=F)
	B = matrix(c(1:12)*10.0,nr=4,nc=3,byr=F)
	C = matrix(c(1:6)*.5,nr=3,nc=2,byr=T)
	D = matrix(c(1:9)*1.0,nc=3)
	E = D*10
	F = matrix(rnorm(16),nc=4)
	G = matrix(rnorm(16),nc=4)
	#
	A; B; C;
	t(A); t(B); t(C)
	ftrans(A); ftrans(B); ftrans(C)
	csum(ftrans(A)-t(A))
	#
	# fprod does not work for non-squares matrices:
	csum(fprod(A,B)-A*B)
	# check:
	AB = A*0; for(i in 1:length(A)){AB[i]=A[i]*B[i]}; AB-A*B
	# but works like %*% for square-matrices:
	fprod(D,E)
	D%*%E
	csum(fprod(D,E)-D%*%E)
	csum(fprod(F,F)-F%*%F)
	csum(fprod(F,ftrans(G))-F%*%t(G))
	#
	csum(crossprod(A,B)-fxprod(A,B))
	csum(tcrossprod(A,B)-ftxprod(A,B))
	#
	csum(crossprod(A)-fxuprod(A))
	csum(tcrossprod(A)-ftxuprod(A))
	#
	#
	n = 100; M = 1000
	# or n = 500; M = 100
	A=matrix(rnorm(n*n),nc=n)
	B=matrix(rnorm(n*n),nc=n)
	system.time({for(i in 1:M){fprod(A,B)}})
	system.time({for(i in 1:M){A%*%B}})
}
