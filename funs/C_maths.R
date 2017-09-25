# dyn.load("./Cpack/libmaths.so")

matmult.C <- function(A,B){
# In: C_maths
# C-version of A%*%B...
# B [n x p] and A [p x n]).
# Returns [p x p] matrix A%*%B.
#
	p = nrow(A)
	n = ncol(A)
	C = numeric(p*p)
	A = t(A)
	B = t(B)
	out <- .C("matmult",X=as.double(c(A)),
						Y=as.double(c(B)),
						Z=as.double(C),
						n=as.integer(n),
						p=as.integer(p))
	return(matrix(out$Z,nc=p,nr=p,byr=T))
}

crossprod.C <- function(A,B=NULL){
# In: C_maths
# C-version of crossprod(), without any checks.
# A and B are expected to be of same dimensions (both [n x p]).
# Returns a [p x p] matrix t(A)%*%B
#
	if(is.null(B)){B=A}
	n = nrow(A)
	p = ncol(A)
	C = numeric(p*p)
	A = t(A)
	B = t(B)
	out <- .C("c_crossprod",X=as.double(c(A)),
						Y=as.double(c(B)),
						Z=as.double(C),
						n=as.integer(n),
						p=as.integer(p))
	return(matrix(out$Z,nc=p,nr=p,byr=T))
}

tcrossprod.C <- function(A,B=NULL){
# In: C_maths
# C-version of tcrossprod(), without any checks.
# A and B are expected to be of same dimensions (both [n x p]).
# Returns a [n x n] matrix A%*%t(B)
#
	if(is.null(B)){B=A}
	n = nrow(A)
	p = ncol(A)
	C = numeric(n*n)
	A = t(A)
	B = t(B)
	out <- .C("c_tcrossprod",X=as.double(c(A)),
						Y=as.double(c(B)),
						Z=as.double(C),
						n=as.integer(n),
						p=as.integer(p))
	return(matrix(out$Z,nc=n,nr=n,byr=T))
}

#----------- testing...
if(0){ ### not run
	n=4
	p=3
	A = matrix(c(1:12),nr=n,nc=p,byr=T)
	B = matrix((2*c(0:11)+1),nr=n,nc=p,byr=T)
	C = t(A)
	
	C%*%B
	matmult.C(C,B)
	
	# t(A)
	t(A)%*%B
	crossprod(A,B)
	A%*%t(B)
	tcrossprod(A,B)
	
	summary(c(crossprod(A,B)-crossprod.C(A,B)))
	summary(c(tcrossprod(A,B)-tcrossprod.C(A,B)))
}