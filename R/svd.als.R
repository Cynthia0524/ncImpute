svd.als=function(x, rank.max=2, lambda=0,thresh = 1e-05, maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE){
  ismiss=is.na(x)
  if(any(ismiss))stop("NAs in x; use softImpute instead")
  this.call=match.call()
  out=ncimpute.als(x,J=rank.max,thresh,lambda,maxit,trace.it,warm.start,final.svd)
  attr(out,"call")=this.call
  attr(out,"lambda")=lambda
  out
}

svd.als.sparse=function(x, rank.max=2, lambda=0,thresh = 1e-05, maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE){
  this.call=match.call()
  out=Ssvd.als(x,J=rank.max,thresh,lambda,maxit,trace.it,warm.start,final.svd)
  attr(out,"call")=this.call
  attr(out,"lambda")=lambda
  out
}

setGeneric("svd.als",svd.als)
setMethod("svd.als","sparseMatrix",svd.als.sparse)
setMethod("svd.als","SparseplusLowRank",svd.als.sparse)
  
