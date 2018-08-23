ncImpute=function(x, rank.max=2, lambda=0, g=30, type=c("svd","als"), type.thresh=c("MC+","SOFT"), thresh=1e-05, maxit=100, trace.it=FALSE, warm.start=NULL, final.svd=TRUE){
   this.call=match.call()
   type=match.arg(type)
   type.thresh=match.arg(type.thresh)
   fit=ncImpute.x(x,J=rank.max,lambda,g,type,type.thresh,thresh,maxit,trace.it,warm.start,final.svd)
   attr(fit,"call")=this.call
   fit
}
# For objects x of type matrix
ncImpute.x.matrix=function(x,J,lambda,g,type,type.thresh,thresh,maxit,trace.it,warm.start,final.svd){
  if(type=="als") stop("Type 'als' is only available for matrices of class 'Incomplete'. Use type 'svd' for matrices of class 'matrix'.")
  ncimpute.svd(x,J,type.thresh,thresh,lambda,g,maxit,trace.it,warm.start,final.svd)
}
# For objects x of type Incomplete
ncImpute.x.Incomplete=function(x,J,lambda,g,type,type.thresh,thresh,maxit,trace.it,warm.start,final.svd){
  switch(type,
         "als"=ncimputeI.als(x,J,type.thresh,thresh,lambda,g,maxit,trace.it,warm.start,final.svd),
         "svd"=ncimputeI.svd(x,J,type.thresh,thresh,lambda,g,maxit,trace.it,warm.start,final.svd)
         )
}
setGeneric("ncImpute.x",ncImpute.x.matrix)
setMethod("ncImpute.x","Incomplete",ncImpute.x.Incomplete)
