ncimpute.svd <- function (x, J = 2, type.thresh = c("MC+","SOFT"), thresh = 1e-05, lambda = 0, g, maxit = 100, trace.it = FALSE, warm.start = NULL,...){
  this.call=match.call()
  a=names(attributes(x))
  binames=c("biScale:row","biScale:column")
  if(all(match(binames,a,FALSE))){
    biats=attributes(x)[binames]
  } else biats=NULL
  n <- dim(x)
  m <- n[2]
  n <- n[1]
  xnas <- is.na(x)

  nz=m*n-sum(xnas)
  xfill <- x
  xfill[xnas] <- 0
  if(!is.null(warm.start)){
  ### Must have u,d and v components
    if(!all(match(c("u","d","v"),names(warm.start),0)>0)) stop("warm.start does not have components u, d and v")
    D=warm.start$d
    nzD=sum(D>0)
    JD=min(nzD,J)
    U=warm.start$u[,seq(JD)]
    V=warm.start$v[,seq(JD)]
    D=D[seq(JD)]
    xhat=U%*%(D*t(V))
    xfill[xnas] <- xhat[xnas]
  }
  svd.xfill=svd(xfill)
  ratio <- 1
  iter <- 0
  while ((ratio > thresh)&(iter<maxit)) {
    iter <- iter + 1
    svd.old=svd.xfill
    d=svd.xfill$d
    d=Eigen.Threshold(d,lambda,g,type.thresh)
    xhat <- svd.xfill$u[, seq(J)] %*% (d[seq(J)] * t(svd.xfill$v[,seq(J)]))
    xfill[xnas] <- xhat[xnas]
    svd.xfill=svd(xfill)
    ratio=Frob(svd.old$u[, seq(J)],d[seq(J)],svd.old$v[, seq(J)],svd.xfill$u[, seq(J)],Eigen.Threshold(svd.xfill$d,lambda,g,type.thresh)[seq(J)],svd.xfill$v[, seq(J)])
    obj=(.5*sum( (xfill-xhat)[!xnas]^2)+concave.penalty(d,lambda,g,type.thresh))/nz  
    if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
  } 
  d=Eigen.Threshold(svd.xfill$d[seq(J)],lambda,g,type.thresh)
  J=min(sum(d>0)+1,J)
  svd.xfill=list(u=svd.xfill$u[, seq(J)], d=d[seq(J)], v=svd.xfill$v[,seq(J)])
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
  attributes(svd.xfill)=c(attributes(svd.xfill),list(lambda=lambda,call=this.call),biats)
  return(list(fit=svd.xfill, obj=obj, est.rank=length(which(d>1e-09))))
}