ncimputeI.als <- function (x, J = 2, type.thresh = c("MC+","SOFT"), thresh = 1e-05, lambda = 0, g, maxit = 100, trace.it = FALSE, warm.start = NULL,...){
### This function expects an object of class "Incomplete" which inherits from "sparseMatrix", where the missing entries
### are replaced with zeros. If it was centered, then it carries the centering info with it. Corresponding singular value 
### decompositions are carried out using the rank-restricted (soft) SVD method Ssvd.als from the R package 'softImpute'
### with a default threshold value of lambda=0
  this.call=match.call()
  a=names(attributes(x))
  binames=c("biScale:row","biScale:column")
  if(all(match(binames,a,FALSE))){
    biats=attributes(x)[binames]
  } else biats=NULL

  if(!inherits(x,"dgCMatrix"))x=as(x,"dgCMatrix")
  ### Row-indices for Incomplete Matrix
  irow=x@i
  ### Column indices for Incomplete Matrix in column-compressed format
  pcol=x@p
  n <- dim(x)
  m <- n[2]
  n <- n[1]
  ### Number of nonzero values
  nz=nnzero(x)
  ### Copy of Incomplete Matrix
  xres=x
   if(!is.null(warm.start)){
    ###must have u,d and v components
    if(!all(match(c("u","d","v"),names(warm.start),0)>0))stop("warm.start does not have components u, d and v")
    D=warm.start$d
    nzD=sum(D>0)
    ### Minimum between specified operating rank and rank provided by the warm started solution
    JD=min(nzD,J)
    U=warm.start$u[,seq(JD)]
    V=warm.start$v[,seq(JD)]
    D=D[seq(JD)]
    if (JD==1) U = matrix(U, length(U), 1)
    ### Forms matrix V*D
    BD=UD(V,D,m)
    ### Computes (U %*% t(BD))[i,j] = sum(u[i,]*v[j,]) for all pairs in irow,pcol
    xhat.Omega=suvC(U,BD,irow,pcol)
    ### Forming Sparse Part of a Sparse plus Low-Rank Matrix
    xres@x=x@x-xhat.Omega
    ### Forming object of type Sparse plus Low-Rank. Sparse: xres; Low-Rank: a=U, b=BD=VD
    xfill=splr(xres,U,BD)
    ### Computes rank-restricted (soft-thresholded) SVD of xfill with warm-start (U,D,V) and default lambda=0
    svd.xfill=Ssvd.als(xfill,J,thresh,lambda=0,maxit,trace.it=FALSE,warm.start=list(u=U,d=D,v=V))
    }
  else{
    xfill=x
    ### Computes rank-restricted (soft-thresholded) SVD of xfill without warm-starts and default lambda=0
    svd.xfill=Ssvd.als(xfill,J,thresh,lambda=0,maxit,trace.it=FALSE)
  }
  ratio <- 1
  iter <- 0
  while ((ratio > thresh)&(iter<maxit)) {
    iter <- iter + 1
    svd.old=svd.xfill
    ### MCP-thresholding needed here, since SVD from Ssvd.als is not properly thresholded
    d=svd.xfill$d
    d=Eigen.Threshold(d,lambda,g,type.thresh)
    ### Forms matrix svd.xfill$v*d
    BD=UD(svd.xfill$v,d,m)
    U=svd.xfill$u
    ### Computes (U %*% t(BD))[i,j]= sum(u[i,]*v[j,]) for all pairs in irow,pcol
    xhat.Omega=suvC(U,BD,irow,pcol)
    ### Forming Sparse Part of a Sparse plus Low-Rank Matrix
    xres@x=x@x-xhat.Omega
    ### Forming object of type Sparse plus Low-Rank. Sparse: xres; Low-Rank: a=U, b=BD=VD    
    xfill=splr(xres,U,BD)
    ### Computes rank-restricted (soft-thresholded) SVD of xfill with warm-start svd.xfill and default lambda=0    
    svd.xfill=Ssvd.als(xfill,J,thresh,lambda=0,maxit,trace.it=FALSE,warm.start=svd.xfill)
    ### MCP-thresholding needed here, since SVD from Ssvd.als is not properly thresholded   
    ratio=Frob(svd.old$u[,seq(J)], d[seq(J)], svd.old$v[,seq(J)], svd.xfill$u[,seq(J)], Eigen.Threshold(svd.xfill$d,lambda,g,type.thresh)[seq(J)], svd.xfill$v[,seq(J)])   
    obj=(.5*sum(xres@x^2)+concave.penalty(d,lambda,g,type.thresh))/nz
    if(trace.it){
      cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
    }
  }
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
   ### Final cleanup of svd: do even if lambda=0
    A=xfill%*%svd.xfill$v
    Asvd=svd(A)
    U=Asvd$u
    V=svd.xfill$v %*% Asvd$v
   ###
   d=Eigen.Threshold(Asvd$d[seq(J)],lambda,g,type.thresh)
   J=min(sum(d>0)+1,J)
   svd.xfill=list(u=U[, seq(J)], d=d[seq(J)], v=V[,seq(J)])
   attributes(svd.xfill)=c(attributes(svd.xfill),list(lambda=lambda,call=this.call),biats)
   return(list(fit=svd.xfill, obj=obj, est.rank=length(which(d>1e-09))))
}