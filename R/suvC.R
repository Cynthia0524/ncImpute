suvC <-
function(u,v,irow,pcol){
###  computes (u%*%t(v))[i,j]= sum(u[i,]*v[j,]) for all pairs in irow,pcol
  dd=dim(u)
  nnrow=as.integer(dd[1])
  nncol=as.integer(nrow(v))
  nrank=dd[2]
  storage.mode(u)="double"
  storage.mode(v)="double"
  storage.mode(irow)="integer"
  storage.mode(pcol)="integer"
  nomega=as.integer(length(irow))
  .Fortran("suvC",
           nnrow,nncol,nrank,u,v,irow,pcol,nomega,
           r=double(nomega),
           PACKAGE="ncImpute"
           )$r
}
