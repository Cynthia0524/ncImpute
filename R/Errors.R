Test.Error <- function(U,d,V,Z.est,imiss){
   Z = U%*%(d*t(V))
   num = sum((Z-Z.est)[imiss]^2)
   den = sum(Z[imiss]^2)  
   return(num/den)
}

Training.Error <- function(U,d,V,E,object,imiss){
   Z = U%*%(d*t(V)) + E
   n = nrow(Z); p = ncol(Z); np=n*p
   obs = setdiff(seq(np),imiss)
   # Obtains observed indices in (i,j) format
   i=row(Z)[obs]
   j=col(Z)[obs]  
   # Computing the new guess only at the observed entries Omega
   Z.est=impute(object,i,j)  
   num = sum((Z[obs]-Z.est)^2)
   den = sum(Z[obs]^2)
   return(num/den)
}

Round.Ratings <- function(Z){
  Z[Z>5]=5
  Z[Z<1]=1
  return(Z)
}

RMSE <- function(Z,True){
   n=length(Z)
   sum.err=(1/n)*(sum((Z-True)^2))
   return(sqrt(sum.err))
}