### Computes the matrix product U*D, where D is a diagonal matrix given in vector format
UD=function(U,D,n=nrow(U)){
  U*outer(rep(1,n),D,"*")
}
