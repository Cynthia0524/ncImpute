Eigen.Threshold <- function(d,lambda,g,type.thresh=c("MC+","SOFT")){
  switch(type.thresh,
         "MC+"=MCP.thresh(d,lambda,g),   
         "SOFT"=Soft.thresh(d,lambda)
         )
}

MCP.thresh <- function(d,lambda,g){
   id1 = which(d<=lambda)
   id2 = which(d>lambda&d<=g*lambda)      
   d[id1] = 0
   d[id2] = sign(d[id2])*((abs(d[id2])-lambda)/(1-g^-1))
   return(d)
}

Soft.thresh <- function(d,lambda){
   d = pmax(d-lambda,0)
   return(d)
}

concave.penalty <- function(d,lambda,g,type.thresh=c("MC+","SOFT")){
  switch(type.thresh,
         "MC+"=MCP(d,lambda,g),
         "SOFT"=LASSO(d,lambda)         
         )
}

MCP <- function(d,lambda,g){
   id1 = which(d<=lambda*g)
   id2 = which(d>lambda*g)      
   d[id1] = lambda*d[id1] - (d[id1]^2)/(2*g)
   d[id2] = 0.5*g*lambda^2
   return(sum(d))
}

LASSO <- function(d,lambda){
   return(lambda*sum(d))
}