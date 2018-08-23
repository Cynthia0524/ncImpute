#### Matrix-vector products of the form (S + ab') %*% w and (S + ab')' %*% w for PROAPACK SVD ####
f <- function(w){ # Convert Matrix object back to vector
 part1 = S %*% w
 part2 = t(b) %*% w
 part2 = a %*% part2
 as.numeric(part1 + part2)

} 
tf <- function(w){ # Convert Matrix object back to vector
 part1 = t(S) %*% w
 part2 = t(a) %*% w
 part2 = b %*% part2
 as.numeric(part1 + part2)
} 
