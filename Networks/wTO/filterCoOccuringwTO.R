## Removes associations between taxa that do not occur > X% together (default 50).
## (i.e. 50% of the times that at least one of them occur)
## Accepts object from e.g. wTO w Results (The 2 first columns being node names)
## and an OTU table w taxa as rows

filterCooccuring = function(edges, taxa_table, minOverlap=.5){
  tt = taxa_table
  edges = as.data.frame(edges)
  toKeep = rep(TRUE, dim(edges)[1])
  for (ei in c(1:dim(edges)[1])){
    abundances = tt[row.names(tt) %in% c(as.character(edges[ei,1]),as.character(edges[ei,2])),]
    both = 0
    one = 0
    for (ai in c(1:dim(abundances)[2])){
      taxa_occur = abundances[,ai]
      
      if (taxa_occur[1] > 0 && taxa_occur[2]>0) both = both +1
      else {
        if (xor(taxa_occur[1] > 0, taxa_occur[2]>0)) {
          one = one + 1
        }
      }
    }
    if (one > 0 && (both / (both + one)) < minOverlap) toKeep[ei] = FALSE
  }
  return(edges[toKeep,])
}
