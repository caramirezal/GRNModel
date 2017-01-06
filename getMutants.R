



decimalToBinary<-function(value,positions){
  val = value
  output <- rep(0,positions)
  counter <- 1
  while(val>0){
    output[[counter]]<-val%%2
    val <- val %/% 2
    counter <- counter+1
  }
  output
}


getMutants<-function(net,phenotypes,fixed=c(0,1)){
  n<-length(net$genes)
  phenotypes.nb<-length(phenotypes)
  mutants<-matrix(0,n,phenotypes.nb)
  for ( i in 1:n ){
    attractors<-getAttractors(fixGenes(net,i,fixed),
                              type = "synchronous",
                              method = "sat.exhaustive")
    attractors.nb<-length(attractors$attractors)
    binaryAttractors<-lapply(attractors$attractors,
                             function(x) decimalToBinary(x[[1]][[1]],n) )
    for (k in 1:length(binaryAttractors)){
      names(binaryAttractors[[k]])<-net$genes
      attractor<-binaryAttractors[[k]]
      matchPatterns<-sapply(phenotypes,function(x) all(x==attractor[names(x)]) )
      for (l in 1:phenotypes.nb){
        if ( matchPatterns[l] == TRUE ) {
          mutants[i,l] <- mutants[i,l] + 1
        }
      }
    }
  }
  #mutants<-binarizeMatrix(mutants,0)
  rownames(mutants)<-net$genes
  colnames(mutants)<-names(phenotypes)
  return(mutants)
}

#mutants<-getMutants(net,phenotypes,fixed = 1)
#mutants