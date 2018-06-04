# Useful functions called frequentle by various addtional functions, so they are 
# convenientley defined in this script.

#######################################################################################

# convert a decimal value to binary (vector)
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


#######################################################################################

# transform a logical vector to a integer number
bitsToInt<-function(x) {
  packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
}

#######################################################################################
# From the total of attractors
# obtain the number of fixed points
getNumberOfFixedPoints<-function(attractors){ 
  numberOfFixedPoints<-0
  for ( i in 1:length(attractors$attractors) ) {
    periodSize<-length(attractors[2]$attractors[[i]][[1]])
    if ( 1 == periodSize ) {
      numberOfFixedPoints<-numberOfFixedPoints + 1  
    }
  }
  return(numberOfFixedPoints)
}


#######################################################################################

# convert a matrix to one which have binary values
# if the entry is > 0 convert the value to 1
binarizeMatrix<-function(gMatrix,umbral){
  for (i in 1:length(gMatrix[,1])) {
    for (j in 1:length(gMatrix[1,])){
      if ( gMatrix[i,j] > umbral ) { gMatrix[i,j] <- 1 }
      else { gMatrix[i,j] <- 0 }
    }
  }
  return(gMatrix)  
}

###################################################################################################
# Create a random state
randomState <- function(net){
  randomState<-vector(mode="numeric",length(net$genes))
  for (i in 1:length(net$genes)) { 
    randomState[i] <- sample(0:1,1)
  }
  return(randomState)
}


#####################################################################################################

# find all nodes that appears as inputs in the net
findRegulators<-function(net){
  regulators<-c()
  for (i in 1:length(net$genes)){
    nodeIndex<-i
    #cat(i,"\n")
    for (j in 1:length(net$genes)){
      # cat(net$genes[j],"\n")
      if ( nodeIndex %in% net$interactions[[j]][[1]]) {
        # cat(nodeIndex,net$interactions[[j]][[1]],"is a regulator \n")
        if  ( ! ( nodeIndex %in% regulators) ) {
          regulators<-c(regulators,nodeIndex)
        }
      }
    }
  }
  return(regulators)
}
