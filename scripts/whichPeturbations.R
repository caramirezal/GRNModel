
# return a list of deterministic perturbations and the transition derived by them

whichPerturbations<-function(net,attractors,iterations=1,file="yes",attractorsNames){
  # 0. define a list of results
  res<-list()
  # 1. traverse attractors. Take an attractor i, define a new i-th list in the result list
  #      defined in 0.
  if (file=="yes"){
    sink("whichPerturbations.txt")
    cat("Attractor,newAttractor,perturbation \n")
  }
  numberOfFixedPoints<-getNumberOfFixedPoints(attractors)
  for (i in 1:numberOfFixedPoints){
    transitions<-list()
    binaryAttractor<-decimalToBinary(attractors$attractors[[i]][[1]][[1]],length(net$genes))
    # 2. For every node, give a node perturbation to an attractor j
    for (j in 1:length(net$genes)){
      pulsedState<-flipNode(binaryAttractor,j)               # gives a single transient perturbations
      for (n in 1:iterations) {
        newAttractor<-getAttractors(net,method="chosen",startStates=list(pulsedState),type = "asynchronous")
        newAttractor<-decimalToBinary(newAttractor$attractors[[1]][[1]][[1]],length(net$genes) )
        newAttractor<-findAttractor(net,newAttractor,attractors,type = "asynchronous",
                                    numberOfFixedPoints = numberOfFixedPoints)
        # return attractor indexes
        # 3. Verify a transition to another attractor.
        if (i != newAttractor) {
          if ( 0 == length(transitions) ){
            transitions[[1]]<-c(j,newAttractor)
            cat(attractorsNames[i],",",attractorsNames[ transitions[[1]][2] ],",",net$genes[ transitions[[1]][1] ],"\n")
          }
          else { 
            counter<-1
            found<-FALSE
            while ( ( counter <= length(transitions) ) & ( found == FALSE ) ){
              found<- all( c(j,newAttractor) == transitions[[counter]] )
              counter<-counter+1
            }
            if ( found == FALSE ) {
              transitions[[length(transitions)+1]]<-c(j,newAttractor)
              cat(attractorsNames[i],",",attractorsNames[ transitions[[length(transitions)]][2] ],",",net$genes[ transitions[[length(transitions)]][1] ],"\n")
            }
          }
        }
      }
    }
    res[[i]]<-transitions
  }
  if (file=="yes"){
    sink()
  } else {
   return(res)
  }  
}



