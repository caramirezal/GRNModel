# takes a list of nodes with values and performs fixed environment simulations
fixEnvironment <- function(net,fixValues,attractors,phenotypes,...){
        fixNodes <- sapply(names(fixValues),function(x) which(x==net$genes))
        att.wt <-tagAttractorsMatrix(net=net,attractors=attractors,patternsList = phenotypes)
        att.wt <- unique(colnames(att.wt))
        attlist <- list()
        for (i in 1:length(taggedAttractors[1,])){
                attlist[[i]] <- taggedAttractors[,i]
                attlist[[i]][fixNodes] <- fixValues 
        }
        att.fixEnv <- getAttractors( fixGenes(net,fixIndices = fixNodes,
                                              values = fixValues ) , 
                                     startStates = attlist, ...)
        att.fixEnv <- tagAttractorsMatrix(net = net,attractors = att.fixEnv,patternsList = phenotypes)
        att.fixEnv <- colnames(att.fixEnv)
        indexes <- sapply(att.fixEnv,function(x) which(x == att.wt))
        indexes <- unlist(indexes[sapply(indexes,length)>0])
        indexes <- unique(indexes)
        res <- rep(0,length(att.wt))
        names(res) <- att.wt
        res[indexes] <- 1
        res
}

newAtt <- fixEnvironment(net,fixValues,attractors,phenotypes)
newAtt

fixedEnvironments <- list("ProNeu"=c("cebpa"=1,"pu1"=0),
                        "ProMono"=c("cebpa"=0,"pu1"=1),
                        "Promast"=c("il3ra"=1,"ckit"=1),
                        "ProEo"=c("cebpa"=1,"gata1"=1,"gata2"=1),
                        "Probas"=c("cebpa"=1,"gata1"=1,"pu1"=0))

simEnvironments <- function(net,fixedEnvironments,attractors,phenotypes,envNames="default",...){
        
        res <- list()
        
        for (i in 1:length(fixedEnvironments)) {
                
                sim <- fixEnvironment(net,fixValues = fixedEnvironments[[i]],
                                      attractors = attractors,phenotypes = phenotypes)
                res[[i]] <- sim
                
        }
        
        res <- as.data.frame(res)
        
        colnames(res) <- names(fixedEnvironments)
        
        plotAttractors.a(t(as.matrix(res)))
        
        res
        
        
}

res <- simEnvironments(net,fixedEnvironments,attractors,phenotypes)


