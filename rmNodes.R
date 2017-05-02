

## remove the string name of a node (nodeName) 
## for a single boolean expression (boolExp)
rmNode <- function(nodeName, boolExp){
        
        boolExp <- gsub(" ","",boolExp)
        
        ## filtering node with operators on the left
        operators <- c("\\&!","!","\\&","\\|")
        for (i in 1:length(operators)) {
                nodeExpression <- paste(operators[i],nodeName,sep = "")
                boolExp <- gsub(nodeExpression,"",boolExp)
        }
        
        ## filtering node with operators on the right
        for (i in 3:length(operators)) {
                nodeExpression <- paste(nodeName,operators[i],sep = "")
                boolExp <- gsub(nodeExpression,"",boolExp)
        }
        
        boolExp <- gsub(paste("\\(",nodeName,"\\)",sep = ""),"",boolExp)
        boolExp <- gsub(paste("^",nodeName,sep = ""),"",boolExp)
        
        boolExp <- gsub("\\(\\)","",boolExp)
        boolExp <- gsub("^[\\||\\&]","",boolExp)
        boolExp <- gsub("[\\||\\&]$","",boolExp)
        
        if ( length(grep(paste("[[:punct:]]",nodeName,"[[:punct:]]",sep = ""),boolExp)) > 0) {
                
                stop("Something went wrong. The node was not removed: ",boolExp)
                
        }
        
        boolExp
        
}


## test for the rmNode() function
rmNodeTest <- function(net){
        
        nodeExpressions <- sapply(net$interactions,function(x) x$expression) 
        
        for (i in nodeExpressions){
                cat("Node Expression: ",i,"\n")
                
                for (j in net$genes) {
                        
                        if (length(grep(j,i))>0) {
                                cat("Removed node: ",j,"\n")
                                cat("Mutated expression: ", rmNode(j,i),"\n" )
                        } 
                        
                }
                
                cat("\n")
        }
        
}

## This expression parse a node name
##grep("^a[[:punct:]| ]|[[:punct:]| ]a[[:punct:]| ]|[ |[:punct:]]a$",c("a ","aaa"," a ","&a&", " a "))

library(BoolNet)
library(dplyr)

source("tagAttractors.R")
source("usefulFunctions.R")
source("tagAttractorsMatrix.R")
source("getTransitionMatrix.R")
netDirectInt <- loadNetwork("regulatoryNetworkGMPModelDirectInteractions.txt")
netFull <- loadNetwork("regulatoryNetworkGMPModel.txt")

rmInteractions <- function(file,attractors,phenotypes){
        
        net <- loadNetwork(file)
        
        res <- list("affectedNode"=c(),"removedNode"=c(),"percentageConservedAttractors"=c())
        
        netLines <- readLines(file)
        
        taggedAttractors <- tagAttractorsMatrix(net = net, attractors = attractors, patternsList = phenotypes)
        
        attList <- lapply(1:length(taggedAttractors[1,]), function(x) taggedAttractors[,x] )
        
        taggedAttractors <- colnames(taggedAttractors)
        
        booleanExp <- sapply(net$interactions,function(x) x$expression)
        
        for (i in 1:length(booleanExp)) {
                
                regulators <- net$genes[ net$interactions[[i]]$input ]
                
                for (j in regulators) {
                        
                        newExp <- rmNode(j,booleanExp[i])
                        
                        if (newExp == "") {
                                
                                newExp <- net$genes[i]
                                
                        }
                        
                        newExp <- paste(net$genes[i],", ",newExp,sep = "")
                        
                        newNetStrings <- netLines
                        
                        newNetStrings[i+1] <- newExp
                        
                        write(newNetStrings,sep = "\n",file = "newNet.tmpl")
                        
                        newNet <- loadNetwork("newNet.tmpl")

                        
                        newAtt <- getAttractors(newNet,startStates = attList,
                                                type = "asynchronous")
                        
                        newTagAtt <- tagAttractorsMatrix(newNet,newAtt,patternsList = phenotypes)
                        
                        newTagAtt <- colnames(newTagAtt)
                        
                        newTagAtt <- newTagAtt[newTagAtt %in% taggedAttractors]
                        
                        perConAtt <- length(newTagAtt) / length(taggedAttractors) 
                        
                        res$"affectedNode" <- c(res$"affectedNode",net$genes[i])
                        
                        res$"removedNode" <- c(res$"removedNode",j)
                        
                        res$"percentageConservedAttractors" <- c(res$"percentageConservedAttractors",perConAtt)
                        
                
                }
                
        }
        
        file.remove("newNet.tmpl")
        
        res <- as.data.frame(res)
        
        arrange(res, percentageConservedAttractors)
        
}

# loading the full BRN model
net<-loadNetwork("regulatoryNetworkGMPModel.txt")

# definition of the phenotype patterns. In the full model
# CCR3, CEBPb, and IL3Ra signalling is added.
monocytes<-c("mcsfr"=1,"pu1"=1)
neutrophils<-c("lf"=1, "cebpa"=1)
eosinophils<-c("fceRIa"=1,"gata1"=1,"cebpa"=1,"ccr3"=1)
basophils<-c("cebpa"=1,"gata2"=1,"runx1"=1,"ccr3"=0)
mast<-c("mitf"=1,"ckit"=1,"cebpa"=0)
Lne<-c("mcsfr"=0,"lf"=0,"mbp"=0,"ckit"=0,"fceRIa"=0,"ccr3"=0)
phenotypes<-list("monocytes"=monocytes,"neutrophils"=neutrophils,
                 "eosinophils"=eosinophils,"mast cells"=mast,
                 "Lne"=Lne,"basophils"=basophils)

# calculate attractors
#attractors<-getAttractors(net, method = "sat.exhaustive")

#rs <- rmInteractions(file = "regulatoryNetworkGMPModel.txt",attractors = attractors,phenotypes = phenotypes )

