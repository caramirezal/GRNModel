

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


library(BoolNet)

netDirectInt <- loadNetwork("regulatoryNetworkGMPModelDirectInteractions.txt")
netFull <- loadNetwork("regulatoryNetworkGMPModel.txt")

n = 1
netStrings <- readLines("regulatoryNetworkGMPModel.txt")
netStrings[n+1] <- rmNode(netFull$genes[1],netStrings[n+1]) 
write(netStrings,"newNet.txt",sep = "\n")


readLines("newNet.txt")
## This expression parse a node name
##grep("^a[[:punct:]| ]|[[:punct:]| ]a[[:punct:]| ]|[ |[:punct:]]a$",c("a ","aaa"," a ","&a&", " a "))

