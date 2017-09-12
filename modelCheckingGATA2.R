
## model checking

## dependencies
library(BoolNet)
library(dplyr)
source("usefulFunctions.R")

## loading interactions database
interactions.db <- read.csv("additionalInfo/SupInfo/modelInteractions.csv")

## loading network model
net <- loadNetwork("regulatoryNetworkGMPModel.txt")
## consider regulators only
regulators <- findRegulators(net)

## formatting interactions to Griffin input
## filtering interactions
interactions.r <- subset(interactions.db,node%in%net$genes[regulators])
interactions.r <- interactions.r[,c("regulator","node","behavior")]
interactions.s <- character(length(interactions.r$behavior))
for ( i in 1:length(interactions.r$behavior) ) {
         if ( interactions.r$behavior[i] == "positive") {
                 interactions.s[i] = "->"
         }
        if ( interactions.r$behavior[i] == "negative") {
                interactions.s[i] = "-|"
        }
}
interactions.r$"symbol" <- interactions.s
interactions.r <- mutate(interactions.r, interaction=paste(regulator,symbol,node,sep = ""))
## definition of network topology
topology <- paste(interactions.r$interaction,collapse = ",")


