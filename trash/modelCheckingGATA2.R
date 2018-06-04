
## model checking

## dependencies
library(BoolNet)
library(dplyr)
source("usefulFunctions.R")
source("extractAttractors.R")
source("tagAttractorsMatrix.R")

## loading interactions database
interactions.db <- read.csv("additionalInfo/SupInfo/modelInteractions.csv")

## loading network model
net <- loadNetwork("regulatoryNetworkGMPModel.txt")
## consider regulators only
regulators <- findRegulators(net)
## formatting interactions to Griffin input
## filtering interactions
interactions.r <- subset(interactions.db,node%in%net$genes[regulators])
interactions.r <- interactions.r[,c("regulator","node","behavior","type")]

formatInteraction <- function(interactionsMatrix){
        interactions.s <- character(length(interactionsMatrix$behavior))
        for ( i in 1:length(interactionsMatrix$behavior) ) {
                if ( interactionsMatrix$behavior[i] == "positive") {
                        interactions.s[i] = "->"
                }
                if ( interactionsMatrix$behavior[i] == "negative") {
                        interactions.s[i] = "-|"
                }
        }
        interactionsMatrix$"symbol" <- interactions.s
        interactionsMatrix <- mutate(interactionsMatrix, interaction=paste(regulator,symbol,node,sep = ""))
        ## definition of network topology
        topology <- paste(interactionsMatrix$interaction,collapse = ",")
        topology
}


## known interactions
known <- subset(interactions.r,type=="direct")
known <- formatInteraction(known)
known <- paste("known = {",known,"}",sep="")

## hypothetical interactions
hypothesis <- subset(interactions.r,type!="direct")
hypothesis <- formatInteraction(hypothesis)
hypothesis <- paste("hypothetical = {",hypothesis,"}",sep="")

## genes input to griffin 
genesInput <- paste(net$genes[regulators],collapse = ",")
genesInput <- paste("genes = {",genesInput,"}",sep = "")
genesInput

## options
#options
options <- "allow.ambiguity = false
allow.additional.states = true 
allow.additional.cycles = true 
allow.hypotheses = true 
block.steady.a.posteriori = false 
divide.query.by.topology = false"

## getting attractors matrix
attractors <- getAttractors(net,
                            type="synchronous",
                            method = "sat.exhaustive")
attractors.m <- extractAttractors(net,attractors)
attractors.m <- attractors.m[regulators,]


formatAttractors <- function(attractorsMatrix){
        attractorsInput <- sapply(1:length(attractorsMatrix[1,]),
                                  function(x) 
                                          paste(attractorsMatrix[,x],
                                                collapse=""))
        attractorsInput <- paste(attractorsInput,collapse=",")
        attractorsInput
}


## wild type attractors
attractors.wt <- formatAttractors(attractors.m)
attractors.wt <- paste("fixed-points() = {",
                         attractors.wt,"}",
                         sep ="")


## gata2 mutants
# definition of the phenotype patterns 
monocytes <- c("mcsfr"=1,"pu1"=1)
neutrophils <- c("lf"=1, "cebpa"=1)
eosinophils <- c("fceRIa"=1,"gata1"=1,"cebpa"=1,"ccr3"=1)
basophils <- c("cebpa"=1,"gata2"=1,"runx1"=1,"ccr3"=0)
mast <- c("mitf"=1,"ckit"=1,"cebpa"=0)
Lne <- c("mcsfr"=0,"lf"=0,"mbp"=0,"ckit"=0,"fceRIa"=0,"ccr3"=0)
phenotypes <- list("monocytes"=monocytes,"neutrophils"=neutrophils,"eosinophils"=eosinophils,
                 "mast cells"=mast,"Lne"=Lne,"basophils"=basophils)
taggedAttractors <- tagAttractorsMatrix(net,
                                        attractors,
                                        patternsList = phenotypes)
taggedAttractors <- taggedAttractors[regulators,]

## subsetting non mast cells attractors
nonMastCells <- which(colnames(taggedAttractors)!="mas")
expectedGata2 <- taggedAttractors[,nonMastCells]
expectedGata2 <- formatAttractors(expectedGata2)
expectedGata2 <- paste("fixed-points(~gata2) = {",
                       expectedGata2,"}",
                       sep ="")

## cjun mutants
expectedcjun <- formatAttractors(attractors.m)
expectedcjun <- paste("fixed-points(~cjun) = {",
                       expectedcjun,"}",
                       sep ="")

## writting griffing input
writeLines(text =  c(genesInput,
                     known,
                     hypothesis,
                     attractors.wt,
                     expectedGata2,
                     expectedcjun,
                     options),
           con = "model_checking_gata2.grf",sep = "\n")

