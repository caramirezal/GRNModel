library(BoolNet)
library(reshape2)
library(ggplot2)


# Loading local scripts
source("tagAttractorsMatrix.R")
source("plotAttractors.a.R")
source("usefulFunctions.R")
source("getMutants.R")

## getting synchronous attractors
net <- loadNetwork("regulatoryNetworkGMPModel.txt")
attractors <- getAttractors(net,type = "synchronous",
                            method = "sat.exhaustive")

# definition of the phenotype patterns 
monocytes<-c("mcsfr"=1,"pu1"=1)
neutrophils<-c("lf"=1, "cebpa"=1)
eosinophils<-c("fceRIa"=1,"gata1"=1,"cebpa"=1,"ccr3"=1)
basophils<-c("cebpa"=1,"gata2"=1,"runx1"=1,"ccr3"=0)
mast<-c("mitf"=1,"ckit"=1,"cebpa"=0)
Lne<-c("mcsfr"=0,"lf"=0,"mbp"=0,"ckit"=0,"fceRIa"=0,"ccr3"=0)
phenotypes<-list("monocytes"=monocytes,"neutrophils"=neutrophils,"eosinophils"=eosinophils,
                 "mast cells"=mast,"Lne"=Lne,"basophils"=basophils)

## tagging attractors
fixedPoints <- tagAttractorsMatrix(net = net,
                                   attractors = attractors,
                                   patternsList = phenotypes)

#############################################################################

## calculating gain of function mutants
mutantsNull<-getMutants(net, phenotypes,fixed = 0)

# to fix the value of a non-regulator node does not change dynamics, so they are discarded.
# Only regulators are taken into account.
regulators<-findRegulators(net)

## simplifying the matrix. 0 (1) accounts for the absence (presence)
## of the attractor phenotype
mutantsNull.r<-mutantsNull[regulators,]
mutantsNull.r<-binarizeMatrix(mutantsNull.r,1)

#write.csv(mutantsNull.r,file = "nullMutantOutput.csv")
plotAttractors.a(mutantsNull.r)

############################################################################

## calculating gain of function mutants
mutantsgof<-getMutants(net, phenotypes,fixed = 1)

## simplifying the matrix. 0 (1) accounts for the absence (presence)
## of the attractor phenotype
mutantsgof.r<-mutantsgof[regulators,]
mutantsgof.r<-binarizeMatrix(mutantsgof.r,1)

plotAttractors.a(mutantsgof.r)
