# Loading dependencies
library(BoolNet)
library(reshape)
library(ggplot2)
library(igraph)
library(dplyr)


# Loading local scripts (found at https://github.com/caramirezal/RegulatoryNetworkGMPModel)
source("cellFateMap.R")
source("getTransitionMatrix.R")
source("simplifyCellFateMap.R")
source("tagAttractorsMatrix.R")
source("plotAttractors.a.R")
source("usefulFunctions.R")
source("getMutants.R")
source("whichPeturbations.R")
source("fixedEnvironments.R")
source("rmNodes.R")

# Loading the regulatory network model of 
# well grounded direct interactions
net<-loadNetwork("regulatoryNetworkGMPModelDirectInteractions.txt")

# definition of the phenotype patterns. Note that in the full regulatory network model
# the marker ccr3 is added to distinguish between basophils and eosinophils phenotypes.
monocytes<-c("mcsfr"=1,"pu1"=1)
neutrophils<-c("lf"=1, "cebpa"=1)
eosinophils<-c("fceRIa"=1,"gata1"=1,"cebpa"=1,"mbp"=1)
basophils<-c("cebpa"=1,"gata2"=1,"runx1"=1,"mbp"=1)
mast<-c("mitf"=1,"ckit"=1,"cebpa"=0)
Lne<-c("mcsfr"=0,"lf"=0,"mbp"=0,"ckit"=0,"fceRIa"=0)
phenotypes<-list("monocytes"=monocytes,"neutrophils"=neutrophils,
                 "eosinophils"=eosinophils,"mast cells"=mast,
                 "Lne"=Lne,"basophils"=basophils)



# calculating attractors of the BRN of well grounded direct interactions
attractors<-getAttractors(net, method = "sat.exhaustive")

# tag attractors according to the defintions given at the beginning
taggedAttractors<-tagAttractorsMatrix(net,attractors,phenotypes)

# Plotting the attractors
pdf("figures/attractorsFirstModel.pdf")
heatmap(taggedAttractors,Rowv = NA,Colv=NA, col=c("#FF6433","#33FF4C"))
dev.off()

##########################################################################

# Loading dependencies
library(BoolNet)
library(reshape)
library(ggplot2)

# Loading local scripts
source("tagAttractorsMatrix.R")
source("plotAttractors.a.R")
source("usefulFunctions.R")

# Definition of the full regulatory network model 
net<-loadNetwork("regulatoryNetworkGMPModel.txt")

# definition of the phenotype patterns 
monocytes<-c("mcsfr"=1,"pu1"=1)
neutrophils<-c("lf"=1, "cebpa"=1)
eosinophils<-c("fceRIa"=1,"gata1"=1,"cebpa"=1,"ccr3"=1)
basophils<-c("cebpa"=1,"gata2"=1,"runx1"=1,"ccr3"=0)
mast<-c("mitf"=1,"ckit"=1,"cebpa"=0)
Lne<-c("mcsfr"=0,"lf"=0,"mbp"=0,"ckit"=0,"fceRIa"=0,"ccr3"=0)
phenotypes<-list("monocytes"=monocytes,"neutrophils"=neutrophils,"eosinophils"=eosinophils,
                 "mast cells"=mast,"Lne"=Lne,"basophils"=basophils)

# calculate attractors
attractors<-getAttractors(net, method = "sat.exhaustive")

# tag attractors according to the defintions given at the beginning
taggedAttractors<-tagAttractorsMatrix(net,attractors,phenotypes)

attractors.nm<-colnames(taggedAttractors)
# A "neu/GMP" attractor has a neutrophil pattern, hence it is renamed "neu"
attractors.nm[which("neu/GMP"==attractors.nm)]<-"neu"
# A "GMP/Lne" attractor pattern is renamed as GMP, because it has a progenitor molecular phenotype
#attractors.nm[which("GMP/Lne"==attractors.nm)]<-"GMP"
colnames(taggedAttractors)<-attractors.nm

# Giving order to attractors
Lne<-c(1,11,12,15,17)
neu<-c(2,16)
eos<-c(7,8,9,10)
mon<-c(18,19)
mas<-c(5,6,13,14,20,21)
bas<-c(3,4)
order<-c(Lne,neu,eos,mon,mas,bas)

# Plotting attractors
plotAttractors.a(taggedAttractors,colOrder=TRUE,colIndexes=order)

#########################################################################
## fixed environments simulation 
fixedEnvironments <- list("Pro Neu"=c("cebpa"=1,"pu1"=0,"gcsfr"=1),
                          "Pro Mon"=c("cebpa"=0,"pu1"=1,"mcsfr"=1),
                          "Pro Mas"=c("il3ra"=1,"ckit"=1),
                          "Pro Eos"=c("cebpa"=1,"gata1"=1,"fceRIa"=1),
                          "Pro Bas"=c("cebpa"=1,"gata2"=1,"runx1"=1,"il3ra"=1),
                          "Null TFs"=c("cebpa"=0,"pu1"=0,"gata1"=0,"mitf"=0),
                          "Cytokines"=c("mcsfr"=1,"gcsfr"=1,"gmcsfr"=1,"il3ra"=1,"ckit"=1,"fceRIa"=1),
                          "Null cytokines"=c("mcsfr"=0,"gcsfr"=0,"gmcsfr"=0,"il3ra"=0,"ckit"=0,"fceRIa"=0))



fixed_env <- simEnvironments(net,fixedEnvironments,attractors,phenotypes)

heatmap(as.matrix(fixed_env),Rowv = NA,Colv = NA,col=c("#FF6433","#33FF4C"))
