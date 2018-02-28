library(BoolNet)

net <- loadSBML(file = "regulatoryNetworkGMPModel.sbml")
att <- getAttractors(net,type = "synchronous",method="sat.exhaustive")

plotAttractors(att,onColor = "#33FF4C",offColor = "#FF6433")


