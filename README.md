# Regulatory Network model of Granulocyte-Monocyte derived cells

This github project provides all the scripts used to construct the images of the paper named "Phenotypic stability and plasticity in GMP-derived cells as determined by their underlying regulatory network" (work in progress). It also provides documentation.

RegulatoryNetworkGMPModel.ipynb shows how to construct the images of the main text.

RegulatoryNetworkGMPModelSupplementaryInfo.ipynb showS how to construct the images of the supplementary information.

CellFateMap.R Performs all possible bit flips perturbations to model attractors and retrieves a matrix of transitions between stationary fixed points.

getMutants.R sistematically performs all gain (or lost) of function mutants.

getTransitionMatrix.R simulates noise driven transitions between attractors.

plotAttractors.a.R plots attractors.

simplifyCellFateMap.R simplifies the matrix of transitions between attractors merging attractors that belongs to the same class or group, that is, if they have related phenotypes defined as sharing molecular signatures.

tagAttractors.R sistematically tags attractors according to molecular signatures that relates to phenotypes.

fixedEnvironments.R performs simulations on fixed environments. 

rmNodes.R systematically removes every interaction in each node Boolean expression and record whether the attractors of the Wild Type network are presented in the new attractors. It returns a table with the fraction of attractors that match these criteria.

