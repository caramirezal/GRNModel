This directory contains the following files:

i) modelFunctions.csv		contain a table with the Boolean functions of the model along with further detailed information about the experiments that supports the particular election of the function.

ii) modelInteractions.csv	contain a data base of the interactions already present in the regulatory network model.

iii) removedInteractions.csv 	comprise a table with the results of the simulations of the regulatory network model with systematic deletion of interactions (one by one).

iv) steadyStatesTransitions.csv	contain a table with transitions between stationary states derived by single transient node perturbations.

v) regulatoryNetworkGMPModel.txt	contain the GMP regulatory network full model in BoolNet format.

vi) regulatoryNetworkGMPModelDirectInteractions.txt	contain the GMP regulatory network first version model of only direct interactions in BoolNet format.

vii) regulatoryNetworkGMPModel.sbml	contain the GMP regulatory network full model in sbml standard format.

viii) supplementary-information.pdf	comprise the supplementary information of results mentioned in the main text.


Here, we describe the information given in each table.


## modelFunctions.xlsx file contains a table with  the following columns:

 - node: the variable node name as appears in the model implementation in the regulatoryNetworkGMPModel.txt file.

 - factors: refers to the Boolean function for the corresponding node. Where "&", "|", ando "!" represent "AND", "OR", and "NOT" Boolean operators. 
Boolean subexpressions of direct interactions are shown in black, whilst indirect and/or proposed interactions are shown in green.

 - Boolean expression details: give a description of the experimental findings which leaves to the definition of each model Boolean function.


## modelInteractionsDB.xlsx comprise a database with the following column variables:

 - node: the variable node name as appears in the model implementation in the regulatoryNetworkGMPModel.txt file.

 - regulator: a variable name which appears in the Boolean function of a node.

 - type: a cathegorical value that indicates if the interaction is direct or indirect or proposed.

 - behavior: a cathegorical value which indicates if the regulator variable acts as a positive, negative or dual input to the node.

 - Molecular.technique: contain further details about the molecular methods to experimentally infer the interaction.

 - Specie: A cathegorical label that indicates if the interaction was inferred in murine, human or both experimental model species.

 - Cell.line: string values with the names of the cell culture used in the experiments to infer the interaction.

 - PMID: numerical ID from PubMed for the reference that supports the interaction.


## removedInteractions.xlsx file comprise a table with the following columns:

 - netInteraction: a string showing the interactions which is deleted from the model (regulator -> node). 

 - missingAttractors: the number of missing wild type steady states which are not longer solutions in the new model with the removed interaction.

 - missingGMP: the number of wild type GMP phenotypes that does not appear as steady states in the perturbed network.

 - type: cathegorical value that indicates if the interactions are direct, indirect or proposed.


## steadyStatesTranstions.xlsx file contains a table with the following columns:

 - Attractor: a string name for the steady states of the model, labeled with GMP lineage names according to the pressence of the molecular signatures 
given in Supplementary file 5. The 21 steady states shown in Figure 1 are present in this column with indexes to distinguish between steady states
in the same class. GMP derived phenotypes are tagged as follows: Lineage negative (Lne), neutrophils (neu), basophils (bas), eosinophils (eos),
and mast cells (mas).

 - newAttractor: the name of the steady state which was reached after the a node was perturbed. 

 - perturbation: a node name indicating which node was perturbed.


All the scripts used to derive the results shown in this work are available in the https://github.com/caramirezal/RegulatoryNetworkGMPModel repository along with further documentation. Additionally, tables and figures are reproduced in the RegulatoryNetworkGMPModel.ipynb and regulatoryNetworkGMPModelSupplementaryInfo.ipynb in the same repository.










