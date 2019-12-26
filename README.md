# EntroRank
An entropy based method for identifying cancer driver genes
Run the EntroRank should under the environment of R >=2.7
1. load the influencegraph, somatic mutation matrix(prostate,lung or breast) and the compartment information from the offered data file.
2. invoke the main_function(influenceGraph)
The final result return a list of ranking driver genes and its corresponding score.
