# MultiLayerLinkage
Code for simulations and simulated data analysis from the paper A Bayesian Multi-Layered Record Linkage Procedure to Analyze Functional Status of Medicare Patients with Traumatic Brain Injury. 

DiscreteSimulations folder contain the code when assuming that the comparison function for the Income-Level variable is that the variables represent the same entry when the difference is smaller than 500, and 0 otherwise.

ContinuousSimulations folder contain the code when assuming that the comparison function for the Income-Level variable is the Euclidean Distance.

GraphGeneration recreate the heatmap that is provided in the paper


SimulatedRealData provides the algorithm that was run on the real dataset. Because the data cannot be shared we implement it on a simulated dataset. 
