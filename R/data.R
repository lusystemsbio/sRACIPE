#' @title A toggle switch circuit for demonstrations
#' @description This data contains the topology of a 
#' circuit with two genes A and B both of which inhibit each
#' other and have self activations. For further details,
#' see Kohar, V. and Lu, M., Role of noise and parametric 
#' variation in the dynamics of gene regulatory circuits, 
#' npj Systems Biology and Applications 4, Article number: 40 (2018)
#' @format A data frame with 4 rows and 3 variables
"demoCircuit"

#' @title Five coupled toggle switches
#' @description This data contains the topology of a 
#' circuit with ten genes A1...A5 and B1...B5. 
#' Genes with different alphabet inhibit each other 
#' whereas genes from different groups activate the other.
#' For further details,
#' see Kohar, V. and Lu, M., Role of noise and parametric 
#' variation in the dynamics of gene regulatory circuits, 
#' npj Systems Biology and Applications 4, Article number: 40 (2018)
#' @format A data frame with 18 rows and 3 variables
"CoupledToggleSwitchSA"

#' @title A circuit for epithelial to mesenchymal transition
#' @description This data contains the topology of a 
#' circuit with sixteen genes invovled in EMT. 
#' For further details,
#' see Kohar, V. and Lu, M., Role of noise and parametric 
#' variation in the dynamics of gene regulatory circuits, 
#' npj Systems Biology and Applications 4, Article number: 40 (2018)
#' @format A data frame with 59 rows and 3 variables
"EMT1"

#' @title A circuit for epithelial to mesenchymal transition 
#' including microRNAs
#' @description This data contains the topology of a 
#' circuit with twenty two nodes including
#' micro RNAs invovled in EMT. 
#' For further details,
#' see Huang et al.,Interrogating the topological robustness 
#' of gene regulatory circuits by randomization, 
#' PLoS computational biology 13 (3), e1005456
#' @format A data frame with 82 rows and 3 variables
"EMT2"

#' @title Configuration Data
#' @description It contains simulation parameters 
#' like integration method
#' (stepper) and other lists or vectors like simParams, 
#' stochParams, hyperParams, options, thresholds etc. 
#' The list simParams contains values for parameters like the 
#' number of models (numModels), 
#' simulation time (simulationTime), step size for simulations 
#' (integrateStepSize), when to start recording the gene expressions 
#' (printStart), time interval between recordings (printInterval), number of 
#' initial conditions (nIC), output precision (outputPrecision), tolerance for
#' adaptive runge kutta method (rkTolerance), parametric variation (paramRange).
#' The list stochParams contains the parameters for stochastic simulations like
#' the number of noise levels to be simulated (nNoise), the ratio of subsequent
#' noise levels (noiseScalingFactor), maximum noise (initialNoise), whether to
#' use same noise for all genes or to scale it as per the median expression of
#' the genes (scaledNoise), ratio of shot noise to additive noise (shotNoise).
#' The list hyperParams contains the parameters like the minimum and maximum 
#' production and degration of the genes, fold change, hill coefficient etc.
#' The list options includes logical values like annealing (anneal), scaling of 
#' noise (scaledNoise), generation of new initial conditions (genIC), parameters
#' (genParams) and whether to integrate or not (integrate). The user
#' modifiable simulation options can be specified as other arguments. This 
#' list should be used if one wants to modify many settings for multiple 
#' simulations.
#' @format list
"configData"
