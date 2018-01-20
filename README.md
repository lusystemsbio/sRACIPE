# sRACIPE
Stochastic analysis for Random Circuit Perturbation

sRACIPE implements a randomization-based method for gene circuit modeling. It allows us to study the effect of both the gene expression noise and the parametric variations on any gene regulatory circuit (GRC) using only its topology. 

sRACIPE provides a holistic picture to evaluate the effects of both the stochastic nature of cellular processes and the parametric variations. We define a distinct cellular state as one of the clusters of steady state gene expression profiles from random models and evaluate how gene expression noise affects the formation and expression patterns of clusters. The stochastic analysis can quantify the relative stability of the steady states for systems allowing multiple states. 

Installing R package from GitHub

1) Install devtools package. In R, type: 
	install.packages("devtools")
2) Load the devtools package:
	library(devtools)
3) Install sRACIPE from GitHub:
	install_github("lusystemsbio/sRACIPE")
 
