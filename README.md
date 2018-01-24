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

sRACIPE is a self contained package and includes all the dependencies in the packrat folder.
Packrat can be installed using install.packages("packrat").

Using the package: 
In the working directory, create an "input" folder which should contain the topology file "*.tpo".
Call the sRACIPE_stochastic() function. 
Results for various noise levels will be in the results folder in the same directory.
Vignettes and information on other functions will be updated later.
 
