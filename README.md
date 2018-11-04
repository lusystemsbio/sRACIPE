![](/www/racipe.png)
# sRACIPE 

*Stochastic analysis for Random Circuit Perturbation*



sRACIPE implements a randomization-based method for gene circuit modeling. It allows us to study the effect of both the gene expression noise and the parametric variation on any gene regulatory circuit (GRC) using only its topology, and simulates an ensemble of models with random kinetic parameters at multiple noise levels. Statistical analysis of the generated gene expressions reveals the basin of attraction and stability of various phenotypic states and their changes associated with intrinsic and extrinsic noises. sRACIPE provides a holistic picture to evaluate the effects of both the stochastic nature of cellular processes and the parametric variation.   

If you use sRACIPE, please consider citing our paper [Role of noise and parametric variation in the dynamics of gene regulatory circuits](https://www.nature.com/articles/s41540-018-0076-x) published in [npj Systems Biology and Applications](https://www.nature.com/npjsba/articles).

*Installing R package from GitHub*

1) Install devtools package. In R, type: 
	install.packages("devtools")
2) Load the devtools package:
	library(devtools)
3) Install sRACIPE from GitHub:
	install_github("lusystemsbio/sRACIPE")

sRACIPE is a self contained package and includes all the dependencies in the packrat folder.
Packrat can be installed using install.packages("packrat").

*Using the package:* 
Please refer to the  [Using sRACIPE](http://htmlpreview.github.io/?https://github.com/lusystemsbio/sRACIPE/blob/master/man/Using_sRACIPE.html) and vignettes. 

 
