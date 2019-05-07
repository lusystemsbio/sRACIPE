![](/www/racipe.png)
# sRACIPE 

*Stochastic analysis for Random Circuit Perturbation*

sRACIPE implements a randomization-based method for gene circuit modeling. It allows us to study the effect of both the gene expression noise and the parametric variation on any gene regulatory circuit (GRC) using only its topology, and simulates an ensemble of models with random kinetic parameters at multiple noise levels. Statistical analysis of the generated gene expressions reveals the basin of attraction and stability of various phenotypic states and their changes associated with intrinsic and extrinsic noises. sRACIPE provides a holistic picture to evaluate the effects of both the stochastic nature of cellular processes and the parametric variation.   

If you use sRACIPE, please consider citing our paper [Role of noise and parametric variation in the dynamics of gene regulatory circuits](https://www.nature.com/articles/s41540-018-0076-x) published in [npj Systems Biology and Applications](https://www.nature.com/npjsba/articles).

*Installation*

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("sRACIPE")

Or install the development version of the package from Github.

library(devtools)
install_github(“lusystembio/sRACIPE”)

# Webserver

For quick exploration of circuits, one can simulate and visualize them in the browser using our webserver [GeneEx](https://shinyapps.jax.org/5c965c4b284ca029b4aa98483f3da3c5/).
