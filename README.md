# sRACIPE 

*Suite for Random Circuit Perturbation Analysis*


sRACIPE implements a randomization-based method for gene circuit modeling. It allows us to study the effect of both the gene expression noise and the parametric variation on any gene regulatory circuit (GRC) using only its topology, and simulates an ensemble of models with random kinetic parameters at multiple noise levels. Statistical analysis of the generated gene expressions reveals the basin of attraction and stability of various phenotypic states and their changes associated with intrinsic and extrinsic noises. sRACIPE provides a holistic picture to evaluate the effects of both the stochastic nature of cellular processes and the parametric variation.   

If you use sRACIPE, please consider citing our paper [Role of noise and parametric variation in the dynamics of gene regulatory circuits](https://www.nature.com/articles/s41540-018-0076-x) published in [npj Systems Biology and Applications](https://www.nature.com/npjsba/articles).

## Installation ##

```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("sRACIPE")
```
Or install the development version of the package from Github.
```
library(devtools)
install_github(“lusystemsbio/sRACIPE”)
```
## Tutorials ##
More tutorials and applications are available [here](https://vivekkohar.github.io/sRACIPE/) 


## Webserver ##

For quick exploration of circuits, one can simulate and visualize them in the browser using our webserver [GeneEx](https://geneex.jax.org/).

## Build Reports ##
- [![Bioconductor Release](https://bioconductor.org/shields/build/release/bioc/sRACIPE.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/sRACIPE/) – Bioconductor Release
- [![Bioconductor devel](https://bioconductor.org/shields/build/devel/bioc/sRACIPE.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/sRACIPE/) – Bioconductor devel

