# File for generating NAMESPACE and package info using roxygen

#' sRACIPE: A package for stochastic random circuit perturbation.
#'
#' sRACIPE is a systems biology tool to study the role of noise and parameter
#' variation in gene regulatory circuits. It implements a randomization-based
#' method for gene circuit modeling. It allows us to study the effect of both
#' the gene expression noise and the parametric variation on any
#' gene regulatory circuit (GRC) using only its topology, and simulates an
#' ensemble of models with random kinetic parameters at multiple noise levels.
#' Statistical analysis of the generated gene expressions reveals the basin of
#' attraction and stability of various phenotypic states and their changes
#' associated with intrinsic and extrinsic noises. sRACIPE provides a holistic
#' picture to evaluate the effects of both the stochastic nature of cellular
#' processes and the parametric variation.
#'
#' @section sRACIPE functions:
#'
#' \code{\link{sracipeSimulate}}
#' Primary function to simulate a circuit.
#'  Contains options for plotting as well.
#' \code{\link{sracipeKnockDown}}
#' In-silico knockdown analysis of the circuit. Plots the relative changes in
#' different cluster proportions.
#'
#' \code{\link{sracipeOverExp}}
#' In-silico over expression analysis of the circuit.
#' Plots the relative changes in
#' different cluster proportions.
#'
#' \code{\link{sracipePlotData}}
#' Plot the simulated data. Includes options to plot the hierarchichal
#' clustering analysis, principal components, and uniform manifold
#' approximation and projection. Can plot the stochastic as well as the
#' knockout simulations.
#'
#' \code{\link{sracipeSimulate}},  \code{\link{sracipeKnockDown}},
#' \code{\link{sracipeOverExp}},  \code{\link{sracipePlotData}}
#'
#' @import visNetwork reshape2 ggplot2 gplots MASS RColorBrewer
#' @import htmlwidgets Rcpp
#' @useDynLib sRACIPE
#' @docType package
#' @name sRACIPE
NULL
