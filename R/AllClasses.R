#' @export
#' @import methods
#' @import SummarizedExperiment
#' @title RacipeSE
#' @description An S4 class for Random Circuit Perturbation (RACIPE)
#' simulations of networks. Extends the \link{SummarizedExperiment} class.
#' RACIPE can simulate a gene regulatory circuit using the circuit and 
#' a large ensemble of parameters. 
#'
#'
.RacipeSE <- setClass("RacipeSE",
                      contains = "SummarizedExperiment")

#' @export
#' @title RacipeSE constructor
#' @description Create an RacipeSE object. RacipeSE is an S4 class for 
#' Random Circuit Perturbation (RACIPE) simulations of networks in which a large
#' number of models with randomized parameters are used for simulation of the
#' circuit. Each model can be considered as a sample. 
#' It extends the \link{SummarizedExperiment} class to store and access 
#' the circuit, simulated gene expressions, parameters, intial conditions and
#' other meta information. 
#' SummarizedExperiment slot assays is used for storing simulated 
#' gene expressions. The rows of these
#' matrix-like elements correspond to various genes in the circuit and columns
#' correspond to models. 
#' The first element is used for unperturbed deterministic 
#' simulations. The subsequent elements are used for stochastic simulations
#' at different noise levels and/or knockout simulations. 
#' SummarizedExperiment slot rowData stores the circuit topology. It is a square 
#' matrix with dimension equal to the number of genes in the circuit. The values
#' of the matrix represent the type of interaction in the gene pair given by
#' row and column. 1 represents activation, 2 inhibition and 0 no interaction.
#' This should not be set directly and  \code{\link{sracipeCircuit}}
#' accessor should be used instead.
#' SummarizedExperiment slot colData contains the parameters 
#' and initial conditions for each 
#' model. Each gene in the circuit has two parameters, namely, its production
#' rate and its degradation rate. Each interaction in the has three parameters, 
#' namely, threshold of activation, the hill coefficient, and the fold change.
#' Each gene has one or more initial gene expression values as specified
#' by nIC. This should not be modified directly and \code{\link{sracipeParams}}
#' and \code{\link{sracipeIC}} accessors should be used instead.
#' SummarizedExperiment slot metadata Contains metadata 
#' information especially the config list 
#' (containing the simulation settings), annotation, nInteraction (number of
#' interactions in the circuit), normalized (whether the data is normalized or 
#' not), data analysis lists like pca, umap, cluster assignment of the models
#' etc. The config list includes simulation parameters like integration method
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
#' modifiable simulation options can be specified as arguments to 
#' \code{\link{sracipeSimulate}} function.
#'    
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata DataFrame isEmpty
#' @importFrom utils data
#' @param .object (optional) Another RacipeSE object.
#' @param assays (optional) assay object for initialization
#' @param rowData (optional) rowData for initialization
#' @param colData (optional) colData for initialization
#' @param metadata (optional) metadata for initialization 
#' @param ... Arguments passed to SummarizedExperiment
#' @return RacipeSE object
#' @examples
#' rSet <- RacipeSE()
#'

RacipeSE <- function(.object = NULL, assays = SimpleList(),
                     rowData = NULL,
                     colData = DataFrame(),
                     metadata = list(), ...) {
    
    if(is(.object,"RacipeSE")) {
        if(isEmpty(assays)) assays <- assays(.object)
        if(is.null(rowData)) rowData <- rowData(.object)
        if(isEmpty(colData)) colData <- colData(.object)
        if(length(metadata)==0) metadata <- metadata(.object)
    } 
    objectTmp <- SummarizedExperiment(
        assays=assays, rowData=rowData, colData=colData, metadata=metadata)

    if(is.null(metadata(objectTmp)$config))
    {
        configData <- NULL
        data("configData",envir = environment(), package = "sRACIPE")
        metadata(objectTmp)$config <- configData
    }

    ## Validity steps
    metadata(objectTmp)$config$simParams["numModels"] <- dim(objectTmp)[2]
    
    if(!is.null(metadata(objectTmp)$nInteractions)) 
        metadata(objectTmp)$nInteractions <- 
        sum(as.matrix(rowData(objectTmp))>0)
    
    .RacipeSE(objectTmp)
}

