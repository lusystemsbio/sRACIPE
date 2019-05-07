
#' @export
#' @title Simulate a gene regulatory circuit
#' @import SummarizedExperiment
#' @importFrom utils read.table write.table data
#' @importFrom S4Vectors metadata
#' @description Simulate a gene regulatory circuit using its topology as the
#' only input. It will generate an ensemble of random models.
#' @param circuit data.frame or character. The file containing the circuit or
#' @param config  (optional) List. It contains simulation parameters 
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
#' @param anneal (optional) Logical. Default FALSE. Whether to use annealing
#' for stochastic simulations. If TRUE, the gene expressions at higher noise
#' are used as initial conditions for simulations at lower noise.
#' @param numModels (optional) Integer. Default 2000. Number of random models
#' to be simulated.
#' @param thresholdModels (optional)  integer. Default 5000. The number of
#' models to be used for calculating the thresholds for genes.
#' @param paramRange (optional) numeric (0-100). Default 100. The relative
#' range of parameters (production rate, degradation rate, fold change).
#' @param prodRateMin (optional) numeric. Default 1. Minimum production rate.
#' @param prodRateMax (optional) numeric. Default 100. Maximum production rate.
#' @param degRateMin (optional) numeric. Default 0.1. Minimum degradation rate.
#' @param degRateMax (optional) numeric. Default 1. Maximum degradation rate.
#' @param foldChangeMin (optional) numeric. Default 1. Minimum fold change for
#' interactions.
#' @param foldChangeMax (optional) numeric. Default 100. Maximum fold change for
#' interactions.
#' @param hillCoeffMin (optional) integer. Default 1. Minimum hill coefficient.
#' @param hillCoeffMax (optional) integer. Default 6. Maximum hill coefficient.
#' @param integrateStepSize (optional) numeric. Default 0.02. step size for
#' integration using "EM" and "RK4" steppers.
#' @param simulationTime (optional) numeric. Total simulation time.
#' @param nIC (optional) integer. Default 1. Number of initial conditions to be
#' simulated for each model.
#' @param ... Other arguments
#' @param nNoise (optional) integer. Default 0.
#' Number of noise levels at which simulations
#' are to be done. Use nNoise = 1 if simulations are to be carried out at a
#' specific noise. If nNoise > 0, simulations will be carried out at nNoise
#' levels as well as for zero noise. "EM" stepper will be used for simulations
#' and any argument for stepper will be ignoired.
#' @param initialNoise (optional) numeric.
#' Default 50/sqrt(number of genes in the circuit). The initial value of noise
#' for simulations. The noise value will decrease by a factor
#' \code{noiseScalingFactor} at subsequent noise levels.
#' @param noiseScalingFactor (optional) numeric (0-1) Default 0.5.
#' The factor by which noise
#' will be decreased when nNoise > 1.
#' @param shotNoise (optional) numeric. Default 0.
#' The ratio of shot noise to additive
#' noise.
#' @param scaledNoise (optional) logical. Default FALSE. Whether to scale the
#' noise in each gene by its expected median expression across all models. If
#' TRUE the noise in each gene will be proportional to its expression levels.
#' @param outputPrecision (optional) integer. Default 12.
#' The decimal point precison of
#' the output.
#' @param knockOut (optional) List of character or vector of characters.
#' Simulation after knocking out one or more genes. To knock out all the genes
#' in the circuit, use \code{knockOut = "all"}. If it is a vector, then all
#' the genes in the vector will be knocked out simultaneously.
#' @param printStart (optional) numeric (0-\code{simulationTime}). 
#' Default \code{simulationTime}. Its use should be avoided. 
#' The time from which the output should be recorded. Useful for time series
#' analysis and studying the dynamics of a model for a particular initial
#' condition. 
#' @param printInterval (optional) numeric (\code{integrateStepSize}-
#' \code{simulationTime - printStart}). Default 10. The separation between
#' two recorded time points for a given trajectory.
#' @param stepper (optional) Character. Stepper to be used for integrating the
#' differential equations. The options include \code{"EM"} for Euler-Maruyama
#' O(1), \code{"RK4"}
#' for fourth order Runge-Kutta O(4) and \code{"DP"} for adaptive stepper based
#' Dormand-Prince algorithm. The default method is \code{"RK4"}
#'  for deterministic
#' simulations and the method defaults to \code{"EM"}
#' for stochastic simulations.
#' @param plots (optional) logical Default \code{FALSE}.
#' Whether to plot the simuated data.
#' @param plotToFile (optional) Default \code{FALSE}. Whether to save the plots
#' to a file.
#' @param genIC (optional) logical. Default \code{TRUE}. Whether to generate
#' the initial conditions. If \code{FALSE}, the initial conditions must be
#' supplied as a dataframe to \code{circuit$ic}.
#' @param genParams (optional) logical. Default \code{TRUE}. Whether to generate
#' the parameters. If \code{FALSE}, the parameters must be
#' supplied as a dataframe to \code{circuit$params}.
#' @param integrate (optional) logical. Default \code{TRUE}. Whether to
#' integrate the differential equations or not. If \code{FALSE}, the function
#' will only generate the parameters and initial conditions. This can be used
#' iteratively as one can fist generate the parameters and initial conditions
#' and then modify these before using these modified values for integration.
#' For example, this can be used to knockOut genes by changing the production
#' rate and initial condition to zero.
#' @param rkTolerance (optional) numeric. Default \code{0.01}. Error tolerance
#' for adaptive integration method.
#' @return List. A list containing the circuit, parameters, initial conditions,
#' simulated gene expressions, and simulation configuration.
#' @examples
#' data("demoCircuit")
#' rSet <- sRACIPE::sracipeSimulate(circuit = demoCircuit)
#' @section Related Functions:
#'
#' \code{\link{sracipeSimulate}},  \code{\link{sracipeKnockDown}},
#' \code{\link{sracipeOverExp}},  \code{\link{sracipePlotData}}
#'

sracipeSimulate <- function( circuit="inputs/test.tpo", config = config,
                      anneal=FALSE, knockOut = NA_character_,numModels=2000,
                      paramRange=100,
                      prodRateMin=1.,prodRateMax=100, degRateMin=0.1,
                      degRateMax=1.,foldChangeMin=1,foldChangeMax=100,
                      hillCoeffMin=1L,hillCoeffMax=6L,integrateStepSize=0.02,
                      simulationTime=50.0,nIC=1L,nNoise=0L,
                      initialNoise=50.0,noiseScalingFactor=0.5,shotNoise=0,
                      scaledNoise=FALSE,outputPrecision=12L,
                      printStart = 50.0,
                      printInterval=10, stepper = "RK4",
                      thresholdModels = 5000, plots = FALSE, 
                      plotToFile = FALSE,
                      genIC = TRUE, genParams = TRUE,
                      integrate = TRUE, rkTolerance = 0.01, ...){
 rSet <- RacipeSE()
 metadataTmp <- metadata(rSet)
 configData <- NULL
 data("configData",envir = environment())
 configuration <- configData
  if(methods::is(circuit,"RacipeSE"))
  {

      rSet <- circuit
      metadataTmp <- metadata(circuit)
      configuration <- metadata(circuit)$config
  }
   if((methods::is(circuit,"character")) | (methods::is(circuit, "data.frame")))
  {
    sracipeCircuit(rSet) <- circuit
    if(methods::is(circuit, "data.frame")){
      metadataTmp$annotation <- deparse(substitute(circuit))
      metadataTmp$nInteractions <- length(circuit[,1])
      }

   }

 if(missing(circuit)){
    message("Please specify a circuit either as a file or as a data.frame")
    return()
  }
if(!missing(config)){
 if(methods::is(config,"list")){
   configuration <- config
 }
  else{
    message("Incorrect config provided!")
    return()
  }
}

 if(!missing(knockOut)){
   if(!(methods::is(knockOut, "list"))){
     knockOut <- list(knockOut)
   }
 }

  if(!missing(anneal)){
    #print("test1")
 if(anneal)
      configuration$options["anneal"] <- anneal
  }
  if(!missing(numModels)){
    configuration$simParams["numModels"] <- numModels
  }

  if(!missing(paramRange)){
    configuration$simParams["paramRange"] <- paramRange
  }
  if(!missing(prodRateMin)){
    configuration$hyperParams["prodRateMin"] <- prodRateMin
  }
  if(!missing(prodRateMax)){
    configuration$hyperParams["prodRateMax"] <- prodRateMax
  }

  if(!missing(degRateMin)){
    configuration$hyperParams["degRateMin"] <- degRateMin
  }
  if(!missing(degRateMax)){
    configuration$hyperParams["degRateMax"] <- degRateMax
  }
  if(!missing(foldChangeMin)){
    configuration$hyperParams["foldChangeMin"] <- foldChangeMin
  }
  if(!missing(foldChangeMax)){
    configuration$hyperParams["foldChangeMax"] <- foldChangeMax
  }
  if(!missing(hillCoeffMin)){
    configuration$hyperParams["hillCoeffMin"] <- hillCoeffMin
  }
  if(!missing(hillCoeffMax)){
    configuration$hyperParams["hillCoeffMax"] <- hillCoeffMax
  }
  if(!missing(integrateStepSize)){
    configuration$simParams["integrateStepSize"] <- integrateStepSize
  }
  if(!missing(simulationTime)){
    configuration$simParams["simulationTime"] <- simulationTime
  }
  if(!missing(nIC)){
    configuration$simParams["nIC"] <- nIC
  }
  if(!missing(outputPrecision)){
    configuration$simParams["outputPrecision"] <- outputPrecision
  }

  if(!missing(nNoise)){
    configuration$stochParams["nNoise"] <- nNoise
  }
  if(!missing(initialNoise)){
    configuration$stochParams["initialNoise"] <- initialNoise
  }
  if(!missing(noiseScalingFactor)){
    configuration$stochParams["noiseScalingFactor"] <- noiseScalingFactor
  }
  if(!missing(shotNoise)){
    configuration$stochParams["shotNoise"] <- shotNoise
  }
  if(!missing(scaledNoise)){
    configuration$options["scaledNoise"] <-scaledNoise
  }
  if(!missing(printStart)){
    configuration$simParams["printStart"] <- printStart
  }
  if(!missing(printInterval)){
    configuration$simParams["printInterval"] <- printInterval
  }

 if(!missing(genIC)){
   configuration$options["genIC"] <- genIC
 }

 if(!missing(genParams)){
   configuration$options["genParams"] <- genParams
 }

 if(!missing(integrate)){
   configuration$options["integrate"] <- integrate
 }
 if(missing(printStart)){
  configuration$simParams["printStart"] <-
    configuration$simParams["simulationTime"]
 }
 if(!missing(rkTolerance)){
   configuration$simParams["rkTolerance"] <- rkTolerance
 }
 # Apply parameter range
  configuration$hyperParams["prodRateMin"] <- 0.5*(
    configuration$hyperParams["prodRateMin"] +
      configuration$hyperParams["prodRateMax"]) - 0.5*(
      configuration$hyperParams["prodRateMax"] -
        configuration$hyperParams["prodRateMin"])*
    configuration$simParams["paramRange"]/100
  configuration$hyperParams["prodRateMax"] <- 0.5*(
    configuration$hyperParams["prodRateMin"] +
      configuration$hyperParams["prodRateMax"]) + 0.5*(
      configuration$hyperParams["prodRateMax"] -
        configuration$hyperParams["prodRateMin"])*
    configuration$simParams["paramRange"]/100

  configuration$hyperParams["degRateMin"] <- 0.5*(
    configuration$hyperParams["degRateMin"] +
      configuration$hyperParams["degRateMax"]) - 0.5*(
      configuration$hyperParams["degRateMax"] -
        configuration$hyperParams["degRateMin"])*
    configuration$simParams["paramRange"]/100
  configuration$hyperParams["degRateMax"] <- 0.5*(
    configuration$hyperParams["degRateMin"] +
      configuration$hyperParams["degRateMax"]) + 0.5*(
      configuration$hyperParams["degRateMax"] -
        configuration$hyperParams["degRateMin"])*
    configuration$simParams["paramRange"]/100

  configuration$hyperParams["foldChangeMin"] <- 0.5*(
    configuration$hyperParams["foldChangeMin"] +
      configuration$hyperParams["foldChangeMax"]) - 0.5*(
      configuration$hyperParams["foldChangeMax"] -
        configuration$hyperParams["foldChangeMin"])*
    configuration$simParams["paramRange"]/100
  configuration$hyperParams["foldChangeMax"] <- 0.5*(
    configuration$hyperParams["foldChangeMin"] +
      configuration$hyperParams["foldChangeMax"]) + 0.5*(
      configuration$hyperParams["foldChangeMax"] -
        configuration$hyperParams["foldChangeMin"])*
    configuration$simParams["paramRange"]/100

  nGenes <- length(rSet)
  geneInteraction <- as.matrix(rowData(rSet))

  thresholdGene <- rep(0, nGenes)

  geneNames <- names(rSet)
  # outFile <- paste0(Sys.Date(),"_",metadata(rSet)$annotation,"_",
  #                  basename(tempfile()))
  outFileGE <- tempfile(fileext = ".txt")
  outFileParams <- tempfile(fileext = ".txt")
  outFileIC <- tempfile(fileext = ".txt")
  if(genParams){
  message("Generating gene thresholds")
  #Rcpp::sourceCpp("src/thresholdGenerator.cpp")
  threshold_test <- generateThresholds(
    geneInteraction,  thresholdGene,  configuration)
  configuration$thresholds <- thresholdGene
  if(threshold_test!=0){
    stop("Error in threshold generation")
    }
  } else {

    if(is.null(colData(rSet))){
    message("Please specify the parameters as genParams is FALSE")
    return(rSet)
  } else {
    params <- as.data.frame(sracipeParams(rSet))

    utils::write.table(params, file = outFileParams,
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  }

    if(!genIC){
      if(is.null(colData(rSet))){
      message("Please specify the initial conditions
              as genIC is FALSE")
      return(rSet)
    } else {
      ic <- as.data.frame(sracipeIC(rSet))
      utils::write.table(ic, file = outFileIC,
                  sep = "\t", quote = FALSE, row.names = FALSE,
                  col.names = FALSE)
    }
  }
if(missing(nNoise)){
  if(!missing(initialNoise) & !anneal)
    configuration$stochParams["nNoise"] <- 2
}
  if(anneal){
    if(missing(nNoise)) {configuration$stochParams["nNoise"] <- 30L}

  }

  if(configuration$stochParams["nNoise"] > 0) {
    if(stepper != "EM"){
      warnings("Defaulting to EM stepper for stochastic simulations")
      stepper <- "EM"
    }
    if(missing(initialNoise)) { configuration$stochParams["initialNoise"] <-
      as.numeric(50/sqrt(nGenes))}
  }

  stepperInt <- 1L
  if(stepper == "RK4"){ stepperInt <- 4L}
  if(stepper == "DP") {stepperInt <- 5L}


  configuration$stepper <- stepper
  annotationTmp <- outFileGE
  message("Running the simulations")
  Time_evolution_test<- simulateGRCCpp(geneInteraction, configuration,outFileGE,
                                       outFileParams,outFileIC, stepperInt)

    if(configuration$options["integrate"]){
  #  geFile <- paste0("tmp/",outFile,"_geneExpression.txt")
    geneExpression <- utils::read.table(outFileGE, header = FALSE)
    if(configuration$stochParams["nNoise"] == 0){
      geneExpression <- t(geneExpression)
      rownames(geneExpression) <- geneNames
    assayDataTmp <- list(deterministic = geneExpression)}
    # colnames(geneExpression) <- geneNames
    if(configuration$stochParams["nNoise"] > 0){
      noiseLevels <- (configuration$stochParams["initialNoise"]*
                        configuration$stochParams["noiseScalingFactor"]^
                        seq(0,configuration$stochParams["nNoise"]-1))
      stochasticSimulations <- as.list(noiseLevels)
      names(stochasticSimulations) <- noiseLevels
      for(i in seq_len(configuration$stochParams["nNoise"])){
        stochasticSimulations[[i]] <- geneExpression[
          , (1+(i-1)*(nGenes)):(i*nGenes)]
        colnames(stochasticSimulations[[i]]) <- geneNames
      }

      geneExpression <- geneExpression[
        ,(1+(configuration$stochParams["nNoise"]*nGenes)):
          ((configuration$stochParams["nNoise"]+1)*nGenes)]

      stochasticSimulations <- lapply(stochasticSimulations,t)

      geneExpression <- t(geneExpression)
      rownames(geneExpression) <- geneNames
      assayDataTmp <- c(list(deterministic = geneExpression),
                        stochasticSimulations)
      metadataTmp$stochasticSimulations <- noiseLevels


    }

    paramName <- sracipeGenParamNames(rSet)
    # paramFile <- paste0("tmp/",outFile,"_parameters.txt")
    parameters <- utils::read.table(outFileParams, header = FALSE)
    colnames(parameters) <- paramName

    # icFile <- paste0("tmp/",outFile,"_IC.txt")
    ic <- utils::read.table(outFileIC, header = FALSE)
    colnames(ic) <- geneNames

    metadataTmp$normalized <- FALSE

    ## Knockouts

      if(!missing(knockOut)){
        if(knockOut[[1]] == "all"){
          knockOut <- as.list(geneNames)
        }
        knockOutData <- knockOut
        names(knockOutData) <- knockOut

        for(ko in seq_along(knockOut)){
          koGene <- knockOut[[ko]]
          knockOut_number <- which(koGene==geneNames)
          if(length(knockOut_number)==0){
            message("knockOut gene not found in the circuit")
            return(rSet)
          }
          params <- sracipeParams(rSet)
          ic <- sracipeIC(rSet)
          ic[,knockOut_number] <- 0
          params[,knockOut_number] <- 0

          utils::write.table(params, file = outFileParams,
                      sep = "\t", quote = FALSE, row.names = FALSE,
                      col.names = FALSE)
          utils::write.table(ic, file = outFileIC,
                      sep = "\t", quote = FALSE, row.names = FALSE,
                      col.names = FALSE)

          Time_evolution_test<- simulateGRCCpp(geneInteraction, configuration,
            outFileGE,outFileParams,outFileIC, stepperInt)


          # geFile <- paste0("tmp/",outFileKO,"_geneExpression.txt")
          geneExpression <- utils::read.table(outFileGE, header = FALSE)
          colnames(geneExpression) <- geneNames
          knockOutData[[ko]] <- geneExpression

        }
        knockOutData <- lapply(knockOutData, t)
        assayDataTmp <- c(assayDataTmp,knockOutData)
        metadataTmp$knockOutSimulations <- names(knockOutData)

  }
    }
  else {
      paramName <- sracipeGenParamNames(rSet)
      parameters <- utils::read.table(outFileParams, header = FALSE)
      colnames(parameters) <- paramName
      
      ic <- utils::read.table(outFileIC, header = FALSE)
      colnames(ic) <- geneNames
      colData <- (cbind(parameters,ic))
      metadataTmp$config <- configuration
      rSet <- RacipeSE(rowData = geneInteraction, colData = colData,
                       
                       metadata = metadataTmp)
      return(rSet)
  }
  
    colData <- (cbind(parameters,ic))
    metadataTmp$config <- configuration
 # return(list(rowData = geneInteraction, colData = colData,
   #         assays =  assayDataTmp,
  #          metadata = metadataTmp))
    rSet <- RacipeSE(rowData = geneInteraction, colData = colData,
                     assays =  assayDataTmp,
                     metadata = metadataTmp)
##############################
 if(plots){
    rSet <- sracipePlotData(rSet,plotToFile = plotToFile,...)
  }

return(rSet)
}



