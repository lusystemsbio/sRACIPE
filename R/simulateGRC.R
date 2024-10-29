
#' @export
#' @title Simulate a gene regulatory circuit
#' @import SummarizedExperiment
#' @import doRNG
#' @import doFuture
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
#' @param simulationTime (optional) numeric. Default 50. Total simulation time.
#' Only used for stochastic and time series simulations. For adjusting
#' deterministic simulation run times, see numConvergenceIter.
#' @param nIC (optional) integer. Default 1. Number of initial conditions to be
#' simulated for each model.
#' @param ... Other arguments
#' @param nNoise (optional) integer. Default 0.
#' Number of noise levels at which simulations
#' are to be done. Use nNoise = 1 if simulations are to be carried out at a
#' specific noise. If nNoise > 0, simulations will be carried out at nNoise
#' levels as well as for zero noise. "EM" stepper will be used for simulations
#' and any argument for stepper will be ignoired.
#' @param simDet (optional) logical. Default TRUE.
#' Whether to simulate at zero noise as well also when using nNoise > 0.
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
#' Default \code{simulationTime}. To be used only when \code{timeSeries} is
#' \code{TRUE}.
#' The time from which the output should be recorded. Useful for time series
#' analysis and studying the dynamics of a model for a particular initial
#' condition.
#' @param printInterval (optional) numeric (\code{integrateStepSize}-
#' \code{simulationTime - printStart}). Default 10. The separation between
#' two recorded time points for a given trajectory.
#' To be used only when \code{timeSeries} is
#' \code{TRUE}.
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
#' @param timeSeries (optional) logical. Default \code{FALSE}.
#' Whether to generate time series for a single model instead of performing
#' RACIPE simulations.
#' @param signalRate (optional) numeric. Default \code{10}. The factor by which
#' the differential equations for fast genes are multiplied by. Fast genes only
#' have inward signaling interactions of types 5 or 6
#' @param uniqueDigits (optional) integer. Default \code{4}. Deterministic
#' simulations with nIC > 1 count the number of unique steady states per model
#' by truncating the expression by this many digits and counting the number of
#' unique expressions
#' @param convergThresh (optional) numeric. Default \code{1e-12}. The threshold
#' for convergence to a steady state for deterministic simulations.
#' @param numStepConverge (optional) integer Default \code{500}. The number of
#' integration steps between convergence tests for deterministic simulations.
#' @param numConvergenceIter (optional) integer. Default \code{25}. The total
#' number of convergence test iterations to run per model
#' initial condition in deterministic simulations.
#' @param limitcycles (optional) logical. Default \code{FALSE}. Whether to check
#' for limit cycles in deterministic simulations.
#' @param LCSimTime (optional) numeric. Default \code{10}. The length of each
#' iteration for the secondary limit cycle simulation
#' @param LCSimStepSize (optional) numeric. Default \code{0.01}. The integration
#' step size for the secondary limit cycle simulation.
#' @param maxLCs (optional) integer. Default \code{10}. The maximum allowable
#' number of limit cycles that can be detected for each model.
#' @param LCIter (optional) integer. Default \code{20}. The number of iterations
#' to run in the secondary limit cycle simulation
#' @param MaxPeriods (optional) integer. Default \code{100}. The number of periods
#' to count in the distance function constructed by the limit cycle detection
#' algorithm
#' @param NumSampledPeriods (optional) integer. Default \code{3}. The number of
#' times the limit cycle detection algorithm tries to calculate the period of the
#' simulated trajectory using the local minima in a constructed distance function.
#' @param AllowedPeriodError (optional) integer. Default \code{3}. The allowed
#' difference in the sampled periods for a trajectory from the limit cycle
#' detection algorithm. Decrease this value to make for more stringent limit
#' cycle detection
#' @param SamePointProximity (optional) numeric. Default \code{0.1}. The max
#' allowed square euclidean difference between two minima of the distance function
#' constructed by the limit cycle detection algorithm before they are considered
#' different points in phase space.
#' @param LCStepper (optional) Character. Default \code{"RK4"} The integration
#' method used for the limit cycle simulation. The options include \code{"E"}
#' for first order Euler and \code{"RK4"} for fourth order Runge Kutta. Any other
#' input will cause the simulation to use the Euler method
#' @param paramSignalVals (optional) Data Frame. Default data.frame(). The
#' first column must be a vector of time values with the first element as 0 and
#' the last element as simulationTime. The other column names must be valid
#' production or degradation parameter names for the circuit. Use
#' sracipeGenParamNames() to generate valid parameter names. The
#' first column must be a vector of time values with the first element as 0 and
#' the last element as simulationTime.
#' @param geneClamping (optional) Data Frame. Default data.frame(). The column
#' names must be genes in the circuit. The number of columns must either be one
#' or equal to numModels. If the number of columns is one, the selected genes
#' are clamped to those values for every model. Otherwise, the gene is clamped
#' to the value of the corresponding row for a particular model.
#'  @param nCores (optional) integer. Default \code{1}
#' Number of cores to be used for computation. Utilizes \code{multisession} from
#' \code{doFuture} pacakge. Will not work in Rstudio.
#' @return \code{RacipeSE} object. RacipeSE class inherits
#' \code{SummarizedExperiment} and contains the circuit, parameters,
#' initial conditions,
#' simulated gene expressions, and simulation configuration. These can be
#' accessed using correponding getters.
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
                      simulationTime=50.0,nIC=1L,nNoise=0L,simDet = TRUE,
                      initialNoise=50.0,noiseScalingFactor=0.5,shotNoise=0,
                      ouNoise_t=1, scaledNoise=FALSE,outputPrecision=12L,
                      printStart = 50.0,
                      printInterval=10, stepper = "RK4",
                      thresholdModels = 5000, plots = FALSE,
                      plotToFile = FALSE,
                      genIC = TRUE, genParams = TRUE,
                      integrate = TRUE, rkTolerance = 0.01, timeSeries = FALSE,
                      signalRate = 10.0, uniqueDigits = 4, convergThresh = 1e-12,
                      numStepsConverge = 500, numConvergenceIter = 25,
                      limitcycles = FALSE, LCSimTime = 10, LCSimStepSize = 0.01,
                      maxLCs = 10, LCIter = 20, MaxPeriods = 100,
                      NumSampledPeriods = 3, AllowedPeriodError = 3,
                      SamePointProximity = 0.1, LCStepper = "RK4",
                      paramSignalVals = data.frame(), geneClamping = data.frame(),
                      nCores = 1L,
                      ...){
 rSet <- RacipeSE()
 metadataTmp <- metadata(rSet)
 configData <- NULL
 data("configData",envir = environment(), package = "sRACIPE")
 configuration <- configData

  if(methods::is(circuit,"RacipeSE"))
  {

      rSet <- circuit
      metadataTmp <- metadata(rSet)
      configuration <- metadata(rSet)$config
  }
   if((methods::is(circuit,"character")) | (methods::is(circuit, "data.frame")))
  {
    sracipeCircuit(rSet) <- circuit
    circuitTmp <- sracipeCircuit(rSet)
    metadataTmp$geneTypes <- metadata(rSet)$geneTypes
    metadataTmp$nInteractions <- length(circuitTmp[,1])

    if(methods::is(circuit, "data.frame")){
      metadataTmp$annotation <- deparse(substitute(circuit))
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

 if(nCores<1) {
   warning("Number of cores, nCores, is less than 1 or not an interger.
          Using nCores=1")
   nCores=1
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
  if(!missing(signalRate)){
    configuration$hyperParams["signalRate"] <- signalRate
  }
  if(!missing(integrateStepSize)){
    configuration$simParams["integrateStepSize"] <- integrateStepSize
  }
  if(!missing(simulationTime)){
    configuration$simParams["simulationTime"] <- simulationTime
  }
 if(!missing(simDet)){
   configuration$options["simDet"] <- simDet
 } else {
   configuration$options["simDet"] <- TRUE
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
  if(!missing(ouNoise_t)){
    configuration$stochParams["ouNoise_tcorr"] <- ouNoise_t
  }
  if(!missing(scaledNoise)){
    configuration$options["scaledNoise"] <-scaledNoise
  }
  if(!missing(printStart)){
    configuration$simParams["printStart"] <- printStart
  }
  if(!missing(printInterval)){
    configuration$simParams["printInterval"] <- printInterval
    if(configuration$simParams["printInterval"] <
       configuration$simParams["integrateStepSize"]){
      configuration$simParams["printInterval"] <-
      configuration$simParams["integrateStepSize"]
      warnings("Print Interval cannot be smaller than integration step size.
               Setting it to integrate step size.")}
  }
 if(!missing(convergThresh)){
   configuration$simParams["convergThresh"] <- convergThresh
 }
 if(!missing(numStepsConverge)){
   configuration$simParams["numStepsConverge"] <- numStepsConverge
 }
 if(!missing(numConvergenceIter)){
   configuration$simParams["numConvergenceIter"] <- numConvergenceIter
 }
 if(!missing(uniqueDigits)){
   configuration$simParams["uniqueDigits"] <- uniqueDigits
 }

 if(!missing(genIC)){
   configuration$options["genIC"] <- genIC
 }

 if(!missing(genParams)){
   configuration$options["genParams"] <- genParams
 }
# stepper is not included in configdata. This can be changed
 configuration$stepper <- stepper

 if(!missing(integrate)){
   configuration$options["integrate"] <- integrate
 }
 if(!missing(limitcycles)){
   configuration$options["limitcycles"] <- limitcycles
 }
 if(missing(printStart)){
  configuration$simParams["printStart"] <-
    configuration$simParams["simulationTime"]
 }
 if(!missing(rkTolerance)){
   configuration$simParams["rkTolerance"] <- rkTolerance
 }
 #Putting LC params in configuration list
 if(!missing(LCSimTime)){
   configuration$LCParams["LCSimTime"] <- LCSimTime
 }
 if(!missing(LCSimStepSize)){
   configuration$LCParams["LCSimStepSize"] <- LCSimStepSize
 }
 if(!missing(maxLCs)){
   configuration$LCParams["maxLCs"] <- maxLCs
 }
 if(!missing(LCIter)){
   configuration$LCParams["LCIter"] <- LCIter
 }
 if(!missing(MaxPeriods)){
   configuration$LCParams["MaxPeriods"] <- MaxPeriods
 }
 if(!missing(NumSampledPeriods)){
   configuration$LCParams["NumSampledPeriods"] <- NumSampledPeriods
 }
 if(!missing(AllowedPeriodError)){
   configuration$LCParams["AllowedPeriodError"] <- AllowedPeriodError
 }
 if(!missing(SamePointProximity)){
   configuration$LCParams["SamePointProximity"] <- SamePointProximity
 }

 # Storing stepper for limit cycle calculations
  configuration$LCStepper <- LCStepper

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
  outFileConverge <- tempfile(fileext = ".txt")
  outFileLC <- tempfile(fileext = ".txt")
  outFileLCIC <- tempfile(fileext = ".txt") #Records the different ics which produce limitcycles
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
      ic <- as.data.frame(t(sracipeIC(rSet)))
      utils::write.table(ic, file = outFileIC,
                  sep = "\t", quote = FALSE, row.names = FALSE,
                  col.names = FALSE)
    }
  }
if(missing(nNoise)){
  if(!missing(initialNoise) & !anneal )
  {
    if(timeSeries)
    # if(timeSeries)
      {configuration$stochParams["nNoise"] <- 1}
    else {if(simDet) configuration$stochParams["nNoise"] <- 2}
  }

}
  if(anneal){
    if(missing(nNoise)) {configuration$stochParams["nNoise"] <- 30L}

  }
  # For time series, default to 1 model and single noise level.
  if(timeSeries){
    configuration$simParams["numModels"] <- 1
    # configuration$stochParams["nNoise"] <- 0

    if(missing(printInterval))
    configuration$simParams["printInterval"] <- max(
      0.05, configuration$simParams["integrateStepSize"])
    if(missing(printStart)) configuration$simParams["printStart"] <- 0

    if(configuration$stepper == "DP") {configuration$stepper <- "RK4"
    warnings("Using RK4 integration method for time series simulation.")

    }
  }
  #Checking whether or not to do convergence tests
  if(all(simDet, nNoise==0L, !timeSeries, nrow(paramSignalVals)==0)){
    convergTesting <- TRUE
    message("Running with convergence tests")
  }else {convergTesting <- FALSE}

  configuration$options["convergTesting"] <- convergTesting

  if(!convergTesting){
    stepperInt <- 1L
    if(configuration$stepper == "RK4"){ stepperInt <- 4L}
    if(configuration$stepper == "DP") {stepperInt <- 5L}
    if(configuration$stepper == "EM_OU") {stepperInt <- 6L}
  } else { #Two digit stepperInt values correlate to running convergence tests
    stepperInt <- 11L
    if(configuration$stepper == "RK4"){ stepperInt <- 41L}
    if(configuration$stepper == "DP") {stepperInt <- 51L}
  }

  paramSignalTypes <- numeric(nGenes)
  if(nrow(paramSignalVals)>0){ #Determining which genes are being varied, if an input is given
    metadataTmp$paramSignalVals <- paramSignalVals
    paramName <- sracipeGenParamNames(rSet)
    prodParams <- paramName[1:nGenes]
    degParams <- paramName[(nGenes+1):(2*nGenes)]
    variedParams <- colnames(paramSignalVals[,2:ncol(paramSignalVals)]) #First column is time  points

    if(paramSignalVals[1,1] != 0 | paramSignalVals[nrow(paramSignalVals), 1] != simulationTime){
      message("The first and last time points for parameter variation must be 0 and simulationTime respectively")
      return()
    }
    if(!(identical(paramSignalVals[,1], sort(paramSignalVals[,1])))){
      message("Time points for parameter variation must be given in ascending order")
      return()
    }

    for(param in variedParams){
      if(param %in% prodParams){
        geneIdx <- which(prodParams == param)
        paramSignalTypes[geneIdx] <- 1
      }
      else if(param %in% degParams){
        geneIdx <- which(degParams == param)
        if(paramSignalTypes[geneIdx] == 1){
          paramSignalTypes[geneIdx] <- 3 #Both deg and prod being varied
        }
        else {
          paramSignalTypes[geneIdx] <- 2
        }
      }
      else{
        message("Invalid Param names given in paramSignalVals")
        return()
      }

    }
    #Sorting paramSignalVals columns dataframe to be read in C++
    tempIdxArray <- data.frame(value = variedParams,
                               Idx1 = match(variedParams, prodParams),
                               Idx2 = match(variedParams, degParams))
    #Replace NaN vals with Inf for sorting
    tempIdxArray$Idx1[is.na(tempIdxArray$Idx1)] <- Inf
    tempIdxArray$Idx2[is.na(tempIdxArray$Idx2)] <- Inf

    tempIdxArray <- tempIdxArray[order(pmin(tempIdxArray$Idx1, tempIdxArray$Idx2), tempIdxArray$Idx1 < Inf,
                                       tempIdxArray$Idx1, tempIdxArray$Idx2), ]

    variedParams <- tempIdxArray$value
    paramSignalVals <- paramSignalVals[, c(colnames(paramSignalVals)[1],variedParams)] #Sorted for C++ steppers

  }else{
    #If no paramSignalVals dataframe provided, create a dummy data frame to pass to C++
    paramSignalVals <- data.frame(temp1 = c(-1,-1), temp2 = c(-1, -1))
  }
  paramSignalValsTmp <- as.matrix(paramSignalVals)

  #gene clamping
  clampedGenes <- rep(0, nGenes)
  if(nrow(geneClamping)>0){
    clamping <- TRUE
    clampedIdx <- match(colnames(geneClamping), geneNames)
    clampedGenes[clampedIdx] <- 1
    message(paste0("clamped genes: ",paste0(clampedGenes, collapse = ",")))

    if(nrow(geneClamping) == 1){
      clampVals <- as.matrix(do.call(rbind, replicate(numModels, geneClamping, simplify = FALSE)))
      configuration$clampVals <- clampVals
    }else if(nrow(geneClamping) == numModels){
      configuration$clampVals <- as.matrix(geneClamping)
    }else{
      message("geneClamping must have a number of rows equal to 1 or numModels")
      return(rSet)
    }
  }else{
    configuration$clampVals <- as.matrix(c(-1))
    clamping <- FALSE
  }
  configuration$clampedGenes <- clampedGenes

  if(configuration$stochParams["nNoise"] > 0) {
    if(stepper != "EM" && stepper != "EM_OU"){
      warnings("Defaulting to EM stepper for stochastic simulations")
      configuration$stepper <- "EM"
      stepperInt <- 1L
    }
    if(missing(initialNoise)) { configuration$stochParams["initialNoise"] <-
      as.numeric(50/sqrt(nGenes))}
  }





  annotationTmp <- outFileGE
  message("Running the simulations")
  # print(configuration$stochParams["nNoise"])

  if(!configuration$options["integrate"] | (nCores==1) | timeSeries){
    Time_evolution_test<- simulateGRCCpp(geneInteraction, configuration,outFileGE,
                                         outFileParams,outFileIC, outFileConverge,
                                         metadataTmp$geneTypes, paramSignalValsTmp,
                                         paramSignalTypes, stepperInt)

  } else{
      if(nCores>1 & !timeSeries){
        configuration$options["integrate"] <- FALSE
        Time_evolution_test<- simulateGRCCpp(geneInteraction, configuration,outFileGE,
                                             outFileParams,outFileIC,outFileConverge,
                                             metadataTmp$geneTypes, paramSignalValsTmp,
                                             paramSignalTypes, stepperInt)
        configuration$options["integrate"] <- TRUE
        requireNamespace("doFuture")
        registerDoFuture()
        plan(multisession)

        configList <- list()
        parModel <- floor(configuration$simParams["numModels"]/nCores)
        gEFileList <- character(length = nCores)
        paramFileList <- character(length = nCores)
        iCFileList <- character(length = nCores)
        convFileList <- character(length = nCores)

        parameters <- utils::read.table(outFileParams, header = FALSE)
        ic <- utils::read.table(outFileIC, header = FALSE)

        for(i in seq(1,nCores)){
          parConfig <- configuration

          gEFileList[i] <- tempfile(fileext = ".txt")
          paramFileList[i] <- tempfile(fileext = ".txt")
          iCFileList[i] <- tempfile(fileext = ".txt")
          convFileList[i] <- tempfile(fileext = ".txt")
          if(i==nCores){
            parConfig$simParams["numModels"] <-
              configuration$simParams["numModels"] - (nCores-1)*parModel
            paramPar <- parameters[(((i-1)*parModel+1):(nrow(parameters))),]
            icPar <- ic[(((i-1)*parModel+1):(nrow(parameters))),]

            if(clamping){
              parConfig$clampVals <- as.matrix(geneClamping[(((i-1)*parModel+1):(nrow(parameters))),])
            }

          } else {
            parConfig$simParams["numModels"] <- parModel
            paramPar <- parameters[(((i-1)*parModel+1):(i*parModel)),]
            icPar <- ic[(((i-1)*parModel+1):(i*parModel)),]

            if(clamping){
              parConfig$clampVals <- as.matrix(geneClamping[(((i-1)*parModel+1):(i*parModel)),])
            }

          }

          configList[[i]] <- parConfig

          utils::write.table(paramPar, file = paramFileList[i],
                             sep = "\t", quote = FALSE, row.names = FALSE,
                             col.names = FALSE)
          utils::write.table(icPar, file = iCFileList[i],
                             sep = "\t", quote = FALSE, row.names = FALSE,
                             col.names = FALSE)

        }
        print("here1")
        x <- foreach(configurationTmp = configList,outFileGETmp = gEFileList,
                     outFileParamsTmp=paramFileList, outFileICTmp=iCFileList,
                     outFileConvergeTmp=convFileList,
                     .export = c("geneInteraction","stepperInt", "metadataTmp",
                                 "paramSignalValsTmp", "paramSignalTypes")) %dorng% {

                       simulateGRCCpp(
                         geneInteraction, configurationTmp,outFileGETmp, outFileParamsTmp,
                         outFileICTmp, metadataTmp$geneTypes, paramSignalValsTmp,
                         paramSignalTypes, stepperInt)


                                 }
        print("here2")
        geList <- list()
        convList <- list()
        for(i in seq(1,nCores)){

          geList[[i]] <- utils::read.table(iCFileList[i],
                                           header = FALSE)
          convList[[i]] <- utils::read.table(convFileList[i],
                                           header = FALSE)
        }
        geneExpression <- do.call(rbind, geList)
        converge <- do.call(rbind, convList)

      }
    }


    if(configuration$options["integrate"]){


      if(timeSeries){

        timeStamps <- c(seq(
          configuration$simParams["printStart"] +
            configuration$simParams["printInterval"],
          configuration$simParams["simulationTime"],
          configuration$simParams["printInterval"]),
          configuration$simParams["simulationTime"])
        geneExpression <- utils::read.table(outFileGE, header = FALSE)
        # print(dim(geneExpression))
        geneExpression <- matrix(geneExpression, ncol = nGenes, byrow = TRUE)
        # print(dim(geneExpression))
        parameters <- utils::read.table(outFileParams, header = FALSE)
        paramName <- sracipeGenParamNames(rSet)
        colnames(parameters) <- paramName
        ic <- utils::read.table(outFileIC, header = FALSE)
        if(configuration$stochParams["nNoise"] ==1){
          # print(dim(geneExpression))
          detGeneExp <- geneExpression[(1+dim(geneExpression)[1]/2):dim(geneExpression)[1],]
          geneExpression <- geneExpression[(seq_len(dim(geneExpression)[1]/2)),]
          #print(dim(geneExpression))
          detGeneExp <- t(detGeneExp)
          rownames(detGeneExp) <- geneNames
          colnames(detGeneExp) <- timeStamps[seq_len(dim(detGeneExp)[2])]
          metadataTmp$timeSeriesDet <- NULL
          metadataTmp$timeSeriesDet <- detGeneExp
          #print(dim(detGeneExp))

        }

        geneExpression <- t(geneExpression)
        #print(dim(geneExpression))
        colnames(ic) <- geneNames
        metadataTmp$normalized <- FALSE

        colData <- cbind(parameters,ic)
        metadataTmp$config <- configuration
      # colnames(geneExpression) <- seq(
      #  configuration$simParams["integrateStepSize"],
      #  configuration$simParams["simulationTime"],
      # configuration$simParams["integrateStepSize"])
        rownames(geneExpression) <- geneNames

        # print(dim(geneExpression))
        # print(length(timeStamps))
        colnames(geneExpression) <- timeStamps[seq_len(dim(geneExpression)[2])]
        metadataTmp$timeSeries <- NULL
      metadataTmp$timeSeries <- geneExpression

        rSet <- RacipeSE(rowData = geneInteraction,
                         colData = colData,
                         metadata = metadataTmp)
        return(rSet)
      }
    if(!(nCores>1 & !timeSeries)){
      geneExpression <- utils::read.table(outFileGE, header = FALSE)
    }
    if(!all(!(is.na(geneExpression)))){ #Checking if any models had problematic results
      cat("\n")
      message("NaN in expression data. Likely due to stiff equations. Try
               lowering step size to fix it.")
    }

    if(
      (configuration$simParams["printStart"] +
       configuration$simParams["printInterval"]) <
       configuration$simParams["simulationTime"]){

      timeStamps <- seq(
        configuration$simParams["printStart"] +
          configuration$simParams["printInterval"],
        configuration$simParams["simulationTime"],
        configuration$simParams["printInterval"])

      tsSimulations <- as.list(timeStamps)
      names(timeStamps) <- timeStamps

      for(i in seq_along(timeStamps)){
        tsSimulations[[i]] <- geneExpression[
          , (1+(i-1)*(nGenes)):(i*nGenes)]
        colnames(tsSimulations[[i]]) <- geneNames
      }

      geneExpression <- geneExpression[
        ,(1+(length(timeStamps)*nGenes)):
          ((length(timeStamps)+1)*nGenes)]

      tsSimulations <- lapply(tsSimulations,t)

      geneExpression <- t(geneExpression)
      rownames(geneExpression) <- geneNames
      assayDataTmp <- c(list(deterministic = geneExpression),
                        tsSimulations)
      metadataTmp$tsSimulations <- timeStamps


    }
    if((configuration$stochParams["nNoise"] == 0) &
       (configuration$simParams["printStart"] +
        configuration$simParams["printInterval"] >
       configuration$simParams["simulationTime"])){
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
      stochasticSimulations <- lapply(stochasticSimulations,t)
      if(simDet){
      geneExpression <- geneExpression[
        ,(1+(configuration$stochParams["nNoise"]*nGenes)):
          ((configuration$stochParams["nNoise"]+1)*nGenes)]


      geneExpression <- t(geneExpression)
      rownames(geneExpression) <- geneNames
      assayDataTmp <- c(list(deterministic = geneExpression),
                        stochasticSimulations)
      } else {
        assayDataTmp <- stochasticSimulations
      }
      metadataTmp$stochasticSimulations <- noiseLevels


    }

    paramName <- sracipeGenParamNames(rSet)
    # paramFile <- paste0("tmp/",outFile,"_parameters.txt")
    if(!(nCores>1 & !timeSeries)){
      parameters <- utils::read.table(outFileParams, header = FALSE)
    }
    colnames(parameters) <- paramName

    # icFile <- paste0("tmp/",outFile,"_IC.txt")
    if(!(nCores>1 & !timeSeries)){
      ic <- utils::read.table(outFileIC, header = FALSE)
    }
    colnames(ic) <- geneNames
    metadataTmp$normalized <- FALSE


    if(convergTesting) {
      if(!(nCores>1)){
        converge<- utils::read.table(outFileConverge, header = FALSE)
      }
      colnames(converge)<-c("Model Convergence", "Tests Done")

      if(nIC > 1){#Counts unique states per model
        geneExpressionRounded <- round(geneExpression, digits = uniqueDigits)
        uniqueStates <- numeric(numModels) #store number of unique states in a vector
        for(modelCount in seq_len(numModels)){
          #grab ICs and convergence data for each model
          startIdx <- (modelCount - 1)*nIC + 1
          endIdx <- modelCount*nIC

          finalModelExpressions <- geneExpressionRounded[, startIdx:endIdx]
          ICconvergences <- converge[startIdx:endIdx, 1]
          convergedICs <- finalModelExpressions[, as.logical(ICconvergences)]

          if(is.null(ncol(unique(convergedICs, MARGIN = 2)))){
            uniqueStates[modelCount] <- 0
          } else {
            uniqueStates[modelCount] <- ncol(unique(convergedICs, MARGIN = 2))
          }
        }
        StateCounts <- data.frame(modelNo = 1:numModels, UniqueStableStateNo = uniqueStates)
        metadataTmp$uniqueStateCounts <- StateCounts
      }

      if(limitcycles){ #Running limit cycle algorithm
        cat("\n")
        message("Checking for limit cycles")
        if(nCores>1){
          #Making sure outFileGE file has valid expressions
          utils::write.table(geneExpression, file = outFileGE,
                             sep = "\t", quote = FALSE, row.names = FALSE,
                             col.names = FALSE)
        }

        LCStepperInt <- 1L
        if(configuration$LCStepper == "RK4"){ LCStepperInt <- 4L}
        LC_Test <- limitcyclesGRC(geneInteraction, outFileLC, outFileLCIC, configuration, converge[,1],
                                  outFileParams, outFileGE, metadataTmp$geneTypes, LCStepperInt)
        if(LC_Test > 0){
          metadataTmp$totalNumofLCs <- LC_Test
          LCs <- utils::read.table(outFileLC, header = FALSE)
          colnames(LCs) <- c("Model No", "Limit Cycle No", "Period", geneNames)
          metadataTmp$LCData <- LCs
          LCICs <- utils::read.table(outFileLCIC, header = FALSE)
          #metadataTmp$LCICs <- LCICs
          #Recording the initial conditions which produced limit cycles and
          #categorizing them separately from non-converging models
          for(k in 1:nrow(LCICs)){
            modelNo <- LCICs[k, 1]
            positions <- which(LCICs[k, -1] == 1)
            startIdx <- (modelNo - 1)*nIC + 1
            converge[startIdx + positions - 1, 1] <- 2
          }

        }
        else{ #Reports errors in file handling
          message("No limit cycles detected")
        }
      }
      metadataTmp$modelConvergence <- converge
    }

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
          params <- parameters
          params[,knockOut_number] <- 0.0

          icTmp <- ic
          icTmp[,knockOut_number] <- 0.0

          configTmp <- configuration
          configTmp$options["genIC"] <- FALSE
          configTmp$options["genParams"] <- FALSE

          utils::write.table(params, file = outFileParams,
                      sep = "\t", quote = FALSE, row.names = FALSE,
                      col.names = FALSE)
          utils::write.table(icTmp, file = outFileIC,
                      sep = "\t", quote = FALSE, row.names = FALSE,
                      col.names = FALSE)

          Time_evolution_test<- simulateGRCCpp(geneInteraction, configTmp,
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



