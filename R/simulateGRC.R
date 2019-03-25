
#' @export
#' @title Simulate a gene regulatory circuit using RACIPE
#' @description Simulate a gene regulatory circuit using its topology as the
#' only input. It will generate an ensemble of random models.
#' @param circuit List or character. The file containing the circuit or the list
#' returned by the function \code{\link{loadCircuit}}
#' @param config  (optional) data.frame or character. The file containing the
#' configuration or the dataframe
#' returned by the function \code{\link{loadConfig}}
#' @param anneal (optional) Logical. Default FALSE. Whether to use annealing
#' for stochastic simulations. If TRUE, the gene expressions at higher noise
#' are used as initial conditions for simulations at lower noise.
#' @param numModels (optional) Integer. Default 500. Number of random models
#' to be simulated.
#' @param thresholdModels (optional)  integer. Default 10000. The number of
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
#' @param outputPrecison (optional) integer. Default 12.
#' The decimal point precison of
#' the output to be printed in the tmp folder.
#' @param knockOut (optional) List of character or vector of characters.
#' Simulation after knocking out one or more genes. To knock out all the genes
#' in the circuit, use \code{knockOut = "all"}. If it is a vector, then all
#' the genes in the vector will be knocked out simultaneously.
#' @param printStart (optional) numeric (0-\code{simulationTime}).
#'  Default \code{simulationTime}
#' The time from which the output should be recorded. Useful for time series
#' analysis and studying the dynamics of a model for a particular initial
#' condition.
#' @param printInterval (optional) numeric (\code{integrateStepSize}-
#' \code{simulationTime - printStart}). Default 10. The separation between
#' two recorded time points.
#' @param stepper (optional) Character. Stepper to be used for integrating the
#' differential equations. The options include \code{"EM"} for Euler-Maruyama
#' O(1), \code{"RK4"}
#' for fourth order Runge-Kutta O(4) and \code{"DP"} for adaptive stepper based
#' Dormand-Prince algorithm. The default method is \code{"RK4"} for deterministic
#' simulations and the method defaults to \code{"EM"} for stochastic simulations.
#' @param plots (optional) logical Default \code{TRUE}.
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
#'
#' @section Related Functions:
#'
#' \code{\link{simulateGRC}},  \code{\link{knockdownAnalysis}},
#' \code{\link{overExprAnalysis}},  \code{\link{plotData}},
#' \code{\link{calcHeatmapSimilarity}}

simulateGRC <- function( circuit="inputs/test.tpo", config ="inputs/sRACIPE.cfg",
                      anneal=FALSE, knockOut = NA_character_,numModels=500, paramRange=100,
                      prodRateMin=1.,prodRateMax=100, degRateMin=0.1,
                      degRateMax=1.,foldChangeMin=1,foldChangeMax=100,
                      hillCoeffMin=1L,hillCoeffMax=6L,integrateStepSize=0.02,
                      simulationTime=50.0,nIC=1L,nNoise=0L,
                      initialNoise=50.0,noiseScalingFactor=0.5,shotNoise=0,
                      scaledNoise=FALSE,outputPrecision=12L,
                       printStart = 50.0,
                      printInterval=10, stepper = "RK4",
                      thresholdModels = 10000, plots = FALSE, plotToFile = TRUE,
                      genIC = TRUE, genParams = TRUE,
                      integrate = TRUE, rkTolerance = 0.01, ...){
 rSet <- list()
 topology <- list()
  data("configuration", envir=environment(), package = "sRACIPE")
   if((class(circuit) == "character") | (class(circuit) == "data.frame"))
  {
    topology <- loadCircuit(circuit = circuit)
    if(class(circuit) == "data.frame"){
      circuit <-  topology$topology_filepath
    }
  }
  if(class(circuit) == "list"){
    topology <- circuit
  }
 if(missing(circuit)){
    message("Please specify a circuit either as a file or list created by loadTopology function.")
    return()
  }
if(!missing(config)){
 if(class(config) == "character"){
   configuration <- loadConfig(config)
 }
 if(class(config) == "list"){
   configuration <- config
 }
}

 # if(missing(config)){
 #   configuration <-  #readRDS("data/config.RDS")
 # }

 if(!missing(knockOut)){
   if(!(class(knockOut) == "list")){
     knockOut <- list(knockOut)
   }
 }

  if(!missing(anneal)){
    #print("test1")
 if(anneal)
      configuration$anneal <- anneal
  }
  if(!missing(numModels)){
    configuration$numModels <- numModels
  }

  if(!missing(paramRange)){
    configuration$paramRange <- paramRange
  }
  if(!missing(prodRateMin)){
    configuration$prodRateMin <- prodRateMin
  }
  if(!missing(prodRateMax)){
    configuration$prodRateMax <- prodRateMax
  }

  if(!missing(degRateMin)){
    configuration$degRateMin <- degRateMin
  }
  if(!missing(degRateMax)){
    configuration$degRateMax <- degRateMax
  }
  if(!missing(foldChangeMin)){
    configuration$foldChangeMin <- foldChangeMin
  }
  if(!missing(foldChangeMax)){
    configuration$foldChangeMax <- foldChangeMax
  }
  if(!missing(hillCoeffMin)){
    configuration$hillCoeffMin <- hillCoeffMin
  }
  if(!missing(hillCoeffMax)){
    configuration$hillCoeffMax <- hillCoeffMax
  }
  if(!missing(integrateStepSize)){
    configuration$integrateStepSize <- integrateStepSize
  }
  if(!missing(simulationTime)){
    configuration$simulationTime <- simulationTime
  }
  if(!missing(nIC)){
    configuration$nIC <- nIC
  }
  if(!missing(outputPrecision)){
    configuration$outputPrecision <- outputPrecision
  }

  if(!missing(nNoise)){
    configuration$nNoise <- nNoise
  }
  if(!missing(initialNoise)){
    configuration$initialNoise <- initialNoise
  }
  if(!missing(noiseScalingFactor)){
    configuration$noiseScalingFactor <- noiseScalingFactor
  }
  if(!missing(shotNoise)){
    configuration$shotNoise <- shotNoise
  }
  if(!missing(scaledNoise)){
    configuration$scaledNoise <-scaledNoise
  }
  if(!missing(printStart)){
    configuration$printStart <- printStart
  }
  if(!missing(printInterval)){
    configuration$printInterval <- printInterval
  }

 if(!missing(genIC)){
   configuration$genIC <- genIC
 }

 if(!missing(genParams)){
   configuration$genParams <- genParams
 }

 if(!missing(integrate)){
   configuration$integrate <- integrate
 }
 if(missing(printStart)){
  configuration$printStart <- configuration$simulationTime
 }
 if(!missing(rkTolerance)){
   configuration$rkTolerance <- rkTolerance
 }
 # Apply parameter range
  configuration$prodRateMin <- 0.5*(
    configuration$prodRateMin + configuration$prodRateMax) - 0.5*(
      configuration$prodRateMax - configuration$prodRateMin)*
    configuration$paramRange/100
  configuration$prodRateMax <- 0.5*(
    configuration$prodRateMin + configuration$prodRateMax) + 0.5*(
      configuration$prodRateMax - configuration$prodRateMin)*
    configuration$paramRange/100

  configuration$degRateMin <- 0.5*(
    configuration$degRateMin + configuration$degRateMax) - 0.5*(
      configuration$degRateMax - configuration$degRateMin)*
    configuration$paramRange/100
  configuration$degRateMax <- 0.5*(
    configuration$degRateMin + configuration$degRateMax) + 0.5*(
      configuration$degRateMax - configuration$degRateMin)*
    configuration$paramRange/100

  configuration$foldChangeMin <- 0.5*(
    configuration$foldChangeMin + configuration$foldChangeMax) - 0.5*(
      configuration$foldChangeMax - configuration$foldChangeMin)*
    configuration$paramRange/100
  configuration$foldChangeMax <- 0.5*(
    configuration$foldChangeMin + configuration$foldChangeMax) + 0.5*(
      configuration$foldChangeMax - configuration$foldChangeMin)*
    configuration$paramRange/100

  # configuration$initialNoise <- configuration$initialNoise/(topology$number_gene)
#  print(configuration)


  results_directory <- ifelse(!dir.exists(file.path(getwd(), "results")),
                              dir.create(file.path(getwd(), "results")), TRUE)
  results_directory <- ifelse(!dir.exists(file.path(getwd(), "tmp")),
                              dir.create(file.path(getwd(), "tmp")), TRUE)
  output_directory <- file.path(getwd(), "results")


  geneInteraction <- matrix(0, nrow = topology$number_gene,
                            ncol = topology$number_gene)

  thresholdGene <- rep(0, topology$number_gene)
  #Rcpp::sourceCpp("src/interaction_reader.cpp")
  geneNames <- character(length = topology$number_gene)


  geneInteraction <- readTopology(geneInteraction,topology$topology_filepath,
                                   topology$filename, geneNames)
  geneNames <- gsub("[[:punct:]]", "", geneNames)

  rSet$geneNames <- geneNames
  outFile <- paste0(Sys.Date(),"_",topology$filename,"_", basename(tempfile()))
  if(genParams){
  message("Generating gene thresholds")
  #Rcpp::sourceCpp("src/thresholdGenerator.cpp")
  threshold_test <- generateThresholds(
    geneInteraction,  thresholdGene,  configuration$prodRateMin,
    configuration$prodRateMax, configuration$degRateMin,
    configuration$degRateMax,  configuration$interactionTypes,
    configuration$numModels,  configuration$thresholdModels,
    configuration$integrateStepSize,  configuration$foldChangeMin,
    configuration$foldChangeMax,  configuration$hillCoeffMin,
    configuration$hillCoeffMax,   configuration$sdFactor)
  if(threshold_test!=0){
    print("Error in threshold generation")
    }
  } else {

    if((!class(circuit)=="list")){
    message("Please specify the parameters as circuit$params as genParams is FALSE")
    return(rSet)
  } else {
    params <- circuit$params
    write.table(params, file = paste0("tmp/",outFile,"_parameters.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  }

    if(!genIC){
      if((!class(circuit)=="list")){
      message("Please specify the initial conditions
              as circuit$ic as genIC is FALSE")
      return(rSet)
    } else {
      ic <- circuit$ic
      write.table(ic, file = paste0("tmp/",outFile,"_IC.txt"),
                  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  }
if(missing(nNoise)){
  if(!missing(initialNoise) & !anneal)
    configuration$nNoise <- 2
}
  if(anneal){
    if(missing(nNoise)) {configuration$nNoise <- 30L}

  }

  if(configuration$nNoise > 0) {
    if(stepper != "EM"){
      warnings("Defaulting to EM stepper for stochastic simulations")
      stepper <- "EM"
    }
    if(missing(initialNoise)) { configuration$initialNoise <- 50/sqrt(topology$number_gene)}
  }

  stepperInt <- 1L
  if(stepper == "RK4"){ stepperInt <- 4L}
  if(stepper == "DP") {stepperInt <- 5L}



configuration$stepper <- stepper
  rSet$fileName <- outFile
  rSet$topology <- topology
  rSet$config <- configuration

  message("Running the simulations")

    Time_evolution_test<- simulateGRCCpp(
    gene_interaction =  geneInteraction, threshold_gene =   thresholdGene,
    g_min =  configuration$prodRateMin,  g_max =  configuration$prodRateMax,
    k_min = configuration$degRateMin,  k_max =  configuration$degRateMax,
    interaction_types =  configuration$interactionTypes,
    model_count_max =  configuration$numModels,
    threshold_max =  configuration$thresholdModels,
    h =  configuration$integrateStepSize,
    lambda_min =  configuration$foldChangeMin,
    lambda_max =  configuration$foldChangeMax,
    n_min =  configuration$hillCoeffMin,
    n_max =   configuration$hillCoeffMax,
    tot_time =  configuration$simulationTime,
    sd_multiplier =  configuration$sdFactor,
    number_gene =  topology$number_gene,
    D_max =   configuration$initialNoise,
    D_shot_scaling =  configuration$shotNoise,
    scaled_noise =  configuration$scaledNoise,
    D_levels =  (1+configuration$nNoise),
    D_scaling =  configuration$noiseScalingFactor,
    output_precision =  configuration$outputPrecision,
    ANNEALING =   configuration$anneal,
    initial_conditions =   configuration$nIC,
    filename =   outFile,
    print_start = configuration$printStart,
    print_interval = configuration$printInterval,
    integrate = configuration$integrate,
    genParams = configuration$genParams,
    genIC = configuration$genIC,
    stepper = stepperInt, rk_tolerance = configuration$rkTolerance)

    if(integrate){
    geFile <- paste0("tmp/",outFile,"_geneExpression.txt")
    geneExpression <- read.table(geFile, header = FALSE)
    colnames(geneExpression) <- geneNames
    rSet$geneExpression <- geneExpression[
      ,(1+ncol(geneExpression) -length(geneNames)):ncol(geneExpression)]
    names(rSet$geneExpression) <- geneNames
    }
    if(configuration$nNoise > 0){
      noiseLevels <- (initialNoise*noiseScalingFactor^seq(0,configuration$nNoise-1))
      stochasticSimulations <- as.list(noiseLevels)
      names(stochasticSimulations) <- noiseLevels
      for(i in 1:configuration$nNoise){
        stochasticSimulations[[i]] <- geneExpression[
          , (1+(i-1)*(topology$number_gene)):(i*topology$number_gene)]
        colnames(stochasticSimulations[[i]]) <- geneNames
      }
      rSet$stochasticSimulations <- stochasticSimulations
      geneExpression <- geneExpression[
        ,(1+(configuration$nNoise)*(topology$number_gene)):((configuration$nNoise+1)*topology$number_gene)]
    }


    paramName <- genParamNames(circuit)
    paramFile <- paste0("tmp/",outFile,"_parameters.txt")
    parameters <- read.table(paramFile, header = FALSE)
    colnames(parameters) <- paramName

    icFile <- paste0("tmp/",outFile,"_IC.txt")
    ic <- read.table(icFile, header = FALSE)
    colnames(ic) <- geneNames


    rSet$normalized <- FALSE
    rSet$params <- parameters
    rSet$ic <- ic

  ##################
  # Handle knockouts here
      if(!missing(knockOut)){
        if(knockOut[[1]] == "all"){
          knockOut <- as.list(geneNames)
        }
        knockOutData <- knockOut
        names(knockOutData) <- knockOut
        for(ko in 1:length(knockOut)){
          koGene <- knockOut[[ko]]
          knockOut_number <- which(koGene==geneNames)
          if(length(knockOut_number)==0){
            message("knockOut gene not found in the circuit")
            return(rSet)
          }
          params <- rSet$params
          ic <- rSet$ic
          ic[,knockOut_number] <- 0
          params[,knockOut_number] <- 0
          outFileKO <- paste0(outFile, "_KO_",paste(knockOut[[ko]], sep="", collapse="_"))
          write.table(params, file = paste0("tmp/",outFileKO,"_parameters.txt"),
                      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
          write.table(ic, file = paste0("tmp/",outFileKO,"_IC.txt"),
                      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

          Time_evolution_test<- simulateGRCCpp(
            gene_interaction =  geneInteraction, threshold_gene =   thresholdGene,
            g_min =  configuration$prodRateMin,  g_max =  configuration$prodRateMax,
            k_min = configuration$degRateMin,  k_max =  configuration$degRateMax,
            interaction_types =  configuration$interactionTypes,
            model_count_max =  configuration$numModels,
            threshold_max =  configuration$thresholdModels,
            h =  configuration$integrateStepSize,
            lambda_min =  configuration$foldChangeMin,
            lambda_max =  configuration$foldChangeMax,
            n_min =  configuration$hillCoeffMin,
            n_max =   configuration$hillCoeffMax,
            tot_time =  configuration$simulationTime,
            sd_multiplier =  configuration$sdFactor,
            number_gene =  topology$number_gene,
            D_max =   configuration$initialNoise,
            D_shot_scaling =  configuration$shotNoise,
            scaled_noise =  configuration$scaledNoise,
            D_levels =  (1+configuration$nNoise),
            D_scaling =  configuration$noiseScalingFactor,
            output_precision =  configuration$outputPrecision,
            ANNEALING =   configuration$anneal,
            initial_conditions =   configuration$nIC,
            filename =   outFileKO,
            print_start = configuration$printStart,
            print_interval = configuration$printInterval,
            integrate = configuration$integrate,
            genParams = FALSE,
            genIC = FALSE,
            stepper = stepperInt, rk_tolerance = configuration$rkTolerance)


          geFile <- paste0("tmp/",outFileKO,"_geneExpression.txt")
          geneExpression <- read.table(geFile, header = FALSE)
          colnames(geneExpression) <- geneNames
          knockOutData[[ko]] <- geneExpression
        }
        rSet$knockOutSimulations <- knockOutData
  }

##############################
 if(plots){
    rSet <- plotData(rSet,plotToFile = plotToFile,...)
  }

    saveRDS(rSet, file = paste0("results/",outFile,".RDS"))
return(rSet)
}



#' @export
#' @title Generate parameter names for a circuit
#'
#'
genParamNames <- function(circuit="inputs/test.tpo"){
  topology <- list()
  configuration <- list()
  if(class(circuit) == "character")
  {
    topology <- loadCircuit(circuit = circuit)
  }
  if(class(circuit) == "list"){
    topology <- circuit
  }
  if(missing(circuit)){
    message("Please specify a circuit either as a file or list created by loadTopology function.")
    return()
  }

  geneInteraction <- matrix(0, nrow = topology$number_gene, ncol = topology$number_gene)

  geneNames <- character(length = topology$number_gene)

  geneInteraction <- readTopology(geneInteraction,topology$topology_filepath,
                                  topology$filename, geneNames)

  paramList <- list()
  tmp <- lapply(geneNames,  function(x) paste("G_",x, sep=""))
  paramList <- append(paramList, tmp)
  tmp <- lapply(geneNames,  function(x) paste("K_",x, sep=""))
  paramList <- append(paramList, tmp)
  tmp <- list()
  for(gene1 in 1:length(geneNames))
  {
    for(gene2 in 1:length(geneNames))
    {
      if(geneInteraction[gene1,gene2]>0)
        tmp <- append(tmp, paste("TH_",geneNames[[gene2]],"_",geneNames[[gene1]], sep=""))
    }
  }
  paramList <- append(paramList, tmp)

  tmp <- list()
  for(gene1 in 1:length(geneNames))
  {
    for(gene2 in 1:length(geneNames))
    {
      if(geneInteraction[gene1,gene2]>0)
        tmp <- append(tmp, paste("N_",geneNames[[gene2]],"_",geneNames[[gene1]],sep=""))
    }
  }
  paramList <- append(paramList, tmp)

  tmp <- list()
  for(gene1 in 1:length(geneNames))
  {
    for(gene2 in 1:length(geneNames))
    {
      if(geneInteraction[gene1,gene2]>0)
        tmp <- append(tmp, paste("FC_",geneNames[[gene2]],"_",geneNames[[gene1]],sep=""))
    }
  }

  paramList <- append(paramList, tmp)
  return(as.character(paramList))

}

sRACIPE_simulate_GRN <- function(
  topology_file="inputs/test.tpo", config_file ="inputs/sRACIPE.cfg",ANNEAL=F,
  NUM_MODELS=100, PARAMETER_RANGE=100, MPR_MIN=1.,MPR_MAX=100, DNR_MIN=0.1,
  DNR_MAX=1.,FCH_MIN=1,FCH_MAX=100,HCO_MIN=1L,hillCoeffMax=6L,STEP_SIZE=0.02,
  SIM_TIME=50.0,MEDIAN_RANGE=100.0,INITIAL_CONDITIONS=1L,NOISE_LEVELS=30L,
  MAX_NOISE=30.0,NOISE_SCALING_FACTOR=0.9,SHOT_NOISE_SCALING=0L,
  GENE_NOISE_SCALING=0L,FILE_WRITING_INTERVAL=1L,OUTPUT_PRECISION=12L,
  PARAMETERS_FILE = 0L, READ_IC = 0L,KNOCKOUT = NA_character_,
  PRINT_START = 50.0, PRINT_INTERVAL=10){

  simulateGRC(circuit = topology_file, config = config_file,
              anneal = ANNEAL, numModels=NUM_MODELS, paramRange=PARAMETER_RANGE,
           prodRateMin=MPR_MIN,prodRateMax=MPR_MAX, degRateMin=DNR_MIN,
           degRateMax=DNR_MAX,foldChangeMin=FCH_MIN,foldChangeMax=FCH_MAX,
           hillCoeffMin=HCO_MIN,hillCoeffMax=HCO_MAX,integrateStepSize=STEP_SIZE,
           simulationTime=SIM_TIME,nIC=INITIAL_CONDITIONS,nNoise=NOISE_LEVELS,
           initialNoise=MAX_NOISE,noiseScalingFactor=NOISE_SCALING_FACTOR,
           shotNoise=SHOT_NOISE_SCALING,scaledNoise=GENE_NOISE_SCALING,
           outputPrecision=OUTPUT_PRECISION, knockOut = NA_character_,
           printStart = PRINT_START, printInterval=PRINT_INTERVAL)

}

