#' @export
#' @title Load the sRACIPE configuration from a file.
#' @description Read the sRACIPE configuration file to initialize the parameters.
#' @param topologyFile Character. Topology file address.
#' @param header (optional) Logical. Whether the input file has a header or not.
#' \code{TRUE} by default.
#' @return List. Contains various parameters and their values.
#' @examples
#' \dontrun{
#' configuration <- loadConfig("inputs/sRACIPE.cfg")
#' }
#'
#' @section Related Functions:
#'
#' \code{\link{simulateGRC}},  \code{\link{knockdownAnalysis}},
#' \code{\link{overExprAnalysis}},  \code{\link{plotData}},
#' \code{\link{calcHeatmapSimilarity}}

loadConfig <- function(configFile="inputs/sRACIPE.cfg", header = TRUE){

  if(missing(configFile)){
    message("Configuration file not specified. Using the default values")
    configuration <- data.frame(ANNEALING=FALSE, NUM_MODELS=1000, PARAMETER_RANGE=100, MPR_MIN=1.,MPR_MAX=100, DNR_MIN=0.1,DNR_MAX=1.,FCH_MIN=1,FCH_MAX=100,HCO_MIN=1L,HCO_MAX=6L,STEP_SIZE=0.02,SIM_TIME=50.0,MEDIAN_RANGE=100.0,INITIAL_CONDITIONS=1L,OUTPUT_PRECISION=12L,THRESHOLD_MODELS=5000, RK_TOLERANCE = 2^(-2), NOISE_LEVELS=30, MAX_NOISE=30, NOISE_SCALING_FACTOR=0.9, standard_deviation_factor=0.5658,SHOT_NOISE_SCALING=0,GENE_NOISE_SCALING=0,FILE_WRITING_INTERVAL=1,CONSTANT_NOISE=1, PRINT_START = 50.0, PRINT_INTERVAL=10, stringsAsFactors = TRUE, PARAMETERS_FILE = 0L, READ_IC = 0L,possible_interactions = 3)

    return(configuration)
}

  if(file.exists(configFile)){

    f_inputs <- read.table(configFile,
                           header = header, stringsAsFactors = TRUE)
  }
  else
  {
    print("Please check the input filename. Make sure that full path is specified. For example, if it is in inputs folder, use inputs/filename as argument to the function. ")
  }
  data("configuration")

  configuration$ANNEALING <- as.logical(f_inputs[f_inputs[,1]=="ANNEALING",2])
  configuration$PARAMETER_RANGE <- as.numeric(f_inputs[f_inputs[,1]=="PARAMETER_RANGE",2])
  configuration$MPR_MIN <- as.numeric(f_inputs[f_inputs[,1]=="MPR_MIN",2])
  configuration$MPR_MAX <- as.numeric(f_inputs[f_inputs[,1]=="MPR_MAX",2])
  configuration$DNR_MIN <- as.numeric(f_inputs[f_inputs[,1]=="DNR_MIN",2])
  configuration$DNR_MAX <- as.numeric(f_inputs[f_inputs[,1]=="DNR_MAX",2])
  configuration$FCH_MIN <- as.numeric(f_inputs[f_inputs[,1]=="FCH_MIN",2])
  configuration$FCH_MAX <- as.numeric(f_inputs[f_inputs[,1]=="FCH_MAX",2])
  configuration$HCO_MIN <- as.numeric(f_inputs[f_inputs[,1]=="HCO_MIN",2])
  configuration$HCO_MAX <- as.numeric(f_inputs[f_inputs[,1]=="HCO_MAX",2])
  configuration$STEP_SIZE <- as.numeric(f_inputs[f_inputs[,1]=="EULER_STEP_SIZE",2])
  configuration$SIM_TIME <- as.numeric(f_inputs[f_inputs[,1]=="EULER_SIM_TIME",2])


  configuration$possible_interactions= as.integer(3);
  configuration$MEDIAN_RANGE <- as.numeric(f_inputs[f_inputs[,1]=="MEDIAN_RANGE",2])
  configuration$NUM_MODELS <- as.integer(f_inputs[f_inputs[,1]=="NUM_MODELS",2])
  configuration$THRESHOLD_MODELS <- as.integer(f_inputs[f_inputs[,1]=="THRESHOLD_MODELS",2])
  configurationTHRESHOLD_MODELS <- as.integer(max(configuration$NUM_MODELS,configuration$THRESHOLD_MODELS,5000))
  configuration$NOISE_LEVELS <- as.integer(f_inputs[f_inputs[,1]=="NOISE_LEVELS",2])
  configuration$MAX_NOISE <- as.numeric(f_inputs[f_inputs[,1]=="MAX_NOISE",2])
  configuration$NOISE_SCALING_FACTOR <- as.numeric(f_inputs[f_inputs[,1]=="NOISE_SCALING_FACTOR",2])
  configuration$standard_deviation_factor= as.numeric(0.5658);
  configuration$SHOT_NOISE_SCALING <- as.logical(f_inputs[f_inputs[,1]=="SHOT_NOISE_SCALING",2])
  configuration$GENE_NOISE_SCALING <- as.logical(f_inputs[f_inputs[,1]=="GENE_NOISE_SCALING",2])
  configuration$OUTPUT_PRECISION <- as.integer(f_inputs[f_inputs[,1]=="OUTPUT_PRECISION",2])
  configuration$INITIAL_CONDITIONS  <- as.integer(f_inputs[f_inputs[,1]=="INITIAL_CONDITIONS",2])
#    message("Configuration file successfully loaded.")
    return(configuration)
}

#' @export
#' @title  (deprecated) Load the sRACIPE configuration from a file.
#' @description Read the sRACIPE configuration file to initialize the parameters.
#' Use the function \code{\link{loadConfig}} instead.
#' @param topologyFile Character. Topology file address.
#' @param header (optional) Logical. Whether the input file has a header or not.
#' \code{TRUE} by default.
#' @return List. Contains various parameters and their values.
#' @examples
#' \dontrun{
#' configuration <- loadConfig("inputs/sRACIPE.cfg")
#' }
#'
sRACIPE_load_configuration <- function(config_file="inputs/sRACIPE.cfg"){
  loadConfig(config_file)
}
