sRACIPE_load_configuration <- function(config_file="inputs/sRACIPE.cfg"){

  if(missing(config_file)){
    message("Configuration file not specified. Using the default values")
    configuration <- data.frame(ANNEALING=FALSE, possible_interactions = 3, NUM_MODELS=1000, PARAMETER_RANGE=100, MPR_MIN=1.,MPR_MAX=100, DNR_MIN=0.1,DNR_MAX=1.,FCH_MIN=1,FCH_MAX=100,HCO_MIN=1L,HCO_MAX=6L,STEP_SIZE=0.02,SIM_TIME=50.0,MEDIAN_RANGE=100.0,INITIAL_CONDITIONS=1L,OUTPUT_PRECISION=12L,THRESHOLD_MODELS=5000, RK_TOLERANCE = 2^(-2), NOISE_LEVELS=2, MAX_NOISE=30, NOISE_SCALING_FACTOR=0.9, standard_deviation_factor=0.5658,SHOT_NOISE_SCALING=0,GENE_NOISE_SCALING=0,FILE_WRITING_INTERVAL=1,CONSTANT_NOISE=1, stringsAsFactors = T)

    return()
}
  else{
    cfg_filename <- basename(config_file)
  }
  if(file.exists(config_file)){

    f_inputs <- read.table(config_file,
                           header = TRUE, stringsAsFactors = T)
  }
  else
  {
    print("Please check the input filename. Make sure that full path is specified. For example, if it is in inputs folder, use inputs/filename as argument to the function. ")
  }
  configuration <- data.frame(ANNEALING=FALSE, possible_interactions = 3, NUM_MODELS=1000, PARAMETER_RANGE=100, MPR_MIN=1.,MPR_MAX=100, DNR_MIN=0.1,DNR_MAX=1.,FCH_MIN=1,FCH_MAX=100,HCO_MIN=1L,HCO_MAX=6L,STEP_SIZE=0.02,SIM_TIME=50.0,MEDIAN_RANGE=100.0,INITIAL_CONDITIONS=1L,OUTPUT_PRECISION=12L,THRESHOLD_MODELS=5000, RK_TOLERANCE = 2^(-2), NOISE_LEVELS=2, MAX_NOISE=30, NOISE_SCALING_FACTOR=0.9, standard_deviation_factor=0.5658,SHOT_NOISE_SCALING=0,GENE_NOISE_SCALING=0,FILE_WRITING_INTERVAL=1,CONSTANT_NOISE=1, stringsAsFactors = T)

  configuration$ANNEALING <- f_inputs[f_inputs[,1]=="ANNEALING",2]
  configuration$PARAMETER_RANGE <- f_inputs[f_inputs[,1]=="PARAMETER_RANGE",2]
  configuration$MPR_MIN <- f_inputs[f_inputs[,1]=="MPR_MIN",2]
  configuration$MPR_MAX <- f_inputs[f_inputs[,1]=="MPR_MAX",2]
  configuration$DNR_MIN <- f_inputs[f_inputs[,1]=="DNR_MIN",2]
  configuration$DNR_MAX <- f_inputs[f_inputs[,1]=="DNR_MAX",2]
  configuration$FCH_MIN <- f_inputs[f_inputs[,1]=="FCH_MIN",2]
  configuration$FCH_MAX <- f_inputs[f_inputs[,1]=="FCH_MAX",2]
  configuration$HCO_MIN <- f_inputs[f_inputs[,1]=="HCO_MIN",2]
  configuration$HCO_MAX <- f_inputs[f_inputs[,1]=="HCO_MAX",2]
  configuration$STEP_SIZE <- f_inputs[f_inputs[,1]=="EULER_STEP_SIZE",2]
  configuration$SIM_TIME <- f_inputs[f_inputs[,1]=="EULER_SIM_TIME",2]


  configuration$possible_interactions=3;
  configuration$MEDIAN_RANGE <- f_inputs[f_inputs[,1]=="MEDIAN_RANGE",2]
  configuration$NUM_MODELS <- f_inputs[f_inputs[,1]=="NUM_MODELS",2]
  configuration$THRESHOLD_MODELS <- f_inputs[f_inputs[,1]=="THRESHOLD_MODELS",2]
  configurationTHRESHOLD_MODELS <- max(configuration$NUM_MODELS,configuration$THRESHOLD_MODELS,5000)
  configuration$NOISE_LEVELS <- f_inputs[f_inputs[,1]=="NOISE_LEVELS",2]
  configuration$MAX_NOISE <- f_inputs[f_inputs[,1]=="MAX_NOISE",2]
  configuration$NOISE_SCALING_FACTOR <-f_inputs[f_inputs[,1]=="NOISE_SCALING_FACTOR",2]
  configuration$standard_deviation_factor=0.5658;
  configuration$SHOT_NOISE_SCALING <-f_inputs[f_inputs[,1]=="SHOT_NOISE_SCALING",2]
  configuration$GENE_NOISE_SCALING <- f_inputs[f_inputs[,1]=="GENE_NOISE_SCALING",2]
  configuration$FILE_WRITING_INTERVAL <-f_inputs[f_inputs[,1]=="FILE_WRITING_INTERVAL",2]
  configuration$OUTPUT_PRECISION <- f_inputs[f_inputs[,1]=="OUTPUT_PRECISION",2]
  configuration$CONSTANT_NOISE  <- 1
  configuration$INITIAL_CONDITIONS  <- f_inputs[f_inputs[,1]=="INITIAL_CONDITIONS",2]

    message("Configuration file successfully loaded.")
    return(configuration)
}
