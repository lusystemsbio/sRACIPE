sRACIPE_load_configuration<-function(){

  working_directory<-getwd()
  if (dir.exists(file.path(working_directory, "inputs"))){
    print("Reading configuration file from inputs folder")
    input_directory<-file.path(working_directory, "inputs")
    setwd(input_directory)

    f_inputs <- read.table(dir(pattern = ".cfg"),
                           header = TRUE)
    ANNEALING <<-f_inputs[f_inputs[,1]=="ANNEALING",2]

    PARAMETER_RANGE <<-f_inputs[f_inputs[,1]=="PARAMETER_RANGE",2]
    MPR_MIN<<-f_inputs[f_inputs[,1]=="MPR_MIN",2]
    MPR_MAX<<-f_inputs[f_inputs[,1]=="MPR_MAX",2]
    DNR_MIN<<-f_inputs[f_inputs[,1]=="DNR_MIN",2]
    DNR_MAX<<-f_inputs[f_inputs[,1]=="DNR_MAX",2]
    FCH_MIN<<-f_inputs[f_inputs[,1]=="FCH_MIN",2]
    FCH_MAX<<-f_inputs[f_inputs[,1]=="FCH_MAX",2]
    HCO_MIN<<-f_inputs[f_inputs[,1]=="HCO_MIN",2]
    HCO_MAX<<-f_inputs[f_inputs[,1]=="HCO_MAX",2]
    EULER_STEP_SIZE<<-f_inputs[f_inputs[,1]=="EULER_STEP_SIZE",2]
    EULER_SIM_TIME <<-f_inputs[f_inputs[,1]=="EULER_SIM_TIME",2]


    possible_interactions=3;
    MEDIAN_RANGE <<-f_inputs[f_inputs[,1]=="MEDIAN_RANGE",2]
    NUM_MODELS <<-f_inputs[f_inputs[,1]=="NUM_MODELS",2]
    THRESHOLD_MODELS <<-f_inputs[f_inputs[,1]=="THRESHOLD_MODELS",2]
    THRESHOLD_MODELS=max(NUM_MODELS,THRESHOLD_MODELS)
    NOISE_LEVELS <<- f_inputs[f_inputs[,1]=="NOISE_LEVELS",2]
    MAX_NOISE <<- f_inputs[f_inputs[,1]=="MAX_NOISE",2]
    NOISE_SCALING_FACTOR <<-f_inputs[f_inputs[,1]=="NOISE_SCALING_FACTOR",2]
    standard_deviation_factor=0.5658;
    SHOT_NOISE_SCALING <<-f_inputs[f_inputs[,1]=="SHOT_NOISE_SCALING",2]
    GENE_NOISE_SCALING <<- f_inputs[f_inputs[,1]=="GENE_NOISE_SCALING",2]
    FILE_WRITING_INTERVAL <<-f_inputs[f_inputs[,1]=="FILE_WRITING_INTERVAL",2]
    OUTPUT_PRECISION <<-f_inputs[f_inputs[,1]=="OUTPUT_PRECISION",2]
    #CONSTANT_NOISE<-f_inputs[f_inputs[,1]=="CONSTANT_NOISE",2]
    CONSTANT_NOISE <<-1
    INITIAL_CONDITIONS <<-f_inputs[f_inputs[,1]=="INITIAL_CONDITIONS",2]
    print("Configuration file successfully loaded.")
    setwd(working_directory)

  } else {
    print("Missing inputs folder! Please make sure that you are in right working directory and the directory contains an input folder with a .cfg file. Using default cfg inputs")}


}
