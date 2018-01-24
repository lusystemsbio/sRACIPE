
#' Function to generate deterministic trajectories


sRACIPE_RK_deterministic <- function( NUM_MODELS=100, PARAMETER_RANGE=100, MPR_MIN=1.,MPR_MAX=100, DNR_MIN=0.1,DNR_MAX=1.,FCH_MIN=1,FCH_MAX=100,HCO_MIN=1L,HCO_MAX=6L,EULER_STEP_SIZE=0.02,EULER_SIM_TIME=50.0,MEDIAN_RANGE=100.0,INITIAL_CONDITIONS=1L,OUTPUT_PRECISION=12L){

  working_directory<-getwd()
  if (dir.exists(file.path(working_directory, "inputs"))){
    #print("Reading configuration file from inputs folder")
    input_directory<-file.path(working_directory, "inputs")
    setwd(input_directory)

    f_inputs <- read.table(dir(pattern = ".cfg"),
                           header = TRUE)
    setwd(working_directory)

    if(missing(NUM_MODELS)){
      NUM_MODELS <- f_inputs[f_inputs[,1]=="NUM_MODELS",2]
      print(NUM_MODELS)
      print(f_inputs[f_inputs[,1]=="NUM_MODELS",2])
    }
    if(missing(PARAMETER_RANGE)){
      PARAMETER_RANGE <-f_inputs[f_inputs[,1]=="PARAMETER_RANGE",2]
    }
    if(missing(MPR_MIN)){
      MPR_MIN <- f_inputs[f_inputs[,1]=="MPR_MIN",2]
    }
    if(missing(MPR_MAX)){
      MPR_MAX <-f_inputs[f_inputs[,1]=="MPR_MAX",2]
    }

    if(missing(DNR_MIN)){
      DNR_MIN <-f_inputs[f_inputs[,1]=="DNR_MIN",2]
    }
    if(missing(DNR_MAX)){
      DNR_MAX <-f_inputs[f_inputs[,1]=="DNR_MAX",2]
    }
    if(missing(FCH_MIN)){
      FCH_MIN <-f_inputs[f_inputs[,1]=="FCH_MIN",2]
    }
    if(missing(DNR_MAX)){
      FCH_MAX <-f_inputs[f_inputs[,1]=="FCH_MAX",2]
    }
    if(missing(HCO_MIN)){
      HCO_MIN <-f_inputs[f_inputs[,1]=="HCO_MIN",2]
    }
    if(missing(HCO_MAX)){
      HCO_MAX <-f_inputs[f_inputs[,1]=="HCO_MAX",2]
    }
    if(missing(EULER_STEP_SIZE)){
      EULER_STEP_SIZE <-f_inputs[f_inputs[,1]=="EULER_STEP_SIZE",2]
    }
    if(missing(EULER_SIM_TIME)){
      EULER_SIM_TIME <-f_inputs[f_inputs[,1]=="EULER_SIM_TIME",2]
    }
    if(missing(MEDIAN_RANGE)){
      MEDIAN_RANGE <-f_inputs[f_inputs[,1]=="MEDIAN_RANGE",2]
    }
    if(missing(INITIAL_CONDITIONS)){
      INITIAL_CONDITIONS <-f_inputs[f_inputs[,1]=="INITIAL_CONDITIONS",2]
    }
    if(missing(OUTPUT_PRECISION)){
      OUTPUT_PRECISION <-f_inputs[f_inputs[,1]=="OUTPUT_PRECISION",2]
    }

    print("Configuration file successfully loaded.")
    possible_interactions<<-3;
    standard_deviation_factor<<-0.5658;
    THRESHOLD_MODELS <<-5000
    THRESHOLD_MODELS <<- max(NUM_MODELS,THRESHOLD_MODELS)

  } else {
    print("Missing inputs folder! Please make sure that you are in right working directory and the directory contains an input folder with a .cfg file. Using default cfg inputs")}

  if (dir.exists(file.path(working_directory, "inputs"))){
    print("Reading topology file from inputs folder")
    input_directory<-file.path(working_directory, "inputs")
    setwd(input_directory)

    tpo_filename<<-dir(pattern = ".tpo")
    f_tpo <<- read.table(tpo_filename, header = TRUE)
    f_tpo

    setwd(working_directory)
    number_gene <<-length(levels(f_tpo$Source))

  } else {
    print("Missing inputs folder! Please make sure that you are in right working directory and the directory contains an input folder with a .tpo file. Using topology of a toggle switch with one self activating gene. You can directly modify f_tpo directly as well.")}

  MPR_MIN <- 0.5*(MPR_MIN+MPR_MAX) - 0.5*(MPR_MAX - MPR_MIN)*PARAMETER_RANGE/100
  MPR_MAX <- 0.5*(MPR_MIN+MPR_MAX) + 0.5*(MPR_MAX - MPR_MIN)*PARAMETER_RANGE/100

  DNR_MIN <- 0.5*(DNR_MIN+DNR_MAX) - 0.5*(DNR_MAX - DNR_MIN)*PARAMETER_RANGE/100
  DNR_MAX <- 0.5*(DNR_MIN+DNR_MAX) + 0.5*(DNR_MAX - DNR_MIN)*PARAMETER_RANGE/100

  FCH_MIN <- 0.5*(FCH_MIN+FCH_MAX) - 0.5*(FCH_MAX - FCH_MIN)*PARAMETER_RANGE/100
  FCH_MAX <- 0.5*(FCH_MIN+FCH_MAX) + 0.5*(FCH_MAX - FCH_MIN)*PARAMETER_RANGE/100


  results_directory<-ifelse(!dir.exists(file.path(working_directory, "results")), dir.create(file.path(working_directory, "results")), TRUE)
  output_directory<-file.path(working_directory, "results")

  gene_interaction <<- matrix(0, nrow = number_gene, ncol = number_gene)
  threshold_gene <<- rep(0, number_gene)
  Rcpp::sourceCpp("interaction_reader.cpp")
  gene_interaction <<- interaction_reader(gene_interaction, tpo_filename, number_gene)



  print("Generating gene thresholds")
  Rcpp::sourceCpp("threshold_generator.cpp")

  threshold_test<- threshold_calculator_uniform( gene_interaction,  threshold_gene,  MPR_MIN,  MPR_MAX, DNR_MIN,  DNR_MAX,  possible_interactions,   NUM_MODELS,  THRESHOLD_MODELS,  EULER_STEP_SIZE,  FCH_MIN, FCH_MAX,  HCO_MIN,  HCO_MAX,  MEDIAN_RANGE,  standard_deviation_factor)
  if(threshold_test!=0)
    print("Error in threshold generation")


  print("Running the determinsitic simulations using Runge-Kutta")

  Rcpp::sourceCpp("multiGeneCircuit_RK_deterministic.cpp")

  Time_evolution_test<- multiGeneCircuit_RK_deterministic( gene_interaction,  threshold_gene,  MPR_MIN,  MPR_MAX,DNR_MIN,  DNR_MAX,possible_interactions,NUM_MODELS, THRESHOLD_MODELS,  EULER_STEP_SIZE,  FCH_MIN,FCH_MAX, HCO_MIN,  HCO_MAX,  EULER_SIM_TIME,  MEDIAN_RANGE, standard_deviation_factor, number_gene,  OUTPUT_PRECISION, INITIAL_CONDITIONS)


  if(Time_evolution_test!=0)
    print("Error in time evolution function")


}
