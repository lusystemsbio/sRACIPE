
#' Function to generate deterministic trajectories

sRACIPE_RK_deterministic_knockout <- function(  topology_file="inputs/test.tpo", config_file ="inputs/sRACIPE.cfg", NUM_MODELS=1000, PARAMETER_RANGE=100, MPR_MIN=1.,MPR_MAX=100, DNR_MIN=0.1,DNR_MAX=1.,FCH_MIN=1,FCH_MAX=100,HCO_MIN=1L,HCO_MAX=6L,STEP_SIZE=0.02,SIM_TIME=50.0,MEDIAN_RANGE=100.0,INITIAL_CONDITIONS=1L,OUTPUT_PRECISION=12L,THRESHOLD_MODELS=5000, KNOCKOUT = NA_character_){
  working_directory<-getwd()
  if(missing(topology_file)){
    if(!exists("topology")){
      message("Please specify a topology file first")
      return()
    }
  }
  else
  {
    topology <- sRACIPE_load_topology(topology_file = topology_file)
    message("Topology file loaded.")
    #print(configuration)
  }
  if(missing(KNOCKOUT)){
    print("Missing knockout")
    return()
  }


  if(missing(config_file)){
    if(!exists("configuration")){
      #configuration <- sRACIPE_load_configuration()
      message("Configuration file not found. Using default configuration settings")
      if(!exists("configuration")) configuration <- data.frame(NUM_MODELS=1000, PARAMETER_RANGE=100, MPR_MIN=1.,MPR_MAX=100, DNR_MIN=0.1,DNR_MAX=1.,FCH_MIN=1,FCH_MAX=100,HCO_MIN=1L,HCO_MAX=6L,STEP_SIZE=0.02,SIM_TIME=50.0,MEDIAN_RANGE=100.0,INITIAL_CONDITIONS=1L,OUTPUT_PRECISION=12L,THRESHOLD_MODELS=5000, stringsAsFactors = T)
      # print(configuration)
    }
  }
  else
  {
    configuration <- sRACIPE_load_configuration(config_file = config_file)
    message("Topology file loaded.")
    #print(configuration)
  }


  if(!missing(NUM_MODELS)){
    configuration$NUM_MODELS <- NUM_MODELS
  }
  if(!missing(PARAMETER_RANGE)){
    configuration$PARAMETER_RANGE <- PARAMETER_RANGE
  }
  if(!missing(MPR_MIN)){
    configuration$MPR_MIN <- MPR_MIN
  }
  if(!missing(MPR_MAX)){
    configuration$MPR_MAX <- MPR_MAX
  }

  if(!missing(DNR_MIN)){
    configuration$DNR_MIN <- DNR_MIN
  }
  if(!missing(DNR_MAX)){
    configuration$DNR_MAX <- DNR_MAX
  }
  if(!missing(FCH_MIN)){
    configuration$FCH_MIN <- FCH_MIN
  }
  if(!missing(FCH_MAX)){
    configuration$FCH_MAX <- FCH_MAX
  }
  if(!missing(HCO_MIN)){
    configuration$HCO_MIN <- HCO_MIN
  }
  if(!missing(HCO_MAX)){
    configuration$HCO_MAX <- HCO_MAX
  }
  if(!missing(STEP_SIZE)){
    configuration$STEP_SIZE <- STEP_SIZE
  }
  if(!missing(SIM_TIME)){
    configuration$SIM_TIME <- SIM_TIME
  }
  if(!missing(MEDIAN_RANGE)){
    configuration$MEDIAN_RANGE <- MEDIAN_RANGE
  }
  if(!missing(INITIAL_CONDITIONS)){
    configuration$INITIAL_CONDITIONS <- INITIAL_CONDITIONS
  }
  if(!missing(OUTPUT_PRECISION)){
    configuration$OUTPUT_PRECISION <- OUTPUT_PRECISION
  }
  if(!missing(THRESHOLD_MODELS)){
    configuration$THRESHOLD_MODELS <- THRESHOLD_MODELS
  }


  #print("Configuration file successfully loaded.")
  configuration$possible_interactions <- 3
  configuration$standard_deviation_factor <- 0.5658

  configuration$MPR_MIN <- 0.5*(configuration$MPR_MIN + configuration$MPR_MAX) - 0.5*(configuration$MPR_MAX - configuration$MPR_MIN)*configuration$PARAMETER_RANGE/100
  configuration$MPR_MAX <- 0.5*(configuration$MPR_MIN + configuration$MPR_MAX) + 0.5*(configuration$MPR_MAX - configuration$MPR_MIN)*configuration$PARAMETER_RANGE/100

  configuration$DNR_MIN <- 0.5*(configuration$DNR_MIN + configuration$DNR_MAX) - 0.5*(configuration$DNR_MAX - configuration$DNR_MIN)*configuration$PARAMETER_RANGE/100
  configuration$DNR_MAX <- 0.5*(configuration$DNR_MIN + configuration$DNR_MAX) + 0.5*(configuration$DNR_MAX - configuration$DNR_MIN)*configuration$PARAMETER_RANGE/100

  configuration$FCH_MIN <- 0.5*(configuration$FCH_MIN + configuration$FCH_MAX) - 0.5*(configuration$FCH_MAX - configuration$FCH_MIN)*configuration$PARAMETER_RANGE/100
  configuration$FCH_MAX <- 0.5*(configuration$FCH_MIN + configuration$FCH_MAX) + 0.5*(configuration$FCH_MAX - configuration$FCH_MIN)*configuration$PARAMETER_RANGE/100

  print(configuration)

  results_directory<-ifelse(!dir.exists(file.path(working_directory, "results")), dir.create(file.path(working_directory, "results")), TRUE)
  output_directory<-file.path(working_directory, "results")

  gene_interaction <- matrix(0, nrow = topology$number_gene, ncol = topology$number_gene)
  threshold_gene <- rep(0, topology$number_gene)
  #Rcpp::sourceCpp("src/interaction_reader.cpp")
  gene_interaction <- interaction_reader(gene_interaction,topology$topology_filepath,  topology$filename, topology$number_gene)
  knockout_number <- NA_integer_
  if(!missing(KNOCKOUT)){
    working_directory <- getwd()
    name_genes <- read.table(paste(working_directory,"/results/gene_interaction_topology_",topology$filename,".txt",sep=""), header = T, stringsAsFactors = F)
    name_genes <- t(as.matrix(name_genes))
    knockout_number <- as.integer(which(name_genes==KNOCKOUT))-1
    if(length(knockout_number)==0)
      {
      print("The knockout gene not found in the list of genes.")
      return()
    }
  }


print(knockout_number)
  message("Generating gene thresholds")
  #Rcpp::sourceCpp("src/threshold_generator.cpp")

  threshold_test<- threshold_calculator_uniform( gene_interaction,  threshold_gene,  configuration$MPR_MIN,  configuration$MPR_MAX, configuration$DNR_MIN,  configuration$DNR_MAX,  configuration$possible_interactions,   configuration$NUM_MODELS,  configuration$THRESHOLD_MODELS,  configuration$STEP_SIZE,  configuration$FCH_MIN, configuration$FCH_MAX,  configuration$HCO_MIN,  configuration$HCO_MAX,  configuration$MEDIAN_RANGE,  configuration$standard_deviation_factor)
  if(threshold_test!=0)
    message("Error in threshold generation")
  #print(topology$number_gene)
  #return()
  message("Running the determinsitic simulations using Runge-Kutta")

  #Rcpp::sourceCpp("src/multiGeneCircuit_RK_deterministic.cpp")

  #print(topology$number_gene)
out_file <- paste(topology$filename,KNOCKOUT,"KO",sep = "_")
  # else
  {
  Time_evolution_test<- multiGeneCircuit_RK_deterministic_knockout( gene_interaction,  threshold_gene,  configuration$MPR_MIN,  configuration$MPR_MAX, configuration$DNR_MIN,  configuration$DNR_MAX,configuration$possible_interactions,configuration$NUM_MODELS, configuration$THRESHOLD_MODELS,  configuration$STEP_SIZE,  configuration$FCH_MIN,configuration$FCH_MAX, configuration$HCO_MIN,  configuration$HCO_MAX,  configuration$SIM_TIME, configuration$MEDIAN_RANGE, configuration$standard_deviation_factor, topology$number_gene,  configuration$OUTPUT_PRECISION, configuration$INITIAL_CONDITIONS, out_file, knockout_number)
}

  if(Time_evolution_test!=0){
    message("Error in time evolution function")
  }

  output_file <- paste(working_directory,"/results/sRACIPE_RK_",out_file,"_g",toString(topology$number_gene), "_output.txt", sep="")
  return(output_file)


}
